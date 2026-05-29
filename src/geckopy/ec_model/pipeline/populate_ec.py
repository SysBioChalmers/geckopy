"""Build and populate the ec substructure on an EcModel.

Covers stages 6, 7, and 8 of GECKO MATLAB `makeEcModel.m`:

- Stage 6 allocates empty per-reaction slots in ec (one per catalyzed
  reaction after stage-5 expansion).
- Stage 7 fills in per-enzyme data (genes, uniprot IDs, MW, sequence,
  concentrations) by looking up each model gene in the UniProt
  database. Genes that do not match UniProt are left out of ec and
  returned as a warning list; their parent reactions are also flagged
  with a note in rxn.notes.
- Stage 8 populates the reaction-to-enzyme coupling matrix
  (ec.rxn_enz_mat) with 1.0 for each (reaction, enzyme) pair implied
  by the reaction's GPR.

These three stages are only meaningful for the full (non-light) ecModel
format. Gecko-light is handled by a separate future module.

All three functions mutate the ec substructure in place. None of them
talk to the network.
"""
from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import numpy as np
from scipy import sparse

if TYPE_CHECKING:
    from ...databases import UniprotDB
    from ...databases.kegg_loader import KeggDB
    from ..ec_model import EcModel

logger = logging.getLogger(__name__)


# --------------------------------------------------------------------------- #
# Stage 6: allocate empty ec slots
# --------------------------------------------------------------------------- #

def allocate_ec_for_catalyzed_reactions(model: "EcModel") -> list[str]:
    """Stage 6: allocate ec.rxns and per-reaction slots.

    Ported from GECKO MATLAB: src/geckomat/change_model/makeEcModel.m
    (stage 6, full-model branch only). Gecko-light uses a different
    allocation scheme and is not yet supported.

    Walks the model (which must already have been through stages 1-5)
    and selects every reaction that has at least one gene associated.
    Each selected reaction gets an entry in ec.rxns; ec.kcat is
    initialized to NaN, and ec.source/notes/eccodes to empty strings,
    ec.concs to NaN. The per-enzyme arrays are left untouched (stage 7
    populates them).

    Parameters
    ----------
    model
        An EcModel, already preprocessed through stage 5. Mutated in place.

    Returns
    -------
    list of str
        The reaction IDs written into ec.rxns, in model order.
    """
    rxn_ids: list[str] = [r.id for r in model.reactions if r.genes]
    n = len(rxn_ids)

    model.ec.rxns = rxn_ids
    # 0 marks "no kcat assigned" (matching MATLAB GECKO).
    model.ec.kcat = np.zeros(n, dtype=float)
    model.ec.source = [""] * n
    model.ec.notes = [""] * n
    model.ec.eccodes = [""] * n
    # ec.concs is a per-enzyme field, not per-reaction. Stage 7 sizes it.

    return rxn_ids


# --------------------------------------------------------------------------- #
# Stage 7: populate per-enzyme data from UniProt
# --------------------------------------------------------------------------- #

_UNCONSTRAINED_NOTE_KEY = "geckopy_warning"


def populate_enzyme_data(
    model: "EcModel",
    uniprot_db: "UniprotDB",
    *,
    kegg_db: "KeggDB | None" = None,
) -> list[str]:
    """Stage 7: fill per-enzyme fields from UniProt (with optional KEGG fallback).

    Ported from GECKO MATLAB: src/geckomat/change_model/makeEcModel.m
    (stage 7), extended with a KEGG fallback that MATLAB does not have.

    MATLAB-COMPAT: GECKO MATLAB iterates model.genes in their original
    model order. geckopy iterates them in alphabetical order so that
    ec.genes / ec.enzymes / ec.mw / ec.sequence are deterministic
    regardless of how the SBML was loaded.

    For every gene in the model (sorted alphabetically for reproducibility),
    a lookup is performed:

    1. ``adapter.get_uniprot_compatible_genes`` transforms model gene IDs
       into a form that should match the UniProt "Gene Names" column.
    2. If an ``uniprotConversion.tsv`` is present, it instead transforms
       model gene IDs directly to UniProt IDs and matches on the "Entry"
       column.

    Matched genes are added to ec.genes / ec.enzymes / ec.mw /
    ec.sequence (same length as ec.genes), and ec.concs is allocated as
    NaN of matching length. Unmatched genes are returned as a list
    (the equivalent of MATLAB's ``noUniprot``).

    When ``kegg_db`` is supplied, each gene UniProt could not match is
    looked up in the KEGG database via
    ``adapter.get_kegg_compatible_genes``. KEGG hits feed ec.enzymes
    with the UniProt accession stored on the KEGG entry; if that is
    empty, the bare KEGG gene ID (``KeggDB.kegg_genes``) is used as a
    fallback identifier. The returned ``no_uniprot`` list then names
    only the genes that neither source could fill.

    MATLAB-COMPAT: GECKO MATLAB only returns the noUniprot list and
    does not annotate the affected reactions. geckopy additionally
    writes a note to rxn.notes['geckopy_warning'] for each affected
    reaction so users can see the issue when inspecting individual
    reactions. MATLAB GECKO could be updated to do the same.

    Reactions whose GPR references an unmatched gene get a warning note
    added to ``rxn.notes[geckopy_warning]``. If such a reaction has at
    least one other matched gene, it will still be enzyme-constrained
    through that gene. If all its genes are unmatched it will be
    unconstrained.

    Parameters
    ----------
    model
        An EcModel. Must have ``adapter`` set. Mutated in place.
    uniprot_db
        Loaded UniProt database.

    Returns
    -------
    list of str
        Model gene IDs (pre-transformation) that did not match.
    """
    from ...adapter import resolve_adapter
    adapter = resolve_adapter(
        model,
        purpose="populate_enzyme_data transforms gene IDs to "
        "UniProt-compatible form via the adapter",
    )

    # Sorted alphabetically so downstream indexing is deterministic.
    model_genes: list[str] = sorted(g.id for g in model.genes)
    if not model_genes:
        model.ec.genes = []
        model.ec.enzymes = []
        model.ec.mw = np.empty(0, dtype=float)
        model.ec.sequence = []
        model.ec.concs = np.empty(0, dtype=float)
        return []

    # Step 1: try gene-based lookup via adapter.get_uniprot_compatible_genes
    transformed = adapter.get_uniprot_compatible_genes(model_genes)
    matches: list[int | None] = [
        uniprot_db.find_by_gene(g) for g in transformed
    ]

    # Step 2: if the adapter has a conversion table, that takes over.
    mapped_ids = adapter.get_uniprot_ids_from_table(transformed)
    if mapped_ids != transformed:
        matches = [
            uniprot_db.find_by_id(i) if i else None for i in mapped_ids
        ]

    matched_genes: list[str] = []
    matched_enzymes: list[str] = []
    matched_mw: list[float] = []
    matched_seq: list[str] = []
    unmatched: list[str] = []
    kegg_filled: list[str] = []
    kegg_used_bare_id: list[tuple[str, str]] = []

    kegg_lookup: dict[str, int] = {}
    kegg_transformed: list[str] = []
    if kegg_db is not None and len(kegg_db) > 0:
        kegg_transformed = adapter.get_kegg_compatible_genes(model_genes)
        for i, g in enumerate(kegg_db.genes):
            if g and g not in kegg_lookup:
                kegg_lookup[g] = i

    for idx, (gene, row_idx) in enumerate(zip(model_genes, matches)):
        if row_idx is not None:
            matched_genes.append(gene)
            matched_enzymes.append(uniprot_db.ids[row_idx])
            matched_mw.append(uniprot_db.mw[row_idx])
            matched_seq.append(uniprot_db.sequences[row_idx])
            continue
        kegg_idx: int | None = None
        if kegg_lookup:
            kegg_idx = kegg_lookup.get(kegg_transformed[idx])
        if kegg_idx is None:
            unmatched.append(gene)
            continue
        # KEGG fallback hit. Prefer the UniProt accession carried on the
        # KEGG row; fall back to the bare KEGG gene id when empty.
        kegg_uniprot = kegg_db.uniprot_ids[kegg_idx]
        if kegg_uniprot:
            enzyme_id = kegg_uniprot
        else:
            enzyme_id = kegg_db.kegg_genes[kegg_idx]
            kegg_used_bare_id.append((gene, enzyme_id))
        matched_genes.append(gene)
        matched_enzymes.append(enzyme_id)
        matched_mw.append(float(kegg_db.mw[kegg_idx]))
        matched_seq.append(kegg_db.sequences[kegg_idx])
        kegg_filled.append(gene)

    if not matched_genes:
        raise ValueError(
            "None of the model genes matched an entry in the UniProt "
            "database (or KEGG, if provided). Check "
            "adapter.params.uniprot.* / adapter.params.kegg.* and, "
            "if needed, provide a data/uniprotConversion.tsv file."
        )

    model.ec.genes = matched_genes
    model.ec.enzymes = matched_enzymes
    model.ec.mw = np.array(matched_mw, dtype=float)
    model.ec.sequence = matched_seq
    model.ec.concs = np.full(len(matched_genes), np.nan, dtype=float)

    # Per-enzyme rows are keyed by gene, so two genes that map to the same
    # UniProt accession produce duplicate ec.enzymes entries. apply_kcat now
    # sums their coupling, but the duplication is worth surfacing because the
    # per-enzyme mw/concs of the second row are otherwise easy to overlook.
    seen: set[str] = set()
    duplicates = sorted({e for e in matched_enzymes if e in seen or seen.add(e)})
    if duplicates:
        logger.warning(
            "populate_enzyme_data: %d UniProt accession(s) are shared by "
            "multiple genes (e.g. %s); their enzyme rows are duplicated in "
            "model.ec.enzymes.",
            len(duplicates), ", ".join(duplicates[:5]),
        )

    if kegg_filled:
        logger.info(
            "populate_enzyme_data: %d gene(s) resolved via KEGG fallback "
            "(UniProt had no entry): %s",
            len(kegg_filled),
            ", ".join(kegg_filled[:10]) + ("..." if len(kegg_filled) > 10 else ""),
        )
        if kegg_used_bare_id:
            preview = ", ".join(
                f"{g}->{eid}" for g, eid in kegg_used_bare_id[:10]
            )
            more = (
                "" if len(kegg_used_bare_id) <= 10
                else f" (and {len(kegg_used_bare_id) - 10} more)"
            )
            logger.warning(
                "populate_enzyme_data: %d KEGG-resolved gene(s) had no "
                "UniProt accession on the KEGG row; the bare KEGG gene id "
                "is used in ec.enzymes instead: %s%s",
                len(kegg_used_bare_id), preview, more,
            )

    # Annotate reactions whose GPR mentions any unmatched gene.
    if unmatched:
        unmatched_set = set(unmatched)
        for rxn in model.reactions:
            if not rxn.genes:
                continue
            missing = sorted(
                g.id for g in rxn.genes if g.id in unmatched_set
            )
            if missing:
                rxn.notes[_UNCONSTRAINED_NOTE_KEY] = (
                    "geckopy: no ec-constraint due to no uniprot/kegg match "
                    "for gene(s): " + ", ".join(missing)
                )

    return unmatched


# --------------------------------------------------------------------------- #
# Stage 8: build ec.rxn_enz_mat
# --------------------------------------------------------------------------- #

def build_rxn_enzyme_coupling(model: "EcModel") -> None:
    """Stage 8: populate ec.rxn_enz_mat.

    Ported from GECKO MATLAB: src/geckomat/change_model/makeEcModel.m
    (stage 8, full-model branch only).

    For each reaction ID in ec.rxns and each gene in ec.genes, writes
    1.0 at (i, j) if the reaction's GPR contains that gene. After
    stage 5, every reaction has at most one AND-clause, so this is
    simply "which enzymes are subunits of the complex catalyzing this
    reaction."

    Downstream functions (e.g. applyComplexData) may later overwrite
    these 1.0 entries with actual subunit stoichiometries.

    Parameters
    ----------
    model
        An EcModel with ec.rxns and ec.genes already populated
        (stages 6 and 7). Mutated in place.
    """
    n_rxns = model.ec.n_rxns
    n_enz = model.ec.n_enzymes

    if n_rxns == 0 or n_enz == 0:
        model.ec.rxn_enz_mat = sparse.csr_matrix((n_rxns, n_enz), dtype=float)
        return

    enzyme_index: dict[str, int] = {g: i for i, g in enumerate(model.ec.genes)}

    mat = sparse.lil_matrix((n_rxns, n_enz), dtype=float)
    for i, rxn_id in enumerate(model.ec.rxns):
        rxn = model.reactions.get_by_id(rxn_id)
        for gene in rxn.genes:
            j = enzyme_index.get(gene.id)
            if j is not None:
                mat[i, j] = 1.0

    model.ec.rxn_enz_mat = mat.tocsr()
