"""SBML reader and writer for ecModels.

SBML is the standard exchange format for constraint-based
modelling. Writing an ecModel as SBML means other tools (MATLAB
GECKO, Memote, escher, BiGG, ...) can read at least the cobra
portion. Tools that don't understand the GECKO extensions still
see a self-consistent FBA model — the enzyme mass-balance
constraints are explicit in the stoichiometry.

The encoding (the "MW_KCAT" format, GECKO 3's default) is:

- Each enzyme is a regular metabolite named ``prot_<uniprot>``.
- These metabolites are grouped in an SBML ``Group`` named
  ``Protein`` so an SBML-aware reader can find them.
- Molecular weight is carried in the species notes (cobra's
  ``metabolite.notes`` dict).
- Proteomics concentration goes in the species
  ``initialAmount``.
- The stoichiometric coefficient of ``prot_<id>`` on a catalysed
  reaction is ``-MW / (3600 * kcat)`` per subunit, matching the
  GECKO 3 MATLAB convention exactly so no rescaling is needed at
  round-trip.

Ported from the legacy geckopy package described in Carrasco et al.
(2023, https://doi.org/10.1128/spectrum.01705-23), file
geckopy/io/sbml.py. Simplified to a single encoding (no
``EcStoichiometry`` enum, no dual KCAT mode).
"""
from __future__ import annotations

import logging
import re
from pathlib import Path
from typing import TYPE_CHECKING, Optional

import cobra
import libsbml
import numpy as np
from scipy import sparse

from ..ec_model.constants import POOL_EXCHANGE_ID, POOL_ID, PROT_PREFIX, USAGE_PREFIX

if TYPE_CHECKING:
    from ..adapter.adapter import ModelAdapter
    from ..ec_model.ec_model import EcModel

logger = logging.getLogger(__name__)

_PROTEIN_GROUP = "Protein"
_COBRA_MET_PREFIX = "M_"


# --------------------------------------------------------------------------- #
# Write
# --------------------------------------------------------------------------- #

def write_sbml_ec_model(model: "EcModel", filename: str | Path) -> None:
    """Write an EcModel as SBML with enzymes in a ``Protein`` Group.

    Coefficients on ``prot_<id>`` metabolites are scaled to
    ``-MW/(3600*kcat)`` so tools reading the file without geckopy
    get a self-consistent FBA model where enzyme balances are
    explicit. The Protein group + per-species MW notes are the
    GECKO-specific extensions and are silently ignored by such
    tools.

    Ported from the legacy geckopy package (Carrasco et al., 2023,
    https://doi.org/10.1128/spectrum.01705-23),
    geckopy/io/sbml.py:799-1190
    (write_sbml_ec_model), simplified to MW_KCAT encoding only.
    """
    filename = str(filename)

    # 1. Tag prot_<id> metabolites with MW/sequence/gene and catalysed
    #    reactions with their kcat source/eccodes/notes/kcat/subunit
    #    counts (cobra serialises notes dicts into SBML <notes>; this
    #    carries through the round-trip cleanly without fighting
    #    libsbml.appendNotes). The FBA coefficient already encodes
    #    subunits*MW/(kcat*3600), but those three are not separable from
    #    the coefficient alone, so they are stored explicitly.
    annotated = _annotate_ec_metadata(model)

    # 2. Write via cobra (handles all standard SBML). The in-memory
    #    coefficient on prot_<id> mets is already MW_KCAT-encoded
    #    (-subunits * MW / (kcat * 3600)), matching the on-disk
    #    convention used by MATLAB GECKO, so no rescaling needed.
    cobra.io.write_sbml_model(annotated, filename)

    # 3. Inject Protein group + MW notes via libsbml.
    reader = libsbml.SBMLReader()
    doc = reader.readSBML(filename)
    sbml_model = doc.getModel()

    # Enable Groups extension (cobra usually loads it).
    if not doc.isPackageEnabled("groups"):
        doc.enablePackage(
            "http://www.sbml.org/sbml/level3/version1/groups/version1",
            "groups",
            True,
        )
    groups_plugin = sbml_model.getPlugin("groups")
    if groups_plugin is None:
        raise RuntimeError(
            "libsbml groups package not available; cannot write Protein group"
        )

    group = groups_plugin.createGroup()
    group.setId(_PROTEIN_GROUP)
    group.setName(_PROTEIN_GROUP)
    group.setKind("collection")

    for u in model.ec.enzymes:
        species_id = f"{_COBRA_MET_PREFIX}{PROT_PREFIX}{u}"
        species = sbml_model.getSpecies(species_id)
        if species is None:
            continue
        idx = model.ec.enzymes.index(u)
        conc = float(model.ec.concs[idx])
        species.setInitialAmount(0.0 if np.isnan(conc) else conc)
        # MW is written into the species notes via cobra's metabolite.notes
        # dict in _scale_for_sbml (avoids libsbml.appendNotes returning -5).
        member = group.createMember()
        member.setIdRef(species_id)

    libsbml.writeSBMLToFile(doc, filename)


def _annotate_ec_metadata(model: "EcModel") -> "EcModel":
    """Return a deep-copy of ``model`` with the ec provenance that is
    not recoverable from the FBA stoichiometry written into notes.

    cobra serialises ``notes`` dicts into SBML ``<notes>`` as
    ``<p>key: value</p>`` entries, which both cobra (on read) and
    geckopy's reader pick up again. This avoids the libsbml
    ``appendNotes`` quirk (returns status -5 even for well-formed
    XHTML).

    Written here:

    - Per ``prot_<id>`` metabolite: ``mw``, ``sequence``, ``gene``.
    - Per catalysed reaction (matched by id to ``ec.rxns``):
      ``ec_source``, ``ec_eccodes``, ``ec_notes``, ``ec_kcat`` and
      ``ec_subunits`` (a ``<enzyme>:<count>`` list). The coefficient
      on ``prot_<id>`` is ``-subunits*MW/(3600*kcat)``; storing kcat
      and subunit counts lets the reader recover both exactly instead
      of assuming one subunit and dividing the kcat.
    """
    # All values written below are read from the original ``model.ec``;
    # ``out.ec`` is never touched, so there is no need to deep-copy it.
    out = model.copy()

    out_met_ids = {m.id for m in out.metabolites}
    out_rxn_ids = {r.id for r in out.reactions}

    # Per-enzyme: MW, sequence, gene.
    for i, enz in enumerate(model.ec.enzymes):
        met_id = f"{PROT_PREFIX}{enz}"
        if met_id not in out_met_ids:
            continue
        met = out.metabolites.get_by_id(met_id)
        mw = float(model.ec.mw[i])
        if np.isfinite(mw):
            met.notes["mw"] = str(mw)
        if i < len(model.ec.sequence) and model.ec.sequence[i]:
            met.notes["sequence"] = model.ec.sequence[i]
        if i < len(model.ec.genes):
            gene = model.ec.genes[i]
            # Only worth storing when it differs from the accession;
            # the reader defaults gene to the accession otherwise.
            if gene and gene != enz:
                met.notes["gene"] = gene

    # Per catalysed reaction: source, eccodes, notes, kcat, subunits.
    mat = model.ec.rxn_enz_mat.tocsr()
    enzymes = model.ec.enzymes
    for i, rid in enumerate(model.ec.rxns):
        if rid not in out_rxn_ids:
            continue
        rxn = out.reactions.get_by_id(rid)
        if i < len(model.ec.source) and model.ec.source[i]:
            rxn.notes["ec_source"] = model.ec.source[i]
        if i < len(model.ec.eccodes) and model.ec.eccodes[i]:
            rxn.notes["ec_eccodes"] = model.ec.eccodes[i]
        if i < len(model.ec.notes) and model.ec.notes[i]:
            rxn.notes["ec_notes"] = model.ec.notes[i]
        kcat = float(model.ec.kcat[i]) if i < len(model.ec.kcat) else 0.0
        if np.isfinite(kcat) and kcat > 0:
            rxn.notes["ec_kcat"] = repr(kcat)
        if i < mat.shape[0]:
            row = mat.getrow(i)
            if row.nnz:
                parts = [f"{enzymes[j]}:{row[0, j]:g}" for j in row.indices]
                rxn.notes["ec_subunits"] = ";".join(parts)
    return out


# --------------------------------------------------------------------------- #
# Read
# --------------------------------------------------------------------------- #

def read_sbml_ec_model(
    filename: str | Path,
    *,
    adapter: Optional["ModelAdapter"] = None,
) -> "EcModel":
    """Read an SBML file written with MW_KCAT encoding back into an EcModel.

    Reconstructs ``model.ec`` from Protein group membership, species
    notes (MW), and reaction stoichiometry. Inverse of
    ``write_sbml_ec_model``; pure-FBA SBMLs without the Protein
    group come back as an EcModel with an empty ``ec`` substructure.

    Parameters
    ----------
    filename
        Path to the SBML file.
    adapter
        ModelAdapter to attach to the returned EcModel, or ``None`` to
        load the model for inspection without a project. Downstream
        functions that need organism parameters will then require an
        explicit ``adapter``.

    Ported from the legacy geckopy package (Carrasco et al., 2023,
    https://doi.org/10.1128/spectrum.01705-23),
    geckopy/io/sbml.py:95-700
    (read_sbml_ec_model), simplified to MW_KCAT only.
    """
    from ..ec_model.ec_data import EcData
    from ..ec_model.ec_model import EcModel

    filename = str(filename)

    # 1. Standard cobra read.
    cobra_model = cobra.io.read_sbml_model(filename)

    # 2. Walk libsbml document for the Protein group.
    reader = libsbml.SBMLReader()
    doc = reader.readSBML(filename)
    sbml_model = doc.getModel()
    groups_plugin = sbml_model.getPlugin("groups")

    enzyme_ids: list[str] = []
    enzyme_mw: dict[str, float] = {}
    enzyme_conc: dict[str, float] = {}
    enzyme_seq: dict[str, str] = {}
    enzyme_gene: dict[str, str] = {}
    cobra_met_ids = {m.id for m in cobra_model.metabolites}
    if groups_plugin is not None:
        for gi in range(groups_plugin.getNumGroups()):
            g = groups_plugin.getGroup(gi)
            if g.getId() != _PROTEIN_GROUP:
                continue
            for mi in range(g.getNumMembers()):
                member = g.getMember(mi)
                sid = member.getIdRef()
                species = sbml_model.getSpecies(sid)
                if species is None:
                    continue
                clean = (
                    sid[len(_COBRA_MET_PREFIX):]
                    if sid.startswith(_COBRA_MET_PREFIX) else sid
                )
                if not clean.startswith(PROT_PREFIX):
                    continue
                uniprot = clean[len(PROT_PREFIX):]
                if uniprot == POOL_ID[len(PROT_PREFIX):]:
                    continue
                enzyme_ids.append(uniprot)
                # Prefer cobra's parsed notes (reliable); fall back to
                # raw libsbml notes (works for hand-written SBMLs that
                # don't follow cobra's `<p>key: value</p>` convention).
                mw = None
                if clean in cobra_met_ids:
                    notes = cobra_model.metabolites.get_by_id(clean).notes
                    if "mw" in notes:
                        try:
                            mw = float(notes["mw"])
                        except (TypeError, ValueError):
                            mw = None
                    if notes.get("sequence"):
                        enzyme_seq[uniprot] = notes["sequence"]
                    if notes.get("gene"):
                        enzyme_gene[uniprot] = notes["gene"]
                if mw is None:
                    mw = _parse_mw_from_notes(species.getNotesString())
                if mw is not None:
                    enzyme_mw[uniprot] = mw
                amt = species.getInitialAmount()
                enzyme_conc[uniprot] = float(amt) if amt > 0 else float("nan")

    # 3. Wrap and populate.
    ec_model = EcModel.from_cobra(cobra_model, adapter)
    if not enzyme_ids:
        ec_model.ec = EcData()
        return ec_model

    _populate_ec_from_sbml(
        ec_model, enzyme_ids, enzyme_mw, enzyme_conc,
        enzyme_seq, enzyme_gene,
    )
    return ec_model


_MW_NOTE_RE = re.compile(r"mw:\s*([0-9.eE+\-]+)")


def _parse_mw_from_notes(notes_string: str) -> Optional[float]:
    if not notes_string:
        return None
    m = _MW_NOTE_RE.search(notes_string)
    return float(m.group(1)) if m else None


def _parse_subunits(value: str) -> dict[str, float]:
    """Parse an ``ec_subunits`` note (``<enzyme>:<count>;...``)."""
    counts: dict[str, float] = {}
    for part in value.split(";"):
        if ":" not in part:
            continue
        enz, cnt = part.rsplit(":", 1)
        try:
            counts[enz] = float(cnt)
        except ValueError:
            continue
    return counts


def _populate_ec_from_sbml(
    ec_model: "EcModel",
    enzyme_ids: list[str],
    enzyme_mw: dict[str, float],
    enzyme_conc: dict[str, float],
    enzyme_seq: Optional[dict[str, str]] = None,
    enzyme_gene: Optional[dict[str, str]] = None,
) -> None:
    """Reconstruct ``model.ec`` from per-enzyme metadata + reaction
    stoichiometry.

    For each metabolic reaction that has a ``prot_<id>`` in its
    stoichiometry, the kcat, subunit count, source, eccodes and notes
    are taken from the reaction's ``ec_*`` notes when present (written
    by ``write_sbml_ec_model``). For files that lack those notes
    (legacy/foreign SBML), kcat falls back to inverting the coefficient
    ``kcat = -mw / (coef * 3600)`` with one subunit assumed.

    Ported in concept from the legacy geckopy package
    (Carrasco et al., 2023, https://doi.org/10.1128/spectrum.01705-23),
    geckopy/io/sbml.py:540-787.
    """
    from ..ec_model.ec_data import EcData

    enzyme_seq = enzyme_seq or {}
    enzyme_gene = enzyme_gene or {}

    # Per-enzyme arrays (order preserved from group iteration).
    n_e = len(enzyme_ids)
    mw_arr = np.array(
        [enzyme_mw.get(u, np.nan) for u in enzyme_ids], dtype=float,
    )
    concs_arr = np.array(
        [enzyme_conc.get(u, np.nan) for u in enzyme_ids], dtype=float,
    )
    sequence = [enzyme_seq.get(u, "") for u in enzyme_ids]
    # gene defaults to the accession when not stored separately.
    genes = [enzyme_gene.get(u, u) for u in enzyme_ids]
    enzyme_index = {u: i for i, u in enumerate(enzyme_ids)}

    # Walk metabolic reactions; collect (rxn_idx, enz_idx, kcat, subunits).
    catalysed_rxns: list[str] = []
    rxn_index: dict[str, int] = {}
    quads: list[tuple[int, int, float, float]] = []
    src_by_idx: dict[int, str] = {}
    ecc_by_idx: dict[int, str] = {}
    notes_by_idx: dict[int, str] = {}

    for rxn in ec_model.reactions:
        if rxn.id == POOL_EXCHANGE_ID or rxn.id.startswith(USAGE_PREFIX):
            continue
        note_kcat: Optional[float] = None
        if "ec_kcat" in rxn.notes:
            try:
                note_kcat = float(rxn.notes["ec_kcat"])
            except (TypeError, ValueError):
                note_kcat = None
        subunit_counts = (
            _parse_subunits(rxn.notes["ec_subunits"])
            if "ec_subunits" in rxn.notes else {}
        )
        for met, coef in rxn.metabolites.items():
            if met.id == POOL_ID or not met.id.startswith(PROT_PREFIX):
                continue
            enzyme = met.id[len(PROT_PREFIX):]
            j = enzyme_index.get(enzyme)
            if j is None:
                continue
            mw = mw_arr[j]
            if not np.isfinite(mw) or coef == 0:
                continue
            if note_kcat is not None:
                kcat = note_kcat
                subunits = subunit_counts.get(enzyme, 1.0)
            else:
                # coef = -mw / (3600 * kcat)  =>  kcat = -mw / (3600 * coef)
                kcat = float(-mw / (3600.0 * coef))
                subunits = 1.0
            if rxn.id not in rxn_index:
                idx = len(catalysed_rxns)
                rxn_index[rxn.id] = idx
                catalysed_rxns.append(rxn.id)
                src_by_idx[idx] = rxn.notes.get("ec_source", "")
                ecc_by_idx[idx] = rxn.notes.get("ec_eccodes", "")
                notes_by_idx[idx] = rxn.notes.get("ec_notes", "")
            quads.append((rxn_index[rxn.id], j, kcat, subunits))

    n_r = len(catalysed_rxns)
    # 0 marks "no kcat assigned".
    kcat_arr = np.zeros(n_r, dtype=float)
    mat = sparse.lil_matrix((n_r, n_e), dtype=float)
    for rxn_idx, enz_idx, kcat, subunits in quads:
        # kcat is per reaction (one value per ec.rxns row); subunits are
        # per (reaction, enzyme).
        kcat_arr[rxn_idx] = kcat
        mat[rxn_idx, enz_idx] = subunits

    ec_model.ec = EcData(
        rxns=catalysed_rxns,
        kcat=kcat_arr,
        source=[src_by_idx.get(i, "") for i in range(n_r)],
        notes=[notes_by_idx.get(i, "") for i in range(n_r)],
        eccodes=[ecc_by_idx.get(i, "") for i in range(n_r)],
        genes=genes,
        enzymes=list(enzyme_ids),
        mw=mw_arr,
        sequence=sequence,
        concs=concs_arr,
        rxn_enz_mat=mat.tocsr(),
    )
