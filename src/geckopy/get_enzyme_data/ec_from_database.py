"""Populate model.ec.eccodes by looking up each reaction's genes in a
protein database.

Ported from GECKO MATLAB:
src/geckomat/get_enzyme_data/getECfromDatabase.m.
"""
from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import TYPE_CHECKING, Iterable, Literal, Optional

from .find_ec_in_db import find_ec_in_db

if TYPE_CHECKING:
    from ..adapter import ModelAdapter
    from ..databases import UniprotDB
    from ..ec_model.ec_model import EcModel

logger = logging.getLogger(__name__)


@dataclass
class _Conflict:
    """One gene with multiple distinct EC strings in the DB, in the
    context of one reaction. Internal aggregation for the 'display'
    action."""
    rxn_idx: int
    gene: str
    protein_indices: list[int]
    db_name: str = "uniprot"


def fill_eccodes_from_database(
    model: "EcModel",
    uniprot_db: "UniprotDB",
    *,
    ec_rxns: Optional[Iterable[str]] = None,
    action: Literal["display", "ignore"] = "display",
) -> None:
    """Populate ``model.ec.eccodes`` by looking each reaction's genes
    up in a UniProt database.

    Ported from GECKO MATLAB:
    src/geckomat/get_enzyme_data/getECfromDatabase.m.

    For each entry of ``model.ec.rxns`` selected by ``ec_rxns``:

    * Gather the gene IDs catalysing the reaction from
      ``model.ec.rxn_enz_mat`` row ``i``.
    * Pass them to ``find_ec_in_db`` along with a pre-built
      gene -> protein-indices map.
    * Store the resulting ``;``-joined EC string in
      ``model.ec.eccodes[i]``.

    The gene -> protein-indices map is built once before the loop. It
    tries the adapter's ``uniprotConversion.tsv`` table first: if that
    table changes any gene IDs, each model gene maps to one UniProt ID
    which uniquely identifies one protein. Otherwise it falls back to
    name-based lookup via ``adapter.get_uniprot_compatible_genes``,
    and a single gene may map to several proteins (which is how
    conflicts arise).

    MATLAB-COMPAT: GECKO MATLAB takes a ``modelAdapter`` arg and
    resolves it via ``ModelAdapterManager.getDefault()``. geckopy
    reads the adapter from ``model.adapter``.

    MATLAB-COMPAT: GECKO MATLAB falls back to a KEGG lookup when the
    UniProt match is empty or has a trailing ``-``. geckopy does not
    yet implement the KEGG fallback (no ``KeggDB`` ported); the
    UniProt result is used as-is.

    MATLAB-COMPAT: GECKO MATLAB exposes three actions: ``'display'``
    (raise), ``'ignore'`` (silent), ``'add'`` (commented-out dead
    code). geckopy drops ``'add'``. ``'display'`` here emits a single
    aggregated ``logger.warning`` rather than raising, since the
    individual conflict warnings from ``find_ec_in_db`` already give
    users enough to act on. Tracked in
    ``docs/future_improvements.md``.

    MATLAB-COMPAT: GECKO MATLAB loops over every reaction even when
    ``ecRxns`` masks most of them, then subsets at the end. geckopy
    subsets up front (matching MATLAB's own TODO comment).

    Parameters
    ----------
    model
        EcModel with ``model.ec.rxns``, ``model.ec.genes``, and
        ``model.ec.rxn_enz_mat`` already populated (typically by
        ``allocate_ec_for_catalyzed_reactions``,
        ``populate_enzyme_data``, and ``build_rxn_enzyme_coupling``
        respectively). Mutated in place.
    uniprot_db
        Pre-loaded UniProt database. Accessed for its ``eccodes``,
        ``mw``, ``genes``, and ``ids`` fields.
    ec_rxns
        Optional iterable of reaction IDs (each must appear in
        ``model.ec.rxns``) selecting which entries of
        ``model.ec.eccodes`` to update. ``None`` means update all.
        Unknown IDs raise ``ValueError``.
    action
        Behaviour when a gene maps to proteins with multiple distinct
        EC strings. ``'display'`` (default) emits an aggregated
        warning at the end of the run listing every conflicting
        ``(rxn, gene, proteins)`` triple. ``'ignore'`` suppresses the
        aggregated warning. Per-gene warnings from ``find_ec_in_db``
        are emitted in both modes.

    Raises
    ------
    ValueError
        If ``model.adapter`` is None; if ``action`` is not
        ``'display'`` or ``'ignore'``; or if any ID in ``ec_rxns`` is
        not present in ``model.ec.rxns``.
    """
    if action not in ("display", "ignore"):
        raise ValueError(
            f"action must be 'display' or 'ignore', got {action!r}"
        )

    if model.adapter is None:
        raise ValueError(
            "EcModel.adapter is None; fill_eccodes_from_database needs an "
            "adapter to transform gene IDs."
        )

    n = model.ec.n_rxns
    if n == 0:
        return

    if ec_rxns is None:
        positions: list[int] = list(range(n))
    else:
        ec_rxns_list = list(ec_rxns)
        index_by_id = {rid: i for i, rid in enumerate(model.ec.rxns)}
        unknown = [rid for rid in ec_rxns_list if rid not in index_by_id]
        if unknown:
            preview = unknown[:5]
            raise ValueError(
                f"{len(unknown)} reaction ID(s) in ec_rxns are not "
                f"present in model.ec.rxns (examples: {preview})"
            )
        positions = [index_by_id[rid] for rid in ec_rxns_list]

    if not positions:
        return

    gene_to_protein_indices = _build_gene_to_protein_indices(
        model.ec.genes, uniprot_db, model.adapter,
    )

    db_eccodes = uniprot_db.eccodes
    db_mw = uniprot_db.mw

    new_eccodes = list(model.ec.eccodes)
    conflicts: list[_Conflict] = []
    rxn_enz_mat = model.ec.rxn_enz_mat.tocsr()

    for i in positions:
        gene_indices = rxn_enz_mat[i].indices
        if len(gene_indices) == 0:
            new_eccodes[i] = ""
            continue
        gene_set = [model.ec.genes[j] for j in gene_indices]

        rxn_conflicts: list[tuple[str, list[int]]] = []
        new_eccodes[i] = find_ec_in_db(
            gene_set, db_eccodes, db_mw, gene_to_protein_indices,
            conflicts=rxn_conflicts,
        )

        for gene, protein_indices in rxn_conflicts:
            conflicts.append(_Conflict(
                rxn_idx=i,
                gene=gene,
                protein_indices=protein_indices,
            ))

    model.ec.eccodes = new_eccodes

    if action == "display" and conflicts:
        logger.warning(_format_conflict_message(conflicts, uniprot_db, model))


def _build_gene_to_protein_indices(
    model_genes: list[str],
    uniprot_db: "UniprotDB",
    adapter: "ModelAdapter",
) -> dict[str, list[int]]:
    """Build the gene -> protein-row-indices map used by ``find_ec_in_db``.

    Tries the conversion table (``uniprotConversion.tsv``) first: if
    that returns different IDs than the gene-name transformation, each
    model gene maps to exactly one UniProt ID -> one protein row.
    Otherwise falls back to gene-name matching, where a single gene
    name may match multiple protein rows.
    """
    if not model_genes:
        return {}

    transformed = adapter.get_uniprot_compatible_genes(model_genes)
    mapped_ids = adapter.get_uniprot_ids_from_table(transformed)

    if mapped_ids != transformed:
        id_to_index = {uid: i for i, uid in enumerate(uniprot_db.ids)}
        result: dict[str, list[int]] = {}
        for original, uid in zip(model_genes, mapped_ids):
            idx = id_to_index.get(uid)
            if idx is not None:
                result[original] = [idx]
        return result

    name_to_indices: dict[str, list[int]] = {}
    for i, g in enumerate(uniprot_db.genes):
        if g:
            name_to_indices.setdefault(g, []).append(i)

    result = {}
    for original, name in zip(model_genes, transformed):
        indices = name_to_indices.get(name)
        if indices:
            result[original] = list(indices)
    return result


def _format_conflict_message(
    conflicts: list[_Conflict],
    uniprot_db: "UniprotDB",
    model: "EcModel",
) -> str:
    """Build the aggregated multi-line message for action='display'."""
    rxn_count = len({c.rxn_idx for c in conflicts})
    lines = [
        f"fill_eccodes_from_database: {len(conflicts)} gene-protein conflict(s) "
        f"found across {rxn_count} reaction(s). Resolve by editing your "
        f"UniProt data, or call fill_eccodes_from_database with "
        f"action='ignore' to silently use the first match per gene:",
    ]
    for c in conflicts:
        rxn_id = model.ec.rxns[c.rxn_idx]
        protein_ids = [uniprot_db.ids[p] for p in c.protein_indices]
        lines.append(
            f"  - rxn {rxn_id!r} / gene {c.gene!r}: "
            f"proteins {', '.join(protein_ids)}"
        )
    return "\n".join(lines)


def get_ec_from_database(
    model: "EcModel",
    uniprot_db: "UniprotDB",
    *,
    ec_rxns: Optional[Iterable[str]] = None,
    action: Literal["display", "ignore"] = "display",
) -> None:
    """Deprecated alias for :func:`fill_eccodes_from_database`.

    Kept for backward compatibility with the original MATLAB name.
    Will be removed in a future release; switch to
    ``fill_eccodes_from_database``.
    """
    import warnings

    warnings.warn(
        "get_ec_from_database is deprecated; use "
        "fill_eccodes_from_database instead. The old name will be "
        "removed in a future release.",
        DeprecationWarning,
        stacklevel=2,
    )
    return fill_eccodes_from_database(
        model, uniprot_db, ec_rxns=ec_rxns, action=action,
    )
