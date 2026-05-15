"""Extract a context-specific ecModel as a subset of a bigger ecModel.

Ported from GECKO MATLAB:
src/geckomat/utilities/getSubsetEcModel.m.
"""
from __future__ import annotations

import copy
import logging
import re
from typing import TYPE_CHECKING

import numpy as np
from scipy import sparse

if TYPE_CHECKING:
    import cobra

    from ..ec_model.ec_model import EcModel


logger = logging.getLogger(__name__)


_REV_SUFFIX = "_REV"
_REV_EXP_INFIX = "_REV_EXP_"
_EXP_RE = re.compile(r"_EXP_\d+")
_USAGE_PREFIX = "usage_prot_"
_POOL_EXCHANGE_ID = "prot_pool_exchange"
_STANDARD_GENE = "standard"
_USAGE_DEFAULT_UB = 1000.0


def get_subset_ec_model(
    big_ec_model: "EcModel",
    small_gem: "cobra.Model",
) -> "EcModel":
    """Trim ``big_ec_model`` to the gene/reaction set of ``small_gem``.

    Both models must have been derived from the same starting GEM.
    Reactions in ``small_gem`` are matched against canonicalised
    (suffix-stripped) reaction ids in ``big_ec_model``; if any
    ``small_gem`` reaction is missing, the function raises.

    Returns a fresh `EcModel`; the input ``big_ec_model`` is not
    mutated.

    Algorithm (mirrors MATLAB):

    1. Copy ``big_ec_model``.
    2. Compute genes to remove (`big.genes` minus `small_gem.genes`,
       preserving the ``standard`` pseudo-gene from
       ``get_standard_kcat``).
    3. Remove reactions blocked by that gene removal (cobra's
       ``GPR.eval(knockouts=...)``).
    4. Remove the genes from `model.genes`.
    5. Trim the corresponding columns of
       ``model.ec.genes/enzymes/mw/sequence/concs/rxn_enz_mat``.
    6. Remove reactions whose canonicalised id is not in
       ``small_gem.reactions``, except ``usage_prot_*`` and
       ``prot_pool_exchange`` which are always preserved.
    7. Trim the corresponding rows of
       ``model.ec.rxns/kcat/source/notes/eccodes/rxn_enz_mat``.

    Ported from GECKO MATLAB:
    src/geckomat/utilities/getSubsetEcModel.m.

    MATLAB-COMPAT: The MATLAB warning checks
    ``lb(usage_prot_*) ~= -1000`` because its usage rxns go reverse;
    geckopy uses forward direction so the equivalent check is
    ``upper_bound != 1000``.

    MATLAB-COMPAT: GECKO MATLAB uses ``dispEM`` which prints AND
    raises; geckopy raises ``ValueError`` directly.

    Parameters
    ----------
    big_ec_model
        The full ecModel covering many genes/reactions.
    small_gem
        The conventional (non-ec) context-specific model whose
        gene/reaction set is the target.

    Returns
    -------
    EcModel
        A new ecModel restricted to `small_gem`'s gene/reaction set,
        with the protein pool / usage rxn machinery preserved.

    Raises
    ------
    ValueError
        If any reaction in ``small_gem`` is missing from
        ``big_ec_model`` (after suffix stripping).
    """
    # --- 1. Validate small_gem reactions are a subset of big's ---
    big_canonical_ids = {
        _canonical(r.id) for r in big_ec_model.reactions
    }
    missing = [
        r.id for r in small_gem.reactions
        if r.id not in big_canonical_ids
    ]
    if missing:
        preview = missing[:5]
        raise ValueError(
            f"{len(missing)} reaction(s) in small_gem not found in "
            f"big_ec_model after suffix stripping (examples: "
            f"{preview}). Were both models derived from the same "
            f"starting GEM?"
        )

    # --- 2. Warn on context-dependent protein constraints ---
    for rxn in big_ec_model.reactions:
        if rxn.id.startswith(_USAGE_PREFIX):
            if rxn.upper_bound != _USAGE_DEFAULT_UB:
                logger.warning(
                    "get_subset_ec_model: big_ec_model has protein-"
                    "concentration constraints (e.g. %s.upper_bound = %g) "
                    "that may not be relevant in the subset model.",
                    rxn.id, rxn.upper_bound,
                )
                break

    # --- 3. Copy the model ---
    # cobra.Model.copy() doesn't know about the `ec` substructure, so
    # deep-copy it explicitly to avoid mutating the input model.
    small = big_ec_model.copy()
    small.ec = copy.deepcopy(big_ec_model.ec)

    # --- 4. Compute genes to remove ---
    small_gem_gene_ids = {g.id for g in small_gem.genes}
    big_gene_ids = {g.id for g in small.genes}
    genes_to_remove = (
        big_gene_ids - small_gem_gene_ids - {_STANDARD_GENE}
    )
    if not genes_to_remove and not _need_rxn_trim(small, small_gem):
        return small  # nothing to do

    # --- 5. Remove reactions blocked by the gene removal ---
    rxns_to_remove_for_genes = []
    for rxn in list(small.reactions):
        if not rxn.gene_reaction_rule:
            continue
        rxn_gene_ids = {g.id for g in rxn.genes}
        if not (rxn_gene_ids & genes_to_remove):
            continue
        if not rxn.gpr.eval(knockouts=genes_to_remove):
            rxns_to_remove_for_genes.append(rxn)
    if rxns_to_remove_for_genes:
        small.remove_reactions(rxns_to_remove_for_genes, remove_orphans=False)

    # --- 6. Remove the genes from model.genes ---
    for gene_id in list(genes_to_remove):
        if gene_id in {g.id for g in small.genes}:
            small.genes.remove(small.genes.get_by_id(gene_id))

    # --- 7. Trim ec per-enzyme fields ---
    if genes_to_remove:
        keep_enz_mask = np.array(
            [g not in genes_to_remove for g in small.ec.genes], dtype=bool,
        )
        small.ec.genes = [
            g for g, k in zip(small.ec.genes, keep_enz_mask) if k
        ]
        small.ec.enzymes = [
            e for e, k in zip(small.ec.enzymes, keep_enz_mask) if k
        ]
        small.ec.mw = small.ec.mw[keep_enz_mask]
        small.ec.sequence = [
            s for s, k in zip(small.ec.sequence, keep_enz_mask) if k
        ]
        if small.ec.concs.size == keep_enz_mask.size:
            small.ec.concs = small.ec.concs[keep_enz_mask]
        if small.ec.rxn_enz_mat.shape[1] == keep_enz_mask.size:
            small.ec.rxn_enz_mat = small.ec.rxn_enz_mat[
                :, keep_enz_mask
            ].tocsr()

    # --- 8. Remove reactions whose canonical id is not in small_gem ---
    small_gem_rxn_ids = {r.id for r in small_gem.reactions}
    rxns_to_remove_for_subset = []
    for rxn in list(small.reactions):
        if rxn.id.startswith(_USAGE_PREFIX) or rxn.id == _POOL_EXCHANGE_ID:
            continue
        if _canonical(rxn.id) not in small_gem_rxn_ids:
            rxns_to_remove_for_subset.append(rxn)
    if rxns_to_remove_for_subset:
        small.remove_reactions(
            rxns_to_remove_for_subset, remove_orphans=True,
        )

    # --- 9. Trim ec per-rxn fields ---
    cobra_rxn_ids = {r.id for r in small.reactions}
    keep_rxn_mask = np.array(
        [rid in cobra_rxn_ids for rid in small.ec.rxns], dtype=bool,
    )
    if not keep_rxn_mask.all():
        small.ec.rxns = [
            r for r, k in zip(small.ec.rxns, keep_rxn_mask) if k
        ]
        small.ec.kcat = small.ec.kcat[keep_rxn_mask]
        small.ec.source = [
            s for s, k in zip(small.ec.source, keep_rxn_mask) if k
        ]
        small.ec.notes = [
            n for n, k in zip(small.ec.notes, keep_rxn_mask) if k
        ]
        small.ec.eccodes = [
            c for c, k in zip(small.ec.eccodes, keep_rxn_mask) if k
        ]
        if small.ec.rxn_enz_mat.shape[0] == keep_rxn_mask.size:
            small.ec.rxn_enz_mat = small.ec.rxn_enz_mat[
                keep_rxn_mask, :
            ].tocsr()

    return small


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #

def _canonical(rxn_id: str) -> str:
    """Strip ``_REV`` and ``_EXP_<N>`` suffixes."""
    rid = rxn_id
    if rid.endswith(_REV_SUFFIX) or _REV_EXP_INFIX in rid:
        rid = rid.replace(_REV_SUFFIX, "")
    return _EXP_RE.sub("", rid)


def _need_rxn_trim(big: "EcModel", small_gem: "cobra.Model") -> bool:
    """True if there's at least one reaction in big that's not in
    small_gem and not part of the protein machinery."""
    small_gem_rxn_ids = {r.id for r in small_gem.reactions}
    for rxn in big.reactions:
        if rxn.id.startswith(_USAGE_PREFIX) or rxn.id == _POOL_EXCHANGE_ID:
            continue
        if _canonical(rxn.id) not in small_gem_rxn_ids:
            return True
    return False
