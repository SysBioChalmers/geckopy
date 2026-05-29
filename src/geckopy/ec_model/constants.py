"""Shared names for the EC layer (metabolites, reactions, subsystem).

Centralises the magic strings that several modules in geckopy use to
identify enzyme metabolites, usage reactions, the shared protein pool,
and the protein-usage subsystem tag. Promoted from private module-level
constants in ``protein_pool.py`` so that other modules
(``enzyme.py``, ``pfba_enzymes``, parallel ``ec_fva``,
``relax_proteomics_greedy``, SBML I/O, ...) can share one definition.
"""

import re

PROT_PREFIX = "prot_"
USAGE_PREFIX = "usage_prot_"
POOL_ID = "prot_pool"
POOL_EXCHANGE_ID = "prot_pool_exchange"
PROTEIN_USAGE_SUBSYSTEM = "Protein usage"

# Reaction-id suffixes added by the preprocessing pipeline:
#   <base>            forward, single isozyme
#   <base>_EXP_<n>    forward, isozyme n (expand_model)
#   <base>_REV        reverse direction (convert_to_irreversible)
#   <base>_REV_EXP_<n> reverse direction, isozyme n
REV_SUFFIX = "_REV"
_EXP_SUFFIX_RE = re.compile(r"_EXP_\d+$")


def canonicalize_rxn_id(rxn_id: str) -> tuple[str, bool]:
    """Map a split reaction id back to its conventional GEM id.

    Returns ``(base_id, is_reverse)``. Strips a trailing ``_EXP_<n>``
    and then a trailing ``_REV`` (so ``_REV_EXP_<n>`` is handled too).
    Reaction ids that merely *contain* ``_REV`` elsewhere — e.g. a real
    id like ``R_REVERTASE`` — are left intact, unlike a blanket
    ``str.replace("_REV", "")``.
    """
    base = _EXP_SUFFIX_RE.sub("", rxn_id)
    is_reverse = base.endswith(REV_SUFFIX)
    if is_reverse:
        base = base[: -len(REV_SUFFIX)]
    return base, is_reverse
