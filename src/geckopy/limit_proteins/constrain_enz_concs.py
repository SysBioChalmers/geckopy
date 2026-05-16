"""Constrain enzyme usage reactions by measured protein concentrations.

Ported from GECKO MATLAB:
src/geckomat/limit_proteins/constrainEnzConcs.m.
"""
from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from ..ec_model.ec_model import EcModel


from ..ec_model.constants import (
    POOL_ID as _PROT_POOL_ID,
    USAGE_PREFIX as _USAGE_PREFIX,
)


def constrain_enz_concs(
    model: "EcModel",
    *,
    remove_constraints: bool = False,
    restrict_to: list[str] | None = None,
) -> None:
    """Constrain ``usage_prot_<enzyme>`` reactions by ``model.ec.concs``.

    For each enzyme in ``model.ec.enzymes`` whose entry in
    ``model.ec.concs`` is non-NaN, the corresponding
    ``usage_prot_<enzyme>`` reaction's upper bound is set to that
    concentration. Enzymes with NaN concentrations have their usage
    reaction reset to the default ``upper_bound = 1000``, drawing
    freely from ``prot_pool``.

    Ported from GECKO MATLAB:
    src/geckomat/limit_proteins/constrainEnzConcs.m.

    MATLAB-COMPAT: GECKO MATLAB's usage reactions go in the reverse
    direction (``prot_<enzyme> -> prot_pool``, bounds ``(-1000, 0)``)
    and the constraint is set as ``lb = -conc``. geckopy's usage
    reactions go forward (``prot_pool -> prot_<enzyme>``, bounds
    ``(0, 1000)``) per the convention established in
    ``protein_pool.py``, so the equivalent constraint is
    ``ub = conc``.

    MATLAB-COMPAT: The MATLAB function applies to full ecModels only
    (gecko-light has no usage reactions). geckopy raises
    ``ValueError`` if any ``usage_prot_<enzyme>`` reaction is
    missing.

    MATLAB-COMPAT: GECKO MATLAB pre-3.2.0 zeroed the ``prot_pool``
    stoichiometry of constrained usage reactions to disconnect them
    from the pool. The current version (and geckopy) just lowers
    the bound; stoichiometry is left at the value set by
    ``add_protein_usage_reactions``.

    Parameters
    ----------
    model
        EcModel with ``model.ec.concs`` populated (typically by
        ``fill_enz_concs``) and the protein pool / usage reaction
        machinery already in place. Mutated in place.
    remove_constraints
        When True, all usage reactions are reset to their default
        ``upper_bound = 1000``, regardless of ``ec.concs``.
        ``ec.concs`` itself is left unchanged.
    restrict_to
        If supplied, only the listed enzyme uniprot IDs are
        processed (others are left untouched). Useful for incremental
        updates from ``Enzyme.concentration`` setters. Unknown IDs
        in the list are silently ignored. Default ``None`` =
        process every enzyme in ``model.ec.enzymes``.

    Raises
    ------
    ValueError
        If the ``prot_pool`` metabolite is missing, or if any
        ``usage_prot_<enzyme>`` reaction is missing for an enzyme
        in ``model.ec.enzymes``.
    """
    if _PROT_POOL_ID not in {m.id for m in model.metabolites}:
        raise ValueError(
            f"Cannot find {_PROT_POOL_ID!r} pseudometabolite. The protein "
            f"pool machinery must be set up before constraining enzyme "
            f"concentrations."
        )

    cobra_rxn_ids = {r.id for r in model.reactions}
    for enz in model.ec.enzymes:
        rxn_id = f"{_USAGE_PREFIX}{enz}"
        if rxn_id not in cobra_rxn_ids:
            raise ValueError(
                f"Usage reaction {rxn_id!r} not found. Usage reactions "
                f"are added by add_protein_usage_reactions; gecko-light "
                f"models do not have them."
            )

    if restrict_to is None:
        index_iter = enumerate(model.ec.enzymes)
    else:
        restrict_set = set(restrict_to)
        index_iter = (
            (i, enz)
            for i, enz in enumerate(model.ec.enzymes)
            if enz in restrict_set
        )

    for i, enz in index_iter:
        rxn = model.reactions.get_by_id(f"{_USAGE_PREFIX}{enz}")
        rxn.upper_bound = 1000.0
        if remove_constraints:
            continue
        conc = model.ec.concs[i]
        if not np.isnan(conc):
            rxn.upper_bound = float(conc)
