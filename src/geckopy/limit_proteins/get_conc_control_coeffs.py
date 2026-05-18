"""Compute enzyme-concentration control coefficients via LP sensitivity.

Ported from GECKO MATLAB:
src/geckomat/limit_proteins/getConcControlCoeffs.m.
"""
from __future__ import annotations

from typing import TYPE_CHECKING, Optional

import numpy as np

if TYPE_CHECKING:
    from ..ec_model.ec_model import EcModel


from ..ec_model.constants import USAGE_PREFIX

_GROWTH_DELTA_THRESHOLD = 1e-10


def get_conc_control_coeffs(
    model: "EcModel",
    proteins: Optional[list[str]] = None,
    fold_change: float = 2.0,
    limit: float = 0.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute control coefficients for enzyme concentrations.

    For each protein in ``proteins``, the coefficient quantifies how
    much the LP optimum (typically growth) increases when the enzyme's
    upper bound is multiplied by ``fold_change``. Proteins whose usage
    flux is at most ``limit`` times their upper bound are skipped (no
    sensitivity to expect there).

    Ported from GECKO MATLAB:
    src/geckomat/limit_proteins/getConcControlCoeffs.m.

    MATLAB-COMPAT: The MATLAB usage-rxn convention is reverse
    (constraint as ``lb = -conc``); geckopy uses the forward
    convention (``ub = conc``). The control coefficient sign is
    identical either way.

    MATLAB-COMPAT: GECKO MATLAB uses ``solveLP`` with hot-starts for
    speed; geckopy uses cobra's ``model.optimize()`` inside a
    ``with model:`` context so each per-protein change is reverted
    cleanly. Hot-starts are not exposed at this level by cobra.

    Parameters
    ----------
    model
        EcModel with a defined LP objective (typically a biomass
        reaction) and ``usage_prot_<enzyme>`` reactions for every
        analysed protein.
    proteins
        UniProt IDs to analyse. Defaults to ``model.ec.enzymes``.
    fold_change
        Factor by which the upper bound is multiplied during the
        sensitivity probe. Default 2.0 matches MATLAB.
    limit
        Skip proteins whose usage flux is at most ``limit`` times
        their upper bound. Default 0.0 includes any protein with
        non-zero flux.

    Returns
    -------
    enz_mask
        Boolean array of length ``len(proteins)``; True where the
        protein was analysed.
    control_coeffs
        Float array of length ``len(proteins)``; 0.0 where the
        protein was either skipped or showed no sensitivity.
    """
    if proteins is None:
        proteins = list(model.ec.enzymes)

    n = len(proteins)
    enz_mask = np.zeros(n, dtype=bool)
    coeffs = np.zeros(n, dtype=float)

    if n == 0:
        return enz_mask, coeffs

    initial_sol = model.optimize()
    initial_growth = initial_sol.objective_value
    if initial_growth is None or np.isnan(initial_growth):
        return enz_mask, coeffs
    initial_fluxes = initial_sol.fluxes

    cobra_rxn_ids = {r.id for r in model.reactions}

    for i, protein in enumerate(proteins):
        rxn_id = f"{USAGE_PREFIX}{protein}"
        if rxn_id not in cobra_rxn_ids:
            continue

        rxn = model.reactions.get_by_id(rxn_id)
        prev_ub = rxn.upper_bound
        if prev_ub <= 0:
            continue

        usage_flux = float(initial_fluxes.get(rxn_id, 0.0))
        if usage_flux / prev_ub <= limit:
            continue

        enz_mask[i] = True
        new_ub = prev_ub * fold_change
        with model:
            rxn.upper_bound = new_ub
            new_sol = model.optimize()
        new_growth = new_sol.objective_value
        if new_growth is None or np.isnan(new_growth):
            continue
        delta_growth = new_growth - initial_growth
        if delta_growth > _GROWTH_DELTA_THRESHOLD:
            coeffs[i] = delta_growth / (new_ub - prev_ub)

    return enz_mask, coeffs
