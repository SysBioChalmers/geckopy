"""Compute enzyme-concentration control coefficients via LP sensitivity.

Ported from GECKO MATLAB:
src/geckomat/limit_proteins/getConcControlCoeffs.m.
"""
from __future__ import annotations

import logging
from typing import TYPE_CHECKING, Optional

import numpy as np

if TYPE_CHECKING:
    from ..ec_model.ec_model import EcModel


from ..ec_model.constants import USAGE_PREFIX

logger = logging.getLogger(__name__)

_GROWTH_DELTA_THRESHOLD = 1e-10

# Solvers we have explicitly confirmed do not expose LP duals via optlang.
# Anything else is assumed to (true for glpk, glpk_exact, gurobi, cplex).
_NO_DUALS_SOLVER_TAGS = ("scipy",)

# One-shot info reporting when the shadow-price path is unavailable.
_FALLBACK_REPORTED: set[str] = set()


def get_conc_control_coeffs(
    model: "EcModel",
    proteins: Optional[list[str]] = None,
    fold_change: float = 2.0,
    limit: float = 0.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute control coefficients for enzyme concentrations.

    For each protein in ``proteins``, the coefficient quantifies how
    much the LP optimum (typically growth) increases per unit increase
    in that enzyme's upper bound. Proteins whose usage flux is at most
    ``limit`` times their upper bound are skipped (no sensitivity to
    expect there).

    The default path uses the ``prot_<id>`` metabolite **shadow price**
    from a single LP solve — equivalent to the per-protein
    finite-difference probe at the margin, but ``n`` times fewer LPs.
    When the active solver is one that does not expose LP duals (notably
    optlang's scipy backend), the function automatically falls back to
    the per-protein finite-difference loop and logs once that it has
    done so. The fold-change probe size is only used in that fallback.

    Ported from GECKO MATLAB:
    src/geckomat/limit_proteins/getConcControlCoeffs.m.

    MATLAB-COMPAT: The MATLAB usage-rxn convention is reverse
    (constraint as ``lb = -conc``); geckopy uses the forward
    convention (``ub = conc``). The control coefficient sign is
    identical either way.

    MATLAB-COMPAT: GECKO MATLAB always uses the 2x finite-difference
    probe; geckopy's default reads the local LP dual (analytic and
    direction-unbiased), with the finite-difference kept as a fallback.

    Parameters
    ----------
    model
        EcModel with a defined LP objective (typically a biomass
        reaction) and ``usage_prot_<enzyme>`` reactions for every
        analysed protein.
    proteins
        UniProt IDs to analyse. Defaults to ``model.ec.enzymes``.
    fold_change
        Probe size for the finite-difference fallback (default 2.0,
        matching MATLAB). Ignored on the shadow-price path.
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
    if n == 0:
        return np.zeros(0, dtype=bool), np.zeros(0, dtype=float)

    if _solver_supports_duals(model):
        return _shadow_price_coeffs(model, proteins, limit)
    _report_finite_difference_fallback(model)
    return _finite_difference_coeffs(model, proteins, fold_change, limit)


def _shadow_price_coeffs(
    model: "EcModel", proteins: list[str], limit: float,
) -> tuple[np.ndarray, np.ndarray]:
    """One-solve control coefficients via the prot_<id> metabolite duals."""
    n = len(proteins)
    enz_mask = np.zeros(n, dtype=bool)
    coeffs = np.zeros(n, dtype=float)

    sol = model.optimize()
    if not _solution_is_optimal(sol):
        return enz_mask, coeffs
    fluxes = sol.fluxes

    cobra_rxn_ids = {r.id for r in model.reactions}
    for i, protein in enumerate(proteins):
        rxn_id = f"{USAGE_PREFIX}{protein}"
        if rxn_id not in cobra_rxn_ids:
            continue
        rxn = model.reactions.get_by_id(rxn_id)
        prev_ub = rxn.upper_bound
        if prev_ub <= 0:
            continue
        usage_flux = float(fluxes.get(rxn_id, 0.0))
        if usage_flux / prev_ub <= limit:
            continue
        enz_mask[i] = True
        # |shadow price| of the prot_<id> metabolite constraint equals the
        # local marginal d(objective)/d(usage_ub). Verified in tests across
        # several scenarios and solvers (glpk, glpk_exact, gurobi).
        sp = abs(model.enzymes.get_by_id(protein).shadow_price)
        if sp > _GROWTH_DELTA_THRESHOLD:
            coeffs[i] = sp

    return enz_mask, coeffs


def _finite_difference_coeffs(
    model: "EcModel",
    proteins: list[str],
    fold_change: float,
    limit: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Per-protein finite-difference fallback (one solve per protein).

    Same contract as :func:`get_conc_control_coeffs`. Kept as the
    fallback for solvers that do not expose LP duals.
    """
    n = len(proteins)
    enz_mask = np.zeros(n, dtype=bool)
    coeffs = np.zeros(n, dtype=float)

    initial_sol = model.optimize()
    if not _solution_is_optimal(initial_sol):
        return enz_mask, coeffs
    initial_growth = initial_sol.objective_value
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
        if not _solution_is_optimal(new_sol):
            continue
        new_growth = new_sol.objective_value
        delta_growth = new_growth - initial_growth
        if delta_growth > _GROWTH_DELTA_THRESHOLD:
            coeffs[i] = delta_growth / (new_ub - prev_ub)

    return enz_mask, coeffs


def _solution_is_optimal(sol) -> bool:
    """Guard against solver states that aren't a genuine LP optimum.

    On an infeasible LP some solvers (notably glpk via optlang) gracefully
    return a non-NaN ``objective_value`` — for our infeasibility-by-binding
    test, glpk reports ``status='infeasible'`` with ``objective_value=1000.0``
    (the lower-bound rhs that made it infeasible) and ``fluxes`` populated
    from the last attempted basis. The previous ``None or isnan`` check
    accepted that, which then caused enz_mask to be set on a non-existent
    optimum. Requiring ``status == "optimal"`` rejects every non-optimal
    state (infeasible, unbounded, suboptimal, ...) cleanly.
    """
    if sol is None or sol.status != "optimal":
        return False
    obj = sol.objective_value
    return obj is not None and not np.isnan(obj)


def _solver_supports_duals(model: "EcModel") -> bool:
    """True iff the model's LP solver exposes constraint duals.

    optlang's scipy backend explicitly raises ``NotImplementedError`` on
    ``constraint.dual``; glpk/glpk_exact/gurobi/cplex all expose duals
    for LP optima. We detect by the solver module name (cheap; no probe
    solve) and assume "supports duals" unless we know otherwise.
    """
    module = type(model.solver).__module__.lower()
    return not any(tag in module for tag in _NO_DUALS_SOLVER_TAGS)


def _report_finite_difference_fallback(model: "EcModel") -> None:
    """Log once per solver that we fell back to the per-protein loop."""
    module = type(model.solver).__module__
    if module in _FALLBACK_REPORTED:
        return
    _FALLBACK_REPORTED.add(module)
    logger.info(
        "get_conc_control_coeffs: solver %r does not expose LP duals; "
        "using the per-protein finite-difference fallback (one solve per "
        "protein, %sx slower). This is reported once per process.",
        module, "n",
    )
