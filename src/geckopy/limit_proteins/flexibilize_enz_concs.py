"""Loosen proteomics constraints until the model can grow.

After ``constrain_enz_concs`` pins each measured enzyme to its
proteomics value, the model often can't reach the experimental
growth rate any more — the measurements are noisy, and a few
enzymes end up too tightly capped to feed the network.

``flexibilize_enz_concs`` fixes that automatically. It computes a
control coefficient per enzyme (how much each enzyme's bound is
limiting growth), picks the worst one, multiplies its upper bound
by ``fold_change``, re-solves, and repeats until growth reaches
the target. The result records which enzymes were loosened and by
how much, so a curator can review.

A related function is ``relax_proteomics_greedy`` in the same
package: same goal but a different algorithm (ranks by absolute
shadow price instead of control coefficient). Pick the one that
converges faster on your model.

Ported from GECKO MATLAB:
src/geckomat/limit_proteins/flexibilizeEnzConcs.m.
"""
from __future__ import annotations

import logging
import math
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Optional

import numpy as np

from .get_conc_control_coeffs import get_conc_control_coeffs

if TYPE_CHECKING:
    from ..ec_model.ec_model import EcModel


logger = logging.getLogger(__name__)

from ..ec_model.constants import (
    POOL_EXCHANGE_ID as _POOL_EXCHANGE_ID,
    USAGE_PREFIX as _USAGE_PREFIX,
)
_POOL_TARGET_NAME = "prot_pool"
_GROWTH_DELTA_THRESHOLD = 1e-3
_CONTROL_COEFF_LIMIT = 0.75


@dataclass
class FlexEnzResult:
    """Information about which enzyme constraints were flexibilized.

    All arrays are sorted by ``ratio_incr`` descending. May contain a
    final entry with ``uniprot_ids[-1] == "prot_pool"`` if the protein
    pool exchange itself was relaxed.
    """

    uniprot_ids: list[str] = field(default_factory=list)
    old_concs: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=float))
    flex_concs: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=float))
    ratio_incr: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=float))
    frequence: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=int))


def flexibilize_enz_concs(
    model: "EcModel",
    *,
    exp_growth: Optional[float] = None,
    fold_change: float = 2.0,
    iter_per_enzyme: int = 5,
    verbose: bool = True,
) -> FlexEnzResult:
    """Iteratively relax enzyme concentration constraints until the
    model reaches ``exp_growth``.

    Each iteration uses ``get_conc_control_coeffs`` to find the most
    limiting enzyme, then multiplies its ``usage_prot_<enzyme>``
    upper bound by ``(1 + fold_change * iteration_count)``. If no
    enzyme is limiting (all control coefficients zero), the function
    tries relaxing ``prot_pool_exchange``.

    After the loop, a refinement pass solves
    ``min prot_pool_exchange  s.t. bio_rxn >= exp_growth`` to find
    the minimal usage that supports growth. For each flexibilized
    enzyme, if the actual usage at this minimum is less than its
    original concentration, the original constraint is restored
    (i.e. that enzyme is dropped from the result). This avoids
    over-relaxing enzymes that the prior iterative pass increased
    more than necessary.

    Ported from GECKO MATLAB:
    src/geckomat/limit_proteins/flexibilizeEnzConcs.m.

    MATLAB-COMPAT: MATLAB's usage rxns go reverse (constraint as
    ``lb = -conc``); geckopy uses forward direction (``ub = conc``).

    MATLAB-COMPAT: MATLAB uses RAVEN's ``setParam('var', ...)`` soft
    constraint when refining the result. geckopy uses a hard
    ``bio_rxn.lower_bound = exp_growth``; if ``exp_growth`` is
    achievable (which the loop ensured) the two are equivalent.

    MATLAB-COMPAT: MATLAB's ``modelAdapter`` arg is dropped. The
    organism-specific ``gR_exp`` defaults to
    ``model.adapter.params.gr_exp``.

    MATLAB-COMPAT: The MATLAB ``var`` constraint with weight 0.5
    allows the LP to slightly miss ``expGrowth`` if needed. geckopy
    enforces ``exp_growth`` exactly.

    Parameters
    ----------
    model
        EcModel with `prot_pool_exchange` and ``usage_prot_*``
        reactions installed and ``model.ec.concs`` populated for at
        least some enzymes (typically by ``fill_enz_concs`` and
        ``constrain_enz_concs``). Mutated in place; ``ec.concs``
        itself is left unchanged.
    exp_growth
        Target growth rate. Defaults to
        ``model.adapter.params.gr_exp``.
    fold_change
        Per-iteration upper-bound multiplier for the limiting enzyme.
    iter_per_enzyme
        Maximum iterations per enzyme before warning and breaking.
        ``0`` means no limit (matches MATLAB).
    verbose
        Whether per-iteration progress is logged at INFO.

    Returns
    -------
    FlexEnzResult

    Raises
    ------
    ValueError
        If no enzyme has a measured concentration (nothing to
        flexibilize), or if ``model.adapter`` is None and
        ``exp_growth`` is also None.
    """
    if exp_growth is None:
        if model.adapter is None or model.adapter.params.gr_exp is None:
            raise ValueError(
                "exp_growth not provided and model.adapter.params.gr_exp "
                "is unavailable."
            )
        exp_growth = float(model.adapter.params.gr_exp)

    if iter_per_enzyme == 0:
        iter_per_enzyme = math.inf

    bio_rxn_id = model.adapter.params.bio_rxn
    bio_rxn = model.reactions.get_by_id(bio_rxn_id)
    if bio_rxn.upper_bound < exp_growth:
        logger.info(
            "flexibilize_enz_concs: bio_rxn upper_bound (%g) was below "
            "exp_growth (%g); raising it.",
            bio_rxn.upper_bound, exp_growth,
        )
        bio_rxn.upper_bound = exp_growth

    measured_idx = np.where(~np.isnan(model.ec.concs))[0]
    if len(measured_idx) == 0:
        raise ValueError(
            "No enzyme concentrations are measured (model.ec.concs is "
            "all-NaN). Run fill_enz_concs and constrain_enz_concs first."
        )
    proteins = [model.ec.enzymes[i] for i in measured_idx]
    n_proteins = len(proteins)
    frequence = np.zeros(n_proteins, dtype=int)

    initial_sol = model.optimize()
    pred_growth = float(initial_sol.objective_value or 0.0)

    flex_break = False
    pool_old: Optional[float] = None
    pool_new: Optional[float] = None

    while pred_growth < exp_growth:
        _, control_coeffs = get_conc_control_coeffs(
            model,
            proteins=proteins,
            fold_change=fold_change,
            limit=_CONTROL_COEFF_LIMIT,
        )

        if not (control_coeffs > 0).any():
            pool_old, pool_new, pool_pred_growth = _try_relax_pool(
                model, bio_rxn, exp_growth, pred_growth, verbose,
            )
            if pool_pred_growth is not None:
                pred_growth = pool_pred_growth
            break

        max_idx = int(np.argmax(control_coeffs))
        frequence[max_idx] += 1
        if frequence[max_idx] > iter_per_enzyme:
            logger.warning(
                "flexibilize_enz_concs: iter limit reached for protein %r; "
                "consider revising its kcat.",
                proteins[max_idx],
            )
            flex_break = True
            break

        increase = fold_change * frequence[max_idx]
        ec_idx = measured_idx[max_idx]
        usage_rxn = model.reactions.get_by_id(
            f"{_USAGE_PREFIX}{proteins[max_idx]}"
        )
        usage_rxn.upper_bound = float(model.ec.concs[ec_idx]) * (1.0 + increase)

        new_sol = model.optimize()
        pred_growth = float(new_sol.objective_value or 0.0)
        if verbose:
            logger.info(
                "flexibilize_enz_concs: protein %r relaxed; growth = %g",
                proteins[max_idx], pred_growth,
            )

    return _build_result(
        model, proteins, measured_idx, frequence,
        flex_break, exp_growth, pool_old, pool_new,
    )


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #

def _try_relax_pool(
    model: "EcModel",
    bio_rxn,
    exp_growth: float,
    pred_growth: float,
    verbose: bool,
) -> tuple[Optional[float], Optional[float], Optional[float]]:
    """Try relaxing prot_pool_exchange. Returns (old_pool, new_pool, growth)
    if a relaxation was applied; (None, None, None) otherwise."""
    pool_rxn = model.reactions.get_by_id(_POOL_EXCHANGE_ID)
    old_pool = pool_rxn.upper_bound

    with model:
        pool_rxn.upper_bound = 1000.0
        bio_rxn.upper_bound = exp_growth
        sol = model.optimize()
        candidate_growth = float(sol.objective_value or 0.0)

    if candidate_growth - pred_growth <= _GROWTH_DELTA_THRESHOLD:
        if verbose:
            logger.info(
                "flexibilize_enz_concs: protein pool was not limiting "
                "either. Maximum growth reachable: %g", candidate_growth,
            )
        return None, None, None

    with model:
        pool_rxn.upper_bound = 1000.0
        bio_rxn.lower_bound = candidate_growth
        model.objective = _POOL_EXCHANGE_ID
        model.objective_direction = "min"
        sol = model.optimize()
        new_pool = float(sol.objective_value or 0.0)

    pool_rxn.upper_bound = new_pool

    if verbose:
        logger.info(
            "flexibilize_enz_concs: prot_pool_exchange upper_bound raised "
            "from %g to %g, enabling growth %g.",
            old_pool, new_pool, candidate_growth,
        )

    return old_pool, new_pool, candidate_growth


def _build_result(
    model: "EcModel",
    proteins: list[str],
    measured_idx: np.ndarray,
    frequence: np.ndarray,
    flex_break: bool,
    exp_growth: float,
    pool_old: Optional[float],
    pool_new: Optional[float],
) -> FlexEnzResult:
    """Build the FlexEnzResult and apply the post-loop UB refinement."""
    flex_mask = frequence > 0
    flex_indices = np.where(flex_mask)[0]
    if len(flex_indices) == 0:
        return _maybe_add_pool(FlexEnzResult(), pool_old, pool_new)

    flex_proteins = [proteins[i] for i in flex_indices]
    usage_rxn_ids = [f"{_USAGE_PREFIX}{p}" for p in flex_proteins]
    old_concs_arr = np.array(
        [float(model.ec.concs[measured_idx[i]]) for i in flex_indices],
        dtype=float,
    )

    if not flex_break:
        # Refinement pass: minimize protein pool with bio_rxn at exp_growth,
        # then drop enzymes that don't actually need extra concentration.
        bio_rxn_id = model.adapter.params.bio_rxn
        with model:
            model.reactions.get_by_id(bio_rxn_id).lower_bound = exp_growth
            model.objective = _POOL_EXCHANGE_ID
            model.objective_direction = "min"
            sol = model.optimize()
            new_concs = np.array(
                [float(sol.fluxes[rid]) for rid in usage_rxn_ids],
                dtype=float,
            )

        keep_old_mask = new_concs < old_concs_arr
        for i, rxn_id in enumerate(usage_rxn_ids):
            ub = (
                old_concs_arr[i] if keep_old_mask[i] else new_concs[i]
            )
            model.reactions.get_by_id(rxn_id).upper_bound = float(ub)

        survive = ~keep_old_mask
        result = FlexEnzResult(
            uniprot_ids=[
                flex_proteins[i] for i in range(len(flex_proteins)) if survive[i]
            ],
            old_concs=old_concs_arr[survive],
            flex_concs=new_concs[survive],
            frequence=frequence[flex_indices][survive],
        )
    else:
        # flex_break: take the relaxed UBs as-is.
        flex_concs = np.array(
            [
                float(model.reactions.get_by_id(rid).upper_bound)
                for rid in usage_rxn_ids
            ],
            dtype=float,
        )
        result = FlexEnzResult(
            uniprot_ids=flex_proteins,
            old_concs=old_concs_arr,
            flex_concs=flex_concs,
            frequence=frequence[flex_indices],
        )

    return _finalize_and_sort(_maybe_add_pool(result, pool_old, pool_new))


def _maybe_add_pool(
    result: FlexEnzResult,
    pool_old: Optional[float],
    pool_new: Optional[float],
) -> FlexEnzResult:
    if pool_old is None or pool_new is None:
        return result
    result.uniprot_ids = list(result.uniprot_ids) + [_POOL_TARGET_NAME]
    result.old_concs = np.concatenate([result.old_concs, [pool_old]])
    result.flex_concs = np.concatenate([result.flex_concs, [pool_new]])
    result.frequence = np.concatenate([result.frequence, [1]])
    return result


def _finalize_and_sort(result: FlexEnzResult) -> FlexEnzResult:
    """Compute ratio_incr and sort all arrays by it descending."""
    if not result.uniprot_ids:
        return result
    ratios = np.where(
        result.old_concs > 0,
        result.flex_concs / result.old_concs,
        np.inf,
    )
    order = np.argsort(-ratios)  # descending
    return FlexEnzResult(
        uniprot_ids=[result.uniprot_ids[i] for i in order],
        old_concs=result.old_concs[order],
        flex_concs=result.flex_concs[order],
        ratio_incr=ratios[order],
        frequence=result.frequence[order],
    )
