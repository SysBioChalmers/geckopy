"""Flux variability analysis for ecModels, mapped back to conventional rxns.

Ported from GECKO MATLAB:
src/geckomat/utilities/ecFVA.m.
"""
from __future__ import annotations

import re
from collections import defaultdict
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

from .map_rxns_to_conv import map_rxns_to_conv

if TYPE_CHECKING:
    import cobra

    from ..ec_model.ec_model import EcModel


_REV_SUFFIX = "_REV"
_REV_EXP_INFIX = "_REV_EXP_"
_EXP_RE = re.compile(r"_EXP_\d+")


def ec_fva(
    ec_model: "EcModel",
    model: "cobra.Model",
) -> pd.DataFrame:
    """Flux variability analysis on an ecModel, mapped to conventional rxns.

    For each canonical reaction (after stripping ``_REV`` and
    ``_EXP_<N>`` suffixes), maximises and minimises the combined
    (forward minus reverse) flux of all its ec variants, then
    aggregates the resulting per-ec-rxn min/max across all
    iterations. Finally translates the (n_ec_rxns, 2) min/max array
    to conventional-rxn space via ``map_rxns_to_conv``; if any
    mapped pair has min > max (a direction flip from the ``_REV``
    handling), they are swapped.

    Ported from GECKO MATLAB:
    src/geckomat/utilities/ecFVA.m.

    MATLAB-COMPAT: GECKO MATLAB uses ``parfor`` to parallelise the
    per-conv-rxn LP solves. geckopy runs the loop serially. Adding
    parallelism (e.g. via ``joblib`` or ``multiprocessing``) would
    require the cobra model to be picklable, which is non-trivial
    for shared solver state. Tracked as a future optimisation.

    MATLAB-COMPAT: GECKO MATLAB returns ``(minFlux, maxFlux)``
    aligned to ``model.rxns``. geckopy returns a ``pd.DataFrame``
    indexed by reaction id with ``min_flux`` / ``max_flux`` columns
    (more directly useful in Python).

    Parameters
    ----------
    ec_model
        EcModel with `_REV` and `_EXP_<N>` suffixes from
        `expand_model` and `convert_to_irreversible`.
    model
        The non-ec model to which the fluxes are mapped. Must be the
        same model that was used to build ``ec_model``.

    Returns
    -------
    pandas.DataFrame
        Index: ``rxn_id`` (in ``model.reactions`` order).
        Columns: ``min_flux``, ``max_flux``.
    """
    n_ec = len(ec_model.reactions)
    if n_ec == 0:
        return pd.DataFrame(columns=["min_flux", "max_flux"])

    # Group ec rxn indices by canonical (stripped) id.
    canonical_groups: dict[str, list[int]] = defaultdict(list)
    is_reverse = [False] * n_ec
    for i, rxn in enumerate(ec_model.reactions):
        rid = rxn.id
        if rid.endswith(_REV_SUFFIX) or _REV_EXP_INFIX in rid:
            is_reverse[i] = True
            rid = rid.replace(_REV_SUFFIX, "")
        rid = _EXP_RE.sub("", rid)
        canonical_groups[rid].append(i)

    canonical_ids = list(canonical_groups.keys())
    n_groups = len(canonical_ids)

    sol_max_all = np.full((n_ec, n_groups), np.nan, dtype=float)
    sol_min_all = np.full((n_ec, n_groups), np.nan, dtype=float)

    for k, conv_id in enumerate(canonical_ids):
        rxn_indices = canonical_groups[conv_id]
        forward = [i for i in rxn_indices if not is_reverse[i]]
        reverse = [i for i in rxn_indices if is_reverse[i]]

        # Maximise the canonical (forward minus reverse) flux.
        sol_max_all[:, k] = _solve_with_objective(
            ec_model, forward_pos=forward, reverse_neg=reverse, sense="max",
        )
        # Minimise (= maximise the negated objective).
        sol_min_all[:, k] = _solve_with_objective(
            ec_model, forward_pos=reverse, reverse_neg=forward, sense="max",
        )

    # Per ec rxn: min over min-direction solutions, max over max-direction.
    with np.errstate(invalid="ignore"):
        if np.all(np.isnan(sol_min_all), axis=1).any():
            # Avoid the all-NaN-slice warning by using nanmin/nanmax with care.
            pass
    min_flux_ec = _safe_nan_reduce(sol_min_all, np.nanmin)
    max_flux_ec = _safe_nan_reduce(sol_max_all, np.nanmax)

    # Map to conventional model.
    combined = np.column_stack([min_flux_ec, max_flux_ec])
    mapped = map_rxns_to_conv(ec_model, model, combined).mapped_flux
    min_flux = mapped[:, 0].copy()
    max_flux = mapped[:, 1].copy()

    # Direction may have flipped during _REV handling; restore min <= max.
    swap = min_flux > max_flux
    if swap.any():
        min_flux[swap], max_flux[swap] = max_flux[swap], min_flux[swap].copy()

    return pd.DataFrame(
        {"min_flux": min_flux, "max_flux": max_flux},
        index=[r.id for r in model.reactions],
    ).rename_axis("rxn_id")


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #

def _solve_with_objective(
    ec_model: "EcModel",
    *,
    forward_pos: list[int],
    reverse_neg: list[int],
    sense: str,
) -> np.ndarray:
    """Solve LP with `forward_pos` rxns at coeff +1, `reverse_neg` at -1.

    Returns the full ec-rxn flux vector, or all-NaN on infeasibility.
    Uses `with ec_model:` so objective changes are reverted.
    """
    n = len(ec_model.reactions)
    out = np.full(n, np.nan, dtype=float)

    if not forward_pos and not reverse_neg:
        return out

    rxns = list(ec_model.reactions)
    with ec_model:
        # Reset all objective coefficients first.
        for rxn in rxns:
            rxn.objective_coefficient = 0.0
        for i in forward_pos:
            rxns[i].objective_coefficient = 1.0
        for i in reverse_neg:
            rxns[i].objective_coefficient = -1.0
        ec_model.objective_direction = sense
        sol = ec_model.optimize()
        if sol.objective_value is None or np.isnan(sol.objective_value):
            return out
        for i, rxn in enumerate(rxns):
            try:
                out[i] = float(sol.fluxes[rxn.id])
            except KeyError:
                out[i] = 0.0

    return out


def _safe_nan_reduce(arr: np.ndarray, reducer) -> np.ndarray:
    """Apply `nanmin` / `nanmax` along axis=1, returning 0 for all-NaN rows
    instead of NaN with a warning."""
    if arr.size == 0:
        return np.empty(arr.shape[0], dtype=float)
    mask_all_nan = np.all(np.isnan(arr), axis=1)
    out = np.zeros(arr.shape[0], dtype=float)
    if (~mask_all_nan).any():
        with np.errstate(invalid="ignore"):
            import warnings
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                out[~mask_all_nan] = reducer(arr[~mask_all_nan, :], axis=1)
    return out
