"""Flux variability analysis (FVA) on an ecModel.

Standard FVA asks, for each reaction in a model, what's the
smallest and largest flux that reaction can carry while keeping
the model feasible at its current objective? The output is a
``(min_flux, max_flux)`` pair per reaction, useful for finding
reactions whose flux is tightly constrained vs reactions with a
lot of slack.

This function does the same thing on an ecModel, with two twists:

- Reactions that were split per isozyme (``_EXP_<N>`` suffix) or
  per direction (``_REV`` suffix) by ``make_ec_model`` are
  collapsed back to their original cobra-reaction id before
  reporting, so the output table has one row per starting-GEM
  reaction.
- The per-reaction LP solves can run in parallel via
  ``multiprocessing`` (``n_proc`` argument). On real-size
  ecModels this matters — yeast-GEM has ~4000 conv reactions,
  each needing two LP solves, so a serial run can be slow.

Ported from GECKO MATLAB: src/geckomat/utilities/ecFVA.m.
"""
from __future__ import annotations

import multiprocessing as mp
import pickle
import sys
from collections import defaultdict
from typing import TYPE_CHECKING, Optional

import cobra
import numpy as np
import pandas as pd

from ..ec_model.constants import canonicalize_rxn_id
from .map_rxns_to_conv import map_rxns_to_conv

try:
    from tqdm import tqdm
except ImportError:  # pragma: no cover - tqdm is a soft dep
    def tqdm(it, **_kwargs):  # type: ignore[no-redef]
        return it

if TYPE_CHECKING:
    from ..ec_model.ec_model import EcModel


# --------------------------------------------------------------------------- #
# Worker globals (process-local; populated by _init_worker)
# --------------------------------------------------------------------------- #

_WORKER_MODEL: "EcModel | None" = None


def _init_worker(model_pickle_bytes: bytes) -> None:
    """Pool initializer: store one EcModel copy per worker process."""
    global _WORKER_MODEL
    _WORKER_MODEL = pickle.loads(model_pickle_bytes)


def _fva_step(arg):
    """Solve max+min objectives for one conv reaction in the worker model.

    arg: (conv_id, forward_indices, reverse_indices).
    Returns: (conv_id, max_vector, min_vector) where vectors have one
    entry per ec reaction in the worker model.
    """
    assert _WORKER_MODEL is not None
    conv_id, forward, reverse = arg
    max_vec = _solve_for_conv_rxn(
        _WORKER_MODEL, forward_pos=forward, reverse_neg=reverse, sense="max",
    )
    min_vec = _solve_for_conv_rxn(
        _WORKER_MODEL, forward_pos=reverse, reverse_neg=forward, sense="max",
    )
    return conv_id, max_vec, min_vec


# --------------------------------------------------------------------------- #
# Public entry point
# --------------------------------------------------------------------------- #

def ec_fva(
    ec_model: "EcModel",
    model: "cobra.Model",
    *,
    progress: bool = True,
    n_proc: Optional[int] = None,
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
    per-conv-rxn LP solves. geckopy now supports parallelism via
    ``multiprocessing`` with the ``n_proc`` parameter; the default
    (``None``) resolves to ``cobra.Configuration().processes``. Pass
    ``n_proc=1`` to force the serial path.

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
    progress
        Show a tqdm progress bar over the per-canonical-rxn LP
        solves. Defaults to True; set False to silence output.
    n_proc
        Number of worker processes. Defaults to
        ``cobra.Configuration().processes``. Set to 1 to run the
        original serial path; values >= 2 parallelise with a spawn
        Pool that pickles ``ec_model`` once per worker.

    Returns
    -------
    pandas.DataFrame
        Index: ``rxn_id`` (in ``model.reactions`` order).
        Columns: ``min_flux``, ``max_flux``.
    """
    n_ec = len(ec_model.reactions)
    if n_ec == 0:
        return pd.DataFrame(columns=["min_flux", "max_flux"])

    canonical_ids, group_forward, group_reverse = _list_conv_rxn_groups(ec_model)
    n_groups = len(canonical_ids)

    sol_max_all = np.full((n_ec, n_groups), np.nan, dtype=float)
    sol_min_all = np.full((n_ec, n_groups), np.nan, dtype=float)

    if n_proc is None:
        n_proc = cobra.Configuration().processes
    n_proc = max(1, int(n_proc))

    if n_proc == 1:
        # Serial path (preserves the prior behaviour for backwards-compat).
        iterator = enumerate(canonical_ids)
        if progress:
            iterator = enumerate(
                tqdm(canonical_ids, desc="ec_fva", unit="rxn")
            )
        for k, conv_id in iterator:
            sol_max_all[:, k] = _solve_for_conv_rxn(
                ec_model,
                forward_pos=group_forward[k],
                reverse_neg=group_reverse[k],
                sense="max",
            )
            sol_min_all[:, k] = _solve_for_conv_rxn(
                ec_model,
                forward_pos=group_reverse[k],
                reverse_neg=group_forward[k],
                sense="max",
            )
    else:
        # Parallel path. "fork" is the most efficient context but is
        # unavailable on Windows; "spawn" is required there and works
        # on POSIX too at the cost of re-importing modules per worker.
        # Some WSL kernels deadlock on spawn-context multiprocessing,
        # so we prefer fork on POSIX where it's available.
        ctx_name = "fork" if sys.platform != "win32" else "spawn"
        ctx = mp.get_context(ctx_name)
        pickled = pickle.dumps(ec_model)
        tasks = [
            (cid, group_forward[k], group_reverse[k])
            for k, cid in enumerate(canonical_ids)
        ]
        index_by_id = {cid: k for k, cid in enumerate(canonical_ids)}
        chunk = max(1, len(tasks) // (n_proc * 4))
        with ctx.Pool(
            n_proc, initializer=_init_worker, initargs=(pickled,),
        ) as pool:
            iterator = pool.imap_unordered(_fva_step, tasks, chunksize=chunk)
            if progress:
                iterator = tqdm(
                    iterator, total=len(tasks), desc="ec_fva", unit="rxn",
                )
            for conv_id, max_vec, min_vec in iterator:
                k = index_by_id[conv_id]
                sol_max_all[:, k] = max_vec
                sol_min_all[:, k] = min_vec

    return _postprocess(sol_min_all, sol_max_all, ec_model, model)


# --------------------------------------------------------------------------- #
# Refactored helpers (top-level so they're picklable for spawn workers)
# --------------------------------------------------------------------------- #

def _list_conv_rxn_groups(
    ec_model: "EcModel",
) -> tuple[list[str], list[list[int]], list[list[int]]]:
    """Group ec rxn indices by canonical (stripped) id.

    Returns
    -------
    canonical_ids : list[str]
        Unique canonical IDs in first-seen order.
    group_forward : list[list[int]]
        Indices of forward variants per canonical id.
    group_reverse : list[list[int]]
        Indices of reverse (``_REV``) variants per canonical id.
    """
    canonical_groups: dict[str, list[int]] = defaultdict(list)
    is_reverse: list[bool] = [False] * len(ec_model.reactions)
    for i, rxn in enumerate(ec_model.reactions):
        rid, rev = canonicalize_rxn_id(rxn.id)
        is_reverse[i] = rev
        canonical_groups[rid].append(i)

    canonical_ids = list(canonical_groups.keys())
    group_forward: list[list[int]] = []
    group_reverse: list[list[int]] = []
    for cid in canonical_ids:
        indices = canonical_groups[cid]
        group_forward.append([i for i in indices if not is_reverse[i]])
        group_reverse.append([i for i in indices if is_reverse[i]])
    return canonical_ids, group_forward, group_reverse


def _solve_for_conv_rxn(
    ec_model: "EcModel",
    *,
    forward_pos: list[int],
    reverse_neg: list[int],
    sense: str,
) -> np.ndarray:
    """Solve LP with ``forward_pos`` rxns at coeff +1, ``reverse_neg``
    at -1.

    Returns the full ec-rxn flux vector, or all-NaN on infeasibility.
    Uses ``with ec_model:`` so objective changes are reverted.
    """
    n = len(ec_model.reactions)
    out = np.full(n, np.nan, dtype=float)

    if not forward_pos and not reverse_neg:
        return out

    rxns = list(ec_model.reactions)
    with ec_model:
        # Assigning an objective dict zeroes every other reaction's
        # coefficient in one step (reverted on context exit), instead of
        # looping over all reactions to reset them.
        objective = {rxns[i]: 1.0 for i in forward_pos}
        objective.update({rxns[i]: -1.0 for i in reverse_neg})
        ec_model.objective = objective
        ec_model.objective_direction = sense
        sol = ec_model.optimize()
        if sol.status != "optimal":
            return out
        for i, rxn in enumerate(rxns):
            try:
                out[i] = float(sol.fluxes[rxn.id])
            except KeyError:
                out[i] = 0.0
    return out


def _postprocess(
    sol_min_all: np.ndarray,
    sol_max_all: np.ndarray,
    ec_model: "EcModel",
    model: "cobra.Model",
) -> pd.DataFrame:
    """Reduce per-ec-rxn min/max across solves, then map to conv space."""
    min_flux_ec = _safe_nan_reduce(sol_min_all, np.nanmin)
    max_flux_ec = _safe_nan_reduce(sol_max_all, np.nanmax)

    combined = np.column_stack([min_flux_ec, max_flux_ec])
    mapped = map_rxns_to_conv(ec_model, model, combined).mapped_flux
    min_flux = mapped[:, 0].copy()
    max_flux = mapped[:, 1].copy()

    swap = min_flux > max_flux
    if swap.any():
        min_flux[swap], max_flux[swap] = max_flux[swap], min_flux[swap].copy()

    return pd.DataFrame(
        {"min_flux": min_flux, "max_flux": max_flux},
        index=[r.id for r in model.reactions],
    ).rename_axis("rxn_id")


def _safe_nan_reduce(arr: np.ndarray, reducer) -> np.ndarray:
    """Apply ``nanmin`` / ``nanmax`` along axis=1, returning 0 for
    all-NaN rows instead of NaN with a warning."""
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
