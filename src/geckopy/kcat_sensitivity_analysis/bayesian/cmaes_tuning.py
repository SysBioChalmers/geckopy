"""CMA-ES kcat tuning: an optimiser in place of ABC-SMC's sampler.

Internal validation (``docs/internal/bayesian_tuning_handover.md``,
"Three methods, one conclusion") found that once a screen has reduced
the problem to a few hundred parameters, ABC-SMC's sample-and-truncate
is the wrong instrument: a separable CMA-ES search over the same
objective reaches a distinctly better fit, several standard errors
better, with far tighter cross-seed reproducibility. This module is
the recommended path for tuning against experimental data --
:func:`~.tuning.bayesian_kcat_tuning` (ABC-SMC) still exists and is
unaffected by anything here, but is no longer what
``docs/bayesian_kcat_tuning.md`` walks through.

Three functions, normally called in this order:

- :func:`screen_kcat_leverage` -- report which kcats the data can
  actually speak to, with no optimisation involved. Perturbs each
  tie-group by +-``fold`` and measures the resulting change in RMSE,
  ranked by that leverage weighted by how much the source is trusted.
  Useful on its own for curation, independent of any tuning run.
- :func:`select_tunable_mask` -- turn a screen report into a boolean
  ``tunable_mask``, by keeping the fewest highest-ranked groups whose
  combined leverage reaches ``target_impact_share`` of the total. A
  *relative* cutoff, so the same ``target_impact_share`` selects a
  comparable quality of parameter set on a model with a different
  leverage scale, unlike an absolute threshold or a fixed count (see
  ``docs/internal/matlab_replication_results.md``, Open items #6).
- :func:`cmaes_kcat_tuning` -- the tuning run. Screens and selects
  automatically when no mask is given.
"""
from __future__ import annotations

import logging
from contextlib import nullcontext
from typing import TYPE_CHECKING, Callable, Optional

import cma
import cobra
import numpy as np
import pandas as pd
from cobra.util import ProcessPool

from ...ec_model.pipeline.apply_kcat import apply_kcat_constraints
from .data import BayesianData, load_bayesian_data
from .distance import compute_excarbon
from .priors import build_sigma0_log, classify_kcat_sources
from .tuning import (
    BayesianTuningResult,
    _init_worker,
    _score_kcat_vector,
    _score_worker,
    kcat_bounds,
)
from .tying import isozyme_tie_map

if TYPE_CHECKING:
    from ...adapter import ModelAdapter
    from ...adapter.params import BayesianParams
    from ...ec_model.ec_model import EcModel

logger = logging.getLogger(__name__)


def _resolve_context(model, adapter, params, bay_data, bio_rxn, okp_method):
    """Same resolution rules as ``tuning.bayesian_kcat_tuning``, so the
    two entry points behave identically when given the same inputs."""
    if adapter is None:
        adapter = getattr(model, "adapter", None)
    if params is None:
        if adapter is None:
            raise ValueError(
                "params not given and model has no adapter to read "
                "adapter.params.bayesian from."
            )
        params = adapter.params.bayesian
    if bay_data is None:
        if adapter is None:
            raise ValueError(
                "bay_data not given and model has no adapter to load it from."
            )
        bay_data = load_bayesian_data(adapter)
    if bay_data.flux_data is None and bay_data.max_grate is None:
        raise ValueError(
            "bay_data has neither flux_data nor max_grate -- nothing to "
            "tune against."
        )
    if bio_rxn is None:
        if adapter is None:
            raise ValueError("bio_rxn not given and model has no adapter.")
        bio_rxn = adapter.params.bio_rxn
    if okp_method is None and adapter is not None:
        okp_method = adapter.params.okp.method
    return params, bay_data, bio_rxn, okp_method


def _tunable_context(
    model, params, bay_data, bio_rxn, okp_method, tunable_mask,
):
    """The tunable subset plus everything derived from it: tie groups,
    trust weights, and the excarbon table. Shared by every function in
    this module so a screen and the run it feeds are always looking at
    exactly the same parameters."""
    is_tunable = model.ec.kcat > 0
    if tunable_mask is not None:
        tunable_mask = np.asarray(tunable_mask, dtype=bool)
        if tunable_mask.shape != is_tunable.shape:
            raise ValueError(
                f"tunable_mask has shape {tunable_mask.shape}; expected "
                f"{is_tunable.shape} to match model.ec.rxns."
            )
        is_tunable = is_tunable & tunable_mask
    tunable_idx = np.flatnonzero(is_tunable)
    if tunable_idx.size == 0:
        raise ValueError(
            "No tunable kcats: model.ec.kcat is all <= 0"
            + (", or tunable_mask excludes every one with a kcat."
               if tunable_mask is not None else ".")
        )
    ec_rxn_ids_tunable = [model.ec.rxns[i] for i in tunable_idx]
    kcat0 = model.ec.kcat[tunable_idx].astype(float).copy()
    sources = [model.ec.source[i] for i in tunable_idx]
    groups = classify_kcat_sources(sources, params, okp_method=okp_method)
    sigma0_log = build_sigma0_log(groups, params)

    excarbon_rxn_ids: set[str] = {bio_rxn}
    for data in (bay_data.flux_data, bay_data.max_grate):
        if data is not None:
            excarbon_rxn_ids.update(data.exch_rxn_ids)
    excarbon_rxn_ids.update(bay_data.zero_flux)
    excarbon = compute_excarbon(model, excarbon_rxn_ids, bio_rxn_id=bio_rxn)

    tie_map = (
        isozyme_tie_map(ec_rxn_ids_tunable, kcat0, list(groups))
        if params.tie_isozymes else np.arange(len(kcat0))
    )
    return (tunable_idx, ec_rxn_ids_tunable, kcat0, groups, sigma0_log,
            excarbon, tie_map)


def screen_kcat_leverage(
    model: "EcModel",
    *,
    adapter: Optional["ModelAdapter"] = None,
    params: Optional["BayesianParams"] = None,
    bay_data: Optional[BayesianData] = None,
    okp_method: Optional[str] = None,
    bio_rxn: Optional[str] = None,
    make_anaerobic: Optional[Callable[["EcModel"], None]] = None,
    change_protein_biomass: Optional[Callable[["EcModel", float], None]] = None,
    tunable_mask: Optional[np.ndarray] = None,
    fold: float = 2.0,
    n_proc: Optional[int] = None,
) -> pd.DataFrame:
    """Report which kcats the data can actually speak to.

    No optimisation happens here -- this only says which parameters
    are worth curating or tuning, and can be read on its own well
    before any tuning run. Each tie-group (isozyme copies sharing a
    prior and a source move together when ``params.tie_isozymes``) is
    perturbed up and down by ``fold``, and its leverage is the largest
    resulting change in plain RMSE -- always the plain fit, never the
    ``prior_penalty_weight``-penalised objective, since leverage asks
    what the data supports, independent of how much a move is charged
    for. Groups are ranked by that leverage weighted by
    ``sigma0_log`` -- :func:`~.parsimony.identifiability_mask`'s
    convention, where a trusted source needs proportionally more
    measured effect to rank alongside an untrusted one. This
    trust-weighting means the reported ``cum_leverage_share`` is a
    different quantity from :func:`~.parsimony.impact_share`, which is
    unweighted.

    Costs one simulation per condition per tie-group probed (two
    probes each, up and down) plus one baseline -- comparable in scale
    to a full tuning run's evaluation budget, not a cheap pre-check.
    FVA-unreachable kcats are not filtered out beforehand: a kcat that
    cannot carry flux under any condition shows zero leverage here by
    construction and simply sorts to the bottom, which is a more
    direct test than a separately computed static reachability mask.

    Returns
    -------
    pandas.DataFrame
        One row per tie-group, sorted by ``rank`` descending:
        ``rxn_id`` (the group's representative), ``n_isozymes``,
        ``source_group``, ``kcat0``, ``leverage``, ``sigma0_log``,
        ``rank`` (``leverage * sigma0_log``), and ``cum_leverage_share``
        (the running total of ``rank``, as a fraction of its sum over
        every row). An internal ``_positions`` column (a tuple of
        indices into ``model.ec.rxns``) is what
        :func:`select_tunable_mask` needs to expand a row back into a
        mask; drop it before showing the table to a person.
    """
    params, bay_data, bio_rxn, okp_method = _resolve_context(
        model, adapter, params, bay_data, bio_rxn, okp_method)
    (tunable_idx, ec_rxn_ids_tunable, kcat0, groups, sigma0_log, excarbon,
     tie_map) = _tunable_context(
        model, params, bay_data, bio_rxn, okp_method, tunable_mask)

    reps = np.flatnonzero(tie_map == np.arange(len(tie_map)))
    members_of = {int(r): np.flatnonzero(tie_map == r) for r in reps}

    if n_proc is None:
        n_proc = cobra.Configuration().processes
    n_proc = max(1, int(n_proc))

    pairs = []
    vectors = [kcat0]
    for rep in reps:
        members = members_of[int(rep)]
        up = kcat0.copy(); up[members] = up[members] * fold
        dn = kcat0.copy(); dn[members] = dn[members] / fold
        pairs.append(int(rep))
        vectors.append(up)
        vectors.append(dn)

    if n_proc == 1:
        rmses = [
            _score_kcat_vector(
                model, tunable_idx, ec_rxn_ids_tunable, bay_data, excarbon,
                bio_rxn, v, make_anaerobic=make_anaerobic,
                change_protein_biomass=change_protein_biomass,
                max_growth_weight=params.max_growth_weight,
            )[1]
            for v in vectors
        ]
    else:
        with ProcessPool(
            n_proc, initializer=_init_worker,
            initargs=(model, tunable_idx, ec_rxn_ids_tunable, bay_data,
                     excarbon, bio_rxn, make_anaerobic, change_protein_biomass,
                     params.max_growth_weight, 0.0, None, None),
        ) as pool:
            chunk = max(1, len(vectors) // (n_proc * 4))
            rmses = [r for _, r in pool.map(_score_worker, vectors, chunksize=chunk)]

    # The n_proc==1 path scores directly against `model`, leaving its
    # ec.kcat at whichever probe vector was scored last; a pool worker
    # scores its own copy instead, so this is a no-op there. Either
    # way, `model` must come out holding the prior it went in with.
    model.ec.kcat[tunable_idx] = kcat0
    apply_kcat_constraints(model, update_rxns=ec_rxn_ids_tunable)

    base_rmse = rmses[0]
    rows = []
    for i, rep in enumerate(pairs):
        up_rmse, dn_rmse = rmses[1 + 2 * i], rmses[2 + 2 * i]
        leverage = max(abs(up_rmse - base_rmse), abs(dn_rmse - base_rmse))
        members = members_of[rep]
        rows.append({
            "rxn_id": ec_rxn_ids_tunable[rep],
            "n_isozymes": int(len(members)),
            "source_group": str(groups[rep]),
            "kcat0": float(kcat0[rep]),
            "leverage": float(leverage),
            "sigma0_log": float(sigma0_log[rep]),
            "_positions": tuple(int(tunable_idx[m]) for m in members),
        })

    df = pd.DataFrame(rows)
    df["rank"] = df["leverage"] * df["sigma0_log"]
    df = df.sort_values("rank", ascending=False, ignore_index=True)
    total = float(df["rank"].sum())
    df["cum_leverage_share"] = (
        (df["rank"].cumsum() / total) if total > 0 else 0.0
    )
    return df


def select_tunable_mask(
    model: "EcModel",
    screen: pd.DataFrame,
    *,
    target_impact_share: float = 0.9,
) -> np.ndarray:
    """Build a ``tunable_mask`` from a :func:`screen_kcat_leverage` report.

    Keeps the fewest highest-ranked groups whose combined leverage
    reaches ``target_impact_share`` of the total -- a relative cutoff,
    so the same ``target_impact_share`` selects a comparable quality of
    parameter set on any model, rather than an absolute leverage value
    or a fixed count calibrated on a different model's scale.
    ``target_impact_share=0.9`` is a starting point, not a validated
    universal constant; check the resulting mask's size against
    :func:`screen_kcat_leverage`'s own curve before trusting it on a
    new model.

    Pure and cheap -- no simulation, only ``screen`` is read.
    """
    mask = np.zeros(len(model.ec.rxns), dtype=bool)
    if screen.empty:
        return mask
    prev_cum = screen["cum_leverage_share"].shift(1, fill_value=0.0)
    keep = prev_cum < target_impact_share
    for positions in screen.loc[keep, "_positions"]:
        mask[list(positions)] = True
    return mask


def cmaes_kcat_tuning(
    model: "EcModel",
    *,
    adapter: Optional["ModelAdapter"] = None,
    params: Optional["BayesianParams"] = None,
    bay_data: Optional[BayesianData] = None,
    okp_method: Optional[str] = None,
    bio_rxn: Optional[str] = None,
    make_anaerobic: Optional[Callable[["EcModel"], None]] = None,
    change_protein_biomass: Optional[Callable[["EcModel", float], None]] = None,
    tunable_mask: Optional[np.ndarray] = None,
    screen: Optional[pd.DataFrame] = None,
    target_impact_share: float = 0.9,
    popsize: Optional[int] = None,
    n_proc: Optional[int] = None,
    seed: Optional[int] = None,
    verbose: bool = True,
) -> BayesianTuningResult:
    """Tune kcats against experimental data with CMA-ES.

    The recommended entry point for fitting kcats to data -- see the
    module docstring for why this replaces ABC-SMC
    (:func:`~.tuning.bayesian_kcat_tuning`). Configuration comes
    entirely from ``BayesianParams``: trust tiers, ``max_growth_weight``,
    ``prior_penalty_weight`` and ``tie_isozymes`` mean exactly what they
    mean there, and ``max_generations``/``rmse_threshold`` double as
    this search's stopping conditions. Nothing new to learn if the
    ABC-SMC path is already familiar.

    Parameters
    ----------
    model
        EcModel with a populated ``ec.kcat``. Mutated in place: on
        return, tunable rows carry the best kcat vector CMA-ES found,
        with kcat constraints already applied.
    tunable_mask, screen, target_impact_share
        Control which kcats are searched over, in order of precedence.
        Pass ``tunable_mask`` to fix the set yourself. Otherwise, a
        mask is built by :func:`select_tunable_mask` at
        ``target_impact_share``, from ``screen`` if given (so a
        previously computed :func:`screen_kcat_leverage` report need
        not be recomputed) or from a fresh screen otherwise.
    popsize
        CMA-ES population size. Defaults to ``cma``'s own
        dimension-scaled default (``4 + floor(3 ln n)``), so it adapts
        to however many free parameters the mask/tying produced rather
        than a value calibrated on one particular model.
    n_proc, seed, verbose, make_anaerobic, change_protein_biomass,
    okp_method, bio_rxn
        As in :func:`~.tuning.bayesian_kcat_tuning`.

    Returns
    -------
    BayesianTuningResult
        ``rxns``/``old_kcat``/``new_kcat``/``groups`` cover the
        selected tunable set, tied groups included (their members share
        one value in ``new_kcat``). ``rmse_trace``/``objective_trace``
        record the best-so-far plain RMSE and optimised objective per
        generation. ``diagnostics_trace`` is always empty here -- CMA-ES
        has no per-generation accepted-particle set to compute
        source-group diagnostics over.

    Raises
    ------
    ValueError
        If there are no tunable kcats, or fewer than two free
        parameters remain after masking and tying -- too little for
        CMA-ES to search over.
    """
    params, bay_data, bio_rxn, okp_method = _resolve_context(
        model, adapter, params, bay_data, bio_rxn, okp_method)

    if tunable_mask is None:
        if screen is None:
            screen = screen_kcat_leverage(
                model, adapter=adapter, params=params, bay_data=bay_data,
                okp_method=okp_method, bio_rxn=bio_rxn,
                make_anaerobic=make_anaerobic,
                change_protein_biomass=change_protein_biomass,
                n_proc=n_proc,
            )
        tunable_mask = select_tunable_mask(
            model, screen, target_impact_share=target_impact_share)

    (tunable_idx, ec_rxn_ids_tunable, kcat0, groups, sigma0_log, excarbon,
     tie_map) = _tunable_context(
        model, params, bay_data, bio_rxn, okp_method, tunable_mask)

    reps = np.flatnonzero(tie_map == np.arange(len(tie_map)))
    if len(reps) < 2:
        raise ValueError(
            f"Only {len(reps)} free parameter(s) after masking and tying; "
            "too few for CMA-ES. Widen tunable_mask/target_impact_share."
        )
    assign = np.searchsorted(reps, tie_map)
    x0 = np.log(kcat0[reps])
    sig = sigma0_log[reps]
    lo, hi = kcat_bounds(kcat0)
    bounds = ([float(np.log(v)) for v in lo[reps]],
              [float(np.log(v)) for v in hi[reps]])

    def expand(xm: np.ndarray) -> np.ndarray:
        """(n_free, n_particles) log-space matrix -> (n_tunable,
        n_particles) kcat matrix, tied rows sharing their rep's value."""
        return np.exp(xm[assign, :])

    if n_proc is None:
        n_proc = cobra.Configuration().processes
    n_proc = max(1, int(n_proc))

    options = {
        "verbose": -9, "maxiter": params.max_generations,
        "CMA_diagonal": True, "bounds": bounds,
        "scaling_of_variables": (sig / np.median(sig)).tolist(),
    }
    if popsize is not None:
        options["popsize"] = int(popsize)
    if seed is not None:
        # cma treats seed=0 as falsy and falls back to an unseeded,
        # time-based draw instead of raising -- verified directly
        # against this cma version, not assumed. +1 keeps every
        # caller-facing seed (0, 1, 2, ...) off that value.
        options["seed"] = int(seed) + 1
    es = cma.CMAEvolutionStrategy(x0, float(np.median(sig)), options)

    pool_cm = (
        nullcontext(None) if n_proc == 1 else
        ProcessPool(
            n_proc, initializer=_init_worker,
            initargs=(model, tunable_idx, ec_rxn_ids_tunable, bay_data,
                     excarbon, bio_rxn, make_anaerobic, change_protein_biomass,
                     params.max_growth_weight, params.prior_penalty_weight,
                     kcat0, sigma0_log),
        )
    )

    def _score_batch(kcat_matrix, pool):
        if pool is None:
            return np.array([
                _score_kcat_vector(
                    model, tunable_idx, ec_rxn_ids_tunable, bay_data,
                    excarbon, bio_rxn, kcat_matrix[:, j],
                    make_anaerobic=make_anaerobic,
                    change_protein_biomass=change_protein_biomass,
                    max_growth_weight=params.max_growth_weight,
                    prior_penalty_weight=params.prior_penalty_weight,
                    kcat0=kcat0, sigma0_log=sigma0_log,
                )
                for j in range(kcat_matrix.shape[1])
            ])
        cols = [kcat_matrix[:, j] for j in range(kcat_matrix.shape[1])]
        chunk = max(1, len(cols) // (n_proc * 4))
        return np.array(pool.map(_score_worker, cols, chunksize=chunk))

    with pool_cm as pool:
        obj0, rmse0 = _score_batch(kcat0.reshape(-1, 1), pool)[0]
        best_obj, best_rmse, best_vec = float(obj0), float(rmse0), kcat0.copy()

        rmse_trace: list[float] = []
        objective_trace: list[float] = []
        generation = 0
        while not es.stop():
            asks = es.ask()
            xm = np.array(asks).T
            scores = _score_batch(expand(xm), pool)
            es.tell(asks, list(scores[:, 0]))
            generation += 1

            j = int(np.argmin(scores[:, 0]))
            if scores[j, 0] < best_obj:
                best_obj = float(scores[j, 0])
                best_rmse = float(scores[j, 1])
                best_vec = expand(xm[:, [j]])[:, 0].copy()
            rmse_trace.append(best_rmse)
            objective_trace.append(best_obj)
            if verbose:
                logger.info(
                    "generation %d: objective %.4f, rmse %.4f",
                    generation, best_obj, best_rmse,
                )
            if best_rmse <= params.rmse_threshold:
                break

    model.ec.kcat[tunable_idx] = best_vec
    apply_kcat_constraints(model, update_rxns=ec_rxn_ids_tunable)

    return BayesianTuningResult(
        rxns=ec_rxn_ids_tunable, old_kcat=kcat0.copy(), new_kcat=best_vec,
        groups=list(groups), rmse_trace=rmse_trace,
        objective_trace=objective_trace, diagnostics_trace=[],
        n_generations=generation,
        converged=best_rmse <= params.rmse_threshold,
    )
