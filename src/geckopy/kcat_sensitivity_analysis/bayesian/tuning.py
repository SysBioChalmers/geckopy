"""End-to-end Bayesian (ABC-SMC) kcat tuning.

Ported from GECKO MATLAB:
src/geckomat/kcat_sensitivity_analysis/Bayesian/bayesianSensitivityTuning.m,
using pyABC's ``Distribution``/``Transition`` classes as components
inside a custom generation loop this module owns (not
``pyabc.ABCSMC.run()`` -- see
``docs/internal/bayesian_tuning_plan.md``'s Architecture decision for
why: MATLAB's fixed-batch-then-truncate selection and pyABC's
streaming-until-N-accepted mechanism don't map onto each other without
fighting one or the other).

Parameterised by one swappable axis:

- ``selection="truncation"`` (MATLAB-faithful: combine this
  generation's new proposals with the previous generation's accepted
  set, keep the top ``min_keep`` fraction by distance) or
  ``"quantile_epsilon"`` (pyABC-native: draw a fresh batch each
  generation -- no carry-over -- and accept against a fixed,
  prospective epsilon derived from the *previous* generation).
Particles carry uniform weights, and per-source regularization is
carried entirely by the priors: each source group's ``sigma0_log``
sets how far a kcat can drift, and
``posterior.update_posterior_shrinkage``'s shrink-weight/
force-to-prior/sparsity-snap blend is computed for the trace only. It
does not feed back into sampling or into the returned model, matching
MATLAB: the final ``ecModel.ec.kcat`` is the best raw accepted
particle, and the next generation samples from the raw accepted
particle set rather than from the blended point estimate.

One MATLAB quirk is deliberately *not* ported: MATLAB drops any
proposal whose RMSE is *exactly* 0.0 ("often signals infeasibility").
This port's infeasibility is an explicit flag
(``simulate.ConditionSimResult.feasible``), never inferred from a
suspicious-looking distance value, so a genuinely perfect-fit proposal
(RMSE == 0) is kept rather than discarded.

Parallel scoring (``n_proc``)
------------------------------
Each generation's new particles are scored independently of each
other, so they parallelise across a process pool -- using
``cobra.util.process_pool.ProcessPool`` (the same primitive cobrapy's
own ``single_gene_deletion``/``single_reaction_deletion``/etc. use for
"repeated FBA on slightly-perturbed model copies", exactly this
module's shape of problem) rather than driving ``multiprocessing``
directly: one persistent ``EcModel`` copy per **worker process**
(deserialised once at pool startup, not once per particle and not one
shared mutable object across processes), which ``ProcessPool`` also
gets right on Windows for free (a temp-file-based initarg handoff that
works around a real performance issue in raw
``multiprocessing.Pool(initializer=...)`` there -- see its docstring).
Each worker then scores every particle it's handed against its own
copy via incremental ``apply_kcat_constraints(update_rxns=...)``
calls, exactly like the serial path -- see the "Spike results" section
of ``docs/internal/bayesian_tuning_plan.md`` for why a per-particle
``EcModel.copy()`` is not used either way (~225x an FBA solve on a
real-scale model).

Sampling itself (``build_kcat_prior``/``GeckoTransition.rvs_single``)
stays single-threaded in the main process -- it's cheap relative to
FBA and, more importantly, doing it there keeps the *sequence* of
particles proposed each generation deterministic given ``seed``,
independent of ``n_proc`` or how work happens to be chunked across
workers. Only the scoring of an already-fixed batch of particles is
parallelised, and scoring is a pure function of the kcat vector --
each particle's solves start from a cold basis
(:func:`_reset_solver_basis`), so a worker's result does not depend on
which particle it scored before. A run with ``n_proc=1`` and the same
``seed`` therefore reproduces bit-for-bit identical results to one
with ``n_proc>1``.

``make_anaerobic``/``change_protein_biomass``, if supplied, must be
importable top-level functions (not lambdas/closures) when running
with ``n_proc>1`` on Windows, where ``multiprocessing`` needs to
pickle them into each worker; POSIX's default start method (typically
``fork``, inherited via copy-on-write with no pickling involved) has
no such restriction.
"""
from __future__ import annotations

import contextlib
import logging
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Callable, Literal, Optional

import cobra
import numpy as np
import pandas as pd
from cobra.util import ProcessPool

from ...ec_model.pipeline.apply_kcat import apply_kcat_constraints
from .data import BayesianData, load_bayesian_data
from .diagnostics import GenerationDiagnostics, compute_generation_diagnostics
from .distance import bayesian_distance, compute_excarbon
from .posterior import PosteriorUpdate, update_posterior_shrinkage
from .priors import (build_kcat_prior, build_proposal_sigma_log,
                     build_sigma0_log, classify_kcat_sources)
from .selection import next_quantile_epsilon, quantile_epsilon_select, truncation_select
from .simulate import simulate_bayesian_dataset
from .transition import GeckoTransition

if TYPE_CHECKING:
    from ...adapter import ModelAdapter
    from ...adapter.params import BayesianParams
    from ...ec_model.ec_model import EcModel

logger = logging.getLogger(__name__)

SelectionVariant = Literal["truncation", "quantile_epsilon"]


@dataclass
class BayesianTuningResult:
    """Outcome of a :func:`bayesian_kcat_tuning` run.

    Matches the ``TunedKcatsResult``/``SigmaFitterResult`` precedent:
    the model is mutated in place (the final generation's best-RMSE
    particle is applied to ``model.ec.kcat``), and this result records
    the run's history for inspection.

    Attributes
    ----------
    rxns
        Tunable ``ec.rxns`` IDs, in the column order used throughout
        (``rxns[i]`` corresponds to ``old_kcat[i]``/``new_kcat[i]``/
        ``groups[i]``).
    old_kcat, new_kcat
        kcat values before/after tuning, parallel to ``rxns``.
    groups
        Source-group name per tunable row (from
        ``priors.classify_kcat_sources``).
    rmse_trace
        Best (lowest) RMSE in the accepted set, per generation. Always
        the plain RMSE, so runs stay comparable across penalties and
        against MATLAB even when selection used a penalised objective.
    objective_trace
        Best (lowest) value of the quantity selection actually
        minimised, per generation. Identical to ``rmse_trace`` when
        ``params.prior_penalty_weight`` is 0.
    posterior_trace
        ``posterior.PosteriorUpdate`` per generation (see the module
        docstring: recorded for inspection, never fed back into
        sampling or the returned model).
    diagnostics_trace
        Per-generation, per-source-group diagnostics (see
        ``diagnostics.py``).
    n_generations
        Number of generations actually run.
    converged
        True if ``rmse_trace[-1] <= params.rmse_threshold``; False if
        ``max_generations`` was reached first.
    """

    rxns: list[str] = field(default_factory=list)
    old_kcat: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=float))
    new_kcat: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=float))
    groups: list[str] = field(default_factory=list)
    rmse_trace: list[float] = field(default_factory=list)
    objective_trace: list[float] = field(default_factory=list)
    posterior_trace: list[PosteriorUpdate] = field(default_factory=list)
    diagnostics_trace: list[GenerationDiagnostics] = field(default_factory=list)
    n_generations: int = 0
    converged: bool = False


def bayesian_kcat_tuning(
    model: "EcModel",
    *,
    adapter: Optional["ModelAdapter"] = None,
    params: Optional["BayesianParams"] = None,
    bay_data: Optional[BayesianData] = None,
    selection: SelectionVariant = "truncation",
    okp_method: Optional[str] = None,
    bio_rxn: Optional[str] = None,
    make_anaerobic: Optional[Callable[["EcModel"], None]] = None,
    change_protein_biomass: Optional[Callable[["EcModel", float], None]] = None,
    tunable_mask: Optional[np.ndarray] = None,
    n_proc: Optional[int] = None,
    seed: Optional[int] = None,
    verbose: bool = True,
) -> BayesianTuningResult:
    """Run ABC-SMC Bayesian kcat tuning.

    Parameters
    ----------
    model
        EcModel with a populated ``ec.kcat``. Mutated in place: on
        return, tunable rows carry the final generation's best-RMSE
        particle.
    adapter
        Used to resolve ``params``/``bay_data``/``bio_rxn``/
        ``okp_method`` when not given explicitly. Defaults to
        ``model.adapter``.
    params
        Hyperparameters. Defaults to ``adapter.params.bayesian``.
    bay_data
        Experimental data. Defaults to ``load_bayesian_data(adapter)``.
    selection
        Selection variant (see module docstring).
    okp_method
        The project's configured OpenKineticsPredictor method, for
        source classification. Defaults to ``adapter.params.okp.method``
        if an adapter is available.
    bio_rxn
        Biomass reaction ID. Defaults to ``adapter.params.bio_rxn``.
    make_anaerobic, change_protein_biomass
        Forwarded to ``simulate.simulate_bayesian_dataset`` -- see its
        docstring; geckopy has no generic organism-agnostic
        implementation of these yet. See the module docstring's
        "Parallel scoring" section for a picklability caveat when
        ``n_proc>1``.
    tunable_mask
        Boolean array over ``model.ec.rxns`` selecting which kcats to
        tune; the rest are held at their prior values and reported as
        unchanged. Use it to exclude parameters the data cannot speak
        to -- reactions that carry no flux in any condition, or whose
        kcat does not measurably move the distance. On the full
        ecYeastGEM dataset 879 of 4834 kcats can never carry flux, and
        an unrestricted run changes essentially all of them, which no
        weighting of the objective can justify. Restricting also
        shrinks the proposal's step length, which grows with the square
        root of the dimension.
    n_proc
        Number of worker processes for scoring each generation's
        particles. Defaults to ``cobra.Configuration().processes``.
        ``1`` runs the original serial path (no ``Pool`` at all). See
        the module docstring's "Parallel scoring" section.
    seed
        If given, seeds ``numpy.random`` before sampling starts, for
        reproducible runs.
    verbose
        Whether per-generation progress is logged at INFO.

    Returns
    -------
    BayesianTuningResult

    Raises
    ------
    ValueError
        If there are no tunable kcats (``ec.kcat`` all ``<= 0``), or
        ``bay_data`` has neither ``flux_data`` nor ``max_grate``.
    """
    if seed is not None:
        np.random.seed(seed)

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
    proposal_sigma_log = build_proposal_sigma_log(groups, params)
    n_params = len(kcat0)
    columns = [f"k{i}" for i in range(n_params)]

    excarbon_rxn_ids: set[str] = {bio_rxn}
    for data in (bay_data.flux_data, bay_data.max_grate):
        if data is not None:
            excarbon_rxn_ids.update(data.exch_rxn_ids)
    excarbon_rxn_ids.update(bay_data.zero_flux)
    excarbon = compute_excarbon(model, excarbon_rxn_ids, bio_rxn_id=bio_rxn)
    kcat_lo, kcat_hi = kcat_bounds(kcat0)

    if n_proc is None:
        n_proc = cobra.Configuration().processes
    n_proc = max(1, int(n_proc))

    if n_proc == 1:
        pool_cm = contextlib.nullcontext(None)
    else:
        # ProcessPool (cobra.util.process_pool) handles serialising the
        # model to each worker -- including a Windows-specific
        # performance workaround -- so `model` is passed as-is, not
        # pre-pickled.
        pool_cm = ProcessPool(
            n_proc, initializer=_init_worker,
            initargs=(
                model, tunable_idx, ec_rxn_ids_tunable, bay_data,
                excarbon, bio_rxn, make_anaerobic, change_protein_biomass,
                params.max_growth_weight, params.prior_penalty_weight,
                kcat0, sigma0_log,
            ),
        )

    with pool_cm as pool:
        if pool is None:
            def score_batch(kcat_matrix: np.ndarray) -> np.ndarray:
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
        else:
            def score_batch(kcat_matrix: np.ndarray) -> np.ndarray:
                cols = [kcat_matrix[:, j] for j in range(kcat_matrix.shape[1])]
                chunk = max(1, len(cols) // (n_proc * 4))
                return np.array(pool.map(_score_worker, cols, chunksize=chunk))

        # Seed the accepted pool with the model's own starting point --
        # MATLAB: `kcats = ecModel.ec.kcat; kcat0 = kcats; rmse =
        # abc_max(...); rmseTop = rmse; kcatTop = kcats;`.
        kcat_top = kcat0.reshape(-1, 1)
        # score_batch returns one (objective, rmse) row per particle:
        # selection runs on the objective, reporting on the plain RMSE.
        seed_scores = score_batch(kcat_top)
        obj_top, rmse_top = seed_scores[:, 0], seed_scores[:, 1]
        weights_top = np.array([1.0])
        epsilon: Optional[float] = None

        result = BayesianTuningResult(
            rxns=ec_rxn_ids_tunable, old_kcat=kcat0.copy(), groups=list(groups),
        )

        # Persistent multiplier on the fitted proposal bandwidth. Stays
        # at 1.0 (a no-op) unless params.adapt_proposal_width is set.
        proposal_scale = 1.0

        generation = 0
        while True:
            generation += 1
            schedule_idx = [
                i for i, sg in enumerate(params.schedule_generations)
                if generation >= sg
            ]
            n_new = params.schedule_samples[schedule_idx[-1] if schedule_idx else 0]

            # One draw per particle, read every column out of that single
            # draw: both samplers return a full parameter vector per call,
            # and the transition kernel's vector is a single parent
            # perturbed as a unit. Cost is linear in particles, not in
            # particles x parameters.
            if generation == 1:
                prior = build_kcat_prior(kcat0, sigma0_log)
                draws = [prior.rvs() for _ in range(n_new)]
                new_particles = np.array(
                    [[draw[c] for c in columns] for draw in draws]
                ).T
                transition: Optional[GeckoTransition] = None
            else:
                X_df = pd.DataFrame(kcat_top.T, columns=columns)
                transition = GeckoTransition(
                    proposal_sigma_log, bandwidth_scale=proposal_scale,
                )
                transition.fit(X_df, weights_top)
                new_particles = transition.rvs_batch(n_new)
            new_particles = np.clip(new_particles, kcat_lo[:, None], kcat_hi[:, None])
            new_scores = score_batch(new_particles)
            new_obj, new_rmse = new_scores[:, 0], new_scores[:, 1]

            if selection == "truncation":
                combined_particles = np.concatenate([new_particles, kcat_top], axis=1)
                combined_obj = np.concatenate([new_obj, obj_top])
                combined_rmse = np.concatenate([new_rmse, rmse_top])
                sel = truncation_select(combined_obj, min_keep=params.min_keep)
            else:
                combined_particles = new_particles
                combined_obj = new_obj
                combined_rmse = new_rmse
                if epsilon is None:  # generation 1: bootstrap from this batch
                    epsilon = next_quantile_epsilon(combined_obj)
                sel = quantile_epsilon_select(combined_obj, epsilon=epsilon)
                if sel.accepted_idx.size == 0:
                    sel.accepted_idx = np.array([int(np.argmin(combined_obj))])

            # Fraction of *this* generation's proposals that survived
            # selection: how many of the n_new draws got in, over n_new.
            # Columns 0..n_new-1 of combined_particles are the new draws
            # and the rest are the carried-over accepted set, so the
            # numerator counts accepted indices below that boundary.
            # Note this is measured at the min_keep threshold, unlike
            # MATLAB's propAccRate, which is measured at its
            # targetAccept percentile before the minKeep clamp -- the
            # two are not on the same scale.
            proposal_accept_rate = proposal_acceptance_rate(
                sel.accepted_idx, new_particles.shape[1]
            )
            if params.adapt_proposal_width:
                proposal_scale = adapt_proposal_scale(
                    proposal_scale, proposal_accept_rate, params,
                )

            n_total_this_gen = len(combined_rmse)
            new_kcat_top = combined_particles[:, sel.accepted_idx]
            new_rmse_top = combined_rmse[sel.accepted_idx]
            new_obj_top = combined_obj[sel.accepted_idx]

            new_weights_top = np.full(
                len(sel.accepted_idx), 1.0 / len(sel.accepted_idx),
            )

            posterior_update = update_posterior_shrinkage(
                new_kcat_top, kcat0, sigma0_log, groups,
                shrink_thr_default=params.shrink_thr_default,
                shrink_thr_source=params.shrink_thr_source,
                force_prior_thr_default=params.force_prior_thr_default,
                force_prior_thr_source=params.force_prior_thr_source,
                sparsity_threshold=params.sparsity_threshold,
            )
            result.posterior_trace.append(posterior_update)

            diag = compute_generation_diagnostics(
                generation, new_kcat_top, new_weights_top, new_rmse_top,
                n_total_this_gen, kcat0, sigma0_log, groups,
            )
            result.diagnostics_trace.append(diag)
            result.rmse_trace.append(diag.best_rmse)

            result.objective_trace.append(float(new_obj_top.min()))

            if selection == "quantile_epsilon":
                epsilon = next_quantile_epsilon(new_obj_top)

            if verbose:
                logger.info(
                    "bayesian_kcat_tuning: generation %d, best RMSE = %g, "
                    "accepted %d/%d, proposal accept %.4f, scale %.4f",
                    generation, diag.best_rmse, diag.n_accepted, diag.n_total,
                    proposal_accept_rate, proposal_scale,
                )

            kcat_top, rmse_top, weights_top = new_kcat_top, new_rmse_top, new_weights_top
            obj_top = new_obj_top

            if diag.best_rmse <= params.rmse_threshold:
                result.converged = True
                break
            if generation >= params.max_generations:
                result.converged = False
                break

    result.n_generations = generation
    # Consistent with what selection optimised: at prior_penalty_weight
    # 0 the objective is the RMSE and this is unchanged.
    best_idx = int(np.argmin(obj_top))
    result.new_kcat = kcat_top[:, best_idx].copy()
    model.ec.kcat[tunable_idx] = result.new_kcat
    apply_kcat_constraints(model, update_rxns=ec_rxn_ids_tunable)

    return result


# --------------------------------------------------------------------------- #
# Scoring core, shared by the serial and parallel paths
# --------------------------------------------------------------------------- #

def kcat_bounds(kcat0: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Biologically plausible bounds for proposed kcats, in 1/s.

    Proposals outside these bounds are not hypotheses worth spending an
    FBA solve on. The window is 1e-2 to 1e4 for an ordinary kcat, and
    always widened far enough to contain the prior with a hundred-fold
    margin either side, so an unusually slow or fast enzyme keeps a
    window around itself rather than being clipped towards the generic
    range.

    That last part matters more than it sounds. MATLAB's
    ``proposeSimple`` widens only upwards, which leaves any prior below
    1e-2 outside its own bounds: on ecYeastGEM that is 350 kcats, and
    every proposal for them is clipped up to 1e-2 before it is scored --
    a 62-fold increase for the slowest, applied before any evidence is
    considered. A prior must always be a proposable value.

    Parameters
    ----------
    kcat0
        Prior kcat per tunable row, 1/s.

    Returns
    -------
    tuple of numpy.ndarray
        Lower and upper bound per row.
    """
    kcat0 = np.asarray(kcat0, dtype=float)
    lo = np.minimum(1e-2, kcat0 / 100.0)
    hi = np.maximum(1e4, kcat0 * 100.0)
    return lo, hi


def _reset_solver_basis(model: "EcModel") -> None:
    """Discard the solver's incumbent basis so the solves below start cold.

    The distance must be a function of the kcat vector alone. Solvers
    warm-start from the previous solve's basis, and these LPs have
    alternate optima: resuming from a basis that is already optimal
    for the new problem returns that vertex in zero iterations, so the
    reported exchange fluxes -- and hence the RMSE -- would otherwise
    depend on which particle this model scored before. That would make
    the accepted set depend on ``n_proc`` and on how the pool happened
    to schedule particles across workers.

    A no-op for solver interfaces whose problem object exposes no
    ``reset`` (only the model's own state then determines the solve,
    which is the property this guarantees for the rest).
    """
    problem = getattr(model.solver, "problem", None)
    reset = getattr(problem, "reset", None)
    if callable(reset):
        reset()


def proposal_acceptance_rate(accepted_idx: np.ndarray, n_new: int) -> float:
    """Fraction of a generation's *new* proposals that survived selection.

    The selection step ranks ``[new_particles, kcat_top]`` as one pool,
    so accepted indices below ``n_new`` are the new draws and the rest
    are carried-over particles. The denominator is ``n_new`` -- the
    number of proposals made -- not the size of the accepted set, which
    would instead measure the accepted set's composition and sits near
    1.0 in early generations whatever the proposals are worth.

    Measured at the ``min_keep`` truncation threshold. MATLAB's
    ``propAccRate`` is measured at its ``targetAccept`` percentile
    before the ``minKeep`` clamp, so the two are not on the same scale
    and their targets are not interchangeable.
    """
    if n_new <= 0:
        return 0.0
    return float(np.count_nonzero(np.asarray(accepted_idx) < n_new) / n_new)


def adapt_proposal_scale(
    scale: float, accept_rate: float, params: "BayesianParams",
) -> float:
    """Steer the proposal bandwidth by how many proposals survive.

    MATLAB's proposal width is ``0.5 * std_obs + 0.5 * sigma0_log``,
    which cannot fall below half the prior width and grows as the
    accepted set spreads. Its own diagnostics show the consequence: on
    the full ecYeastGEM run the width climbs 0.213 -> 1.037 while the
    proposal acceptance rate collapses 0.10 -> 0.005, so the sampler
    widens exactly when its steps are already too large and the last
    third of the run makes no progress.

    This multiplies the fitted bandwidth by a scale that follows the
    acceptance rate instead, in the usual adaptive-MCMC form:
    ``scale *= exp(rate_param * (accept_rate - target))``, clamped to
    ``params.proposal_scale_bounds``. Below target the steps shrink;
    above it they grow back.

    MATLAB has no such feedback, so this is a deliberate departure from
    the faithful path and is off unless ``params.adapt_proposal_width``
    is set.
    """
    lo, hi = params.proposal_scale_bounds
    adjusted = scale * float(np.exp(
        params.proposal_adaptation_rate * (accept_rate - params.target_accept_rate)
    ))
    return float(np.clip(adjusted, lo, hi))


def _score_kcat_vector(
    model: "EcModel",
    tunable_idx: np.ndarray,
    ec_rxn_ids_tunable: list[str],
    bay_data: BayesianData,
    excarbon: dict[str, float],
    bio_rxn_id: str,
    kcat_vec: np.ndarray,
    *,
    make_anaerobic,
    change_protein_biomass,
    max_growth_weight: float = 1.0,
    prior_penalty_weight: float = 0.0,
    kcat0: Optional[np.ndarray] = None,
    sigma0_log: Optional[np.ndarray] = None,
) -> tuple[float, float]:
    """Score one kcat vector against ``model`` (mutated in place).

    Writes the candidate kcat vector into ``model.ec.kcat`` and
    rewrites just the tunable rows' stoichiometry via
    ``apply_kcat_constraints(update_rxns=...)`` -- no per-particle
    ``EcModel.copy()`` (prohibitively expensive; see the "Spike
    results" section of ``docs/internal/bayesian_tuning_plan.md``).
    Called directly (looped in-process) for the serial path, or via
    :func:`_score_worker` from inside a pool worker for the parallel
    path -- either way, ``model`` is one persistent object reused
    across every particle it's asked to score.
    """
    model.ec.kcat[tunable_idx] = kcat_vec
    apply_kcat_constraints(model, update_rxns=ec_rxn_ids_tunable)
    _reset_solver_basis(model)

    flux_sims = None
    if bay_data.flux_data is not None:
        flux_sims = simulate_bayesian_dataset(
            model, bay_data.flux_data,
            constrain=True, zero_flux_rxns=bay_data.zero_flux,
            bio_rxn_id=bio_rxn_id,
            make_anaerobic=make_anaerobic,
            change_protein_biomass=change_protein_biomass,
        )
    max_grate_sims = None
    if bay_data.max_grate is not None:
        max_grate_sims = simulate_bayesian_dataset(
            model, bay_data.max_grate,
            constrain=False, zero_flux_rxns=[],
            bio_rxn_id=bio_rxn_id,
            make_anaerobic=make_anaerobic,
            change_protein_biomass=change_protein_biomass,
        )
    rmse, _ = bayesian_distance(
        bay_data, flux_sims=flux_sims, max_grate_sims=max_grate_sims,
        excarbon=excarbon, bio_rxn_id=bio_rxn_id,
        max_growth_weight=max_growth_weight,
    )
    # Selection may run on a penalised objective, but the plain RMSE is
    # always carried alongside it so a run stays comparable to MATLAB's
    # and to runs at other penalties.
    objective = rmse
    if prior_penalty_weight and kcat0 is not None and sigma0_log is not None:
        dev = np.log(np.asarray(kcat_vec, dtype=float) / kcat0) / sigma0_log
        objective = rmse + prior_penalty_weight * float(np.mean(dev ** 2))
    return objective, rmse


# --------------------------------------------------------------------------- #
# Parallel scoring (n_proc > 1): worker-process globals, populated once
# by _init_worker and reused for every particle that worker scores.
# --------------------------------------------------------------------------- #

_WORKER_MODEL: Optional["EcModel"] = None
_WORKER_TUNABLE_IDX: Optional[np.ndarray] = None
_WORKER_EC_RXN_IDS_TUNABLE: Optional[list[str]] = None
_WORKER_BAY_DATA: Optional[BayesianData] = None
_WORKER_EXCARBON: Optional[dict[str, float]] = None
_WORKER_BIO_RXN: Optional[str] = None
_WORKER_MAKE_ANAEROBIC = None
_WORKER_CHANGE_PROTEIN_BIOMASS = None
_WORKER_MAX_GROWTH_WEIGHT: float = 1.0
_WORKER_PRIOR_PENALTY_WEIGHT: float = 0.0
_WORKER_KCAT0: Optional[np.ndarray] = None
_WORKER_SIGMA0_LOG: Optional[np.ndarray] = None


def _init_worker(
    model: "EcModel",
    tunable_idx: np.ndarray,
    ec_rxn_ids_tunable: list[str],
    bay_data: BayesianData,
    excarbon: dict[str, float],
    bio_rxn_id: str,
    make_anaerobic,
    change_protein_biomass,
    max_growth_weight: float = 1.0,
    prior_penalty_weight: float = 0.0,
    kcat0: Optional[np.ndarray] = None,
    sigma0_log: Optional[np.ndarray] = None,
) -> None:
    """Pool initializer: stash this worker process's own EcModel copy
    (deserialised by ``ProcessPool``, not by us) and everything else
    needed to score a particle."""
    global _WORKER_MODEL, _WORKER_TUNABLE_IDX, _WORKER_EC_RXN_IDS_TUNABLE
    global _WORKER_BAY_DATA, _WORKER_EXCARBON, _WORKER_BIO_RXN
    global _WORKER_MAKE_ANAEROBIC, _WORKER_CHANGE_PROTEIN_BIOMASS
    global _WORKER_MAX_GROWTH_WEIGHT, _WORKER_PRIOR_PENALTY_WEIGHT
    global _WORKER_KCAT0, _WORKER_SIGMA0_LOG
    _WORKER_MODEL = model
    _WORKER_TUNABLE_IDX = tunable_idx
    _WORKER_EC_RXN_IDS_TUNABLE = ec_rxn_ids_tunable
    _WORKER_BAY_DATA = bay_data
    _WORKER_EXCARBON = excarbon
    _WORKER_BIO_RXN = bio_rxn_id
    _WORKER_MAKE_ANAEROBIC = make_anaerobic
    _WORKER_CHANGE_PROTEIN_BIOMASS = change_protein_biomass
    _WORKER_MAX_GROWTH_WEIGHT = max_growth_weight
    _WORKER_PRIOR_PENALTY_WEIGHT = prior_penalty_weight
    _WORKER_KCAT0 = kcat0
    _WORKER_SIGMA0_LOG = sigma0_log


def _score_worker(kcat_vec: np.ndarray) -> tuple[float, float]:
    assert _WORKER_MODEL is not None, "_score_worker called before _init_worker"
    return _score_kcat_vector(
        _WORKER_MODEL, _WORKER_TUNABLE_IDX, _WORKER_EC_RXN_IDS_TUNABLE,
        _WORKER_BAY_DATA, _WORKER_EXCARBON, _WORKER_BIO_RXN, kcat_vec,
        make_anaerobic=_WORKER_MAKE_ANAEROBIC,
        change_protein_biomass=_WORKER_CHANGE_PROTEIN_BIOMASS,
        max_growth_weight=_WORKER_MAX_GROWTH_WEIGHT,
        prior_penalty_weight=_WORKER_PRIOR_PENALTY_WEIGHT,
        kcat0=_WORKER_KCAT0, sigma0_log=_WORKER_SIGMA0_LOG,
    )
