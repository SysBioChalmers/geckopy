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

Parameterised by two independently swappable axes:

- ``selection="truncation"`` (MATLAB-faithful: combine this
  generation's new proposals with the previous generation's accepted
  set, keep the top ``min_keep`` fraction by distance) or
  ``"quantile_epsilon"`` (pyABC-native: draw a fresh batch each
  generation -- no carry-over -- and accept against a fixed,
  prospective epsilon derived from the *previous* generation).
- ``regularization="shrinkage"`` (MATLAB-faithful: uniform particle
  weights; ``posterior.update_posterior_shrinkage``'s shrink-weight/
  force-to-prior/sparsity-snap blend is computed purely for the trace
  -- tracing MATLAB's own data flow shows that computation never
  feeds back into sampling *or* the returned model there either: the
  final ``ecModel.ec.kcat`` comes from the best raw accepted particle,
  and the next generation samples from the raw accepted particle set,
  not from the blended point estimate) or ``"importance_weighting"``
  (proper SMC-ABC: ``importance_weights.compute_importance_weights``'
  prior/transition-density weights feed both the next generation's
  resampling, via ``GeckoTransition``, and the weighted diagnostics).

One MATLAB quirk is deliberately *not* ported: MATLAB drops any
proposal whose RMSE is *exactly* 0.0 ("often signals infeasibility").
This port's infeasibility is an explicit flag
(``simulate.ConditionSimResult.feasible``), never inferred from a
suspicious-looking distance value, so a genuinely perfect-fit proposal
(RMSE == 0) is kept rather than discarded.
"""
from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Callable, Literal, Optional

import numpy as np
import pandas as pd

from ...ec_model.pipeline.apply_kcat import apply_kcat_constraints
from .data import BayesianData, load_bayesian_data
from .diagnostics import GenerationDiagnostics, compute_generation_diagnostics
from .distance import bayesian_distance, compute_excarbon
from .importance_weights import compute_importance_weights
from .posterior import PosteriorUpdate, update_posterior_shrinkage
from .priors import build_kcat_prior, build_sigma0_log, classify_kcat_sources, kcat_prior_logpdf
from .selection import next_quantile_epsilon, quantile_epsilon_select, truncation_select
from .simulate import simulate_bayesian_dataset
from .transition import GeckoTransition

if TYPE_CHECKING:
    from ...adapter import ModelAdapter
    from ...adapter.params import BayesianParams
    from ...ec_model.ec_model import EcModel

logger = logging.getLogger(__name__)

SelectionVariant = Literal["truncation", "quantile_epsilon"]
RegularizationVariant = Literal["shrinkage", "importance_weighting"]


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
        Best (lowest) RMSE in the accepted set, per generation.
    posterior_trace
        ``posterior.PosteriorUpdate`` per generation -- only populated
        when ``regularization="shrinkage"`` (see module docstring: for
        ``"importance_weighting"`` there is no equivalent blended
        point estimate, only the particle population + weights
        already captured in ``diagnostics_trace``).
    diagnostics_trace
        Per-generation, per-source-group diagnostics -- populated for
        both regularization variants, and directly comparable between
        them (see ``diagnostics.py``).
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
    regularization: RegularizationVariant = "shrinkage",
    okp_method: Optional[str] = None,
    bio_rxn: Optional[str] = None,
    make_anaerobic: Optional[Callable[["EcModel"], None]] = None,
    change_protein_biomass: Optional[Callable[["EcModel", float], None]] = None,
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
        Axis 1 variant (see module docstring).
    regularization
        Axis 2 variant (see module docstring).
    okp_method
        The project's configured OpenKineticsPredictor method, for
        source classification. Defaults to ``adapter.params.okp.method``
        if an adapter is available.
    bio_rxn
        Biomass reaction ID. Defaults to ``adapter.params.bio_rxn``.
    make_anaerobic, change_protein_biomass
        Forwarded to ``simulate.simulate_bayesian_dataset`` -- see its
        docstring; geckopy has no generic organism-agnostic
        implementation of these yet.
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

    tunable_idx = np.flatnonzero(model.ec.kcat > 0)
    if tunable_idx.size == 0:
        raise ValueError("No tunable kcats: model.ec.kcat is all <= 0.")
    ec_rxn_ids_tunable = [model.ec.rxns[i] for i in tunable_idx]
    kcat0 = model.ec.kcat[tunable_idx].astype(float).copy()
    sources = [model.ec.source[i] for i in tunable_idx]
    groups = classify_kcat_sources(sources, params, okp_method=okp_method)
    sigma0_log = build_sigma0_log(groups, params)
    n_params = len(kcat0)
    columns = [f"k{i}" for i in range(n_params)]

    excarbon_rxn_ids: set[str] = {bio_rxn}
    for data in (bay_data.flux_data, bay_data.max_grate):
        if data is not None:
            excarbon_rxn_ids.update(data.exch_rxn_ids)
    excarbon_rxn_ids.update(bay_data.zero_flux)
    excarbon = compute_excarbon(model, excarbon_rxn_ids, bio_rxn_id=bio_rxn)

    score = _make_scorer(
        model, tunable_idx, ec_rxn_ids_tunable, bay_data, excarbon, bio_rxn,
        make_anaerobic=make_anaerobic, change_protein_biomass=change_protein_biomass,
    )

    # Seed the accepted pool with the model's own starting point --
    # MATLAB: `kcats = ecModel.ec.kcat; kcat0 = kcats; rmse =
    # abc_max(...); rmseTop = rmse; kcatTop = kcats;`.
    kcat_top = kcat0.reshape(-1, 1)
    rmse_top = np.array([score(kcat0.copy())])
    weights_top = np.array([1.0])
    epsilon: Optional[float] = None
    prior_logpdf = lambda theta: kcat_prior_logpdf(theta, kcat0, sigma0_log)  # noqa: E731

    result = BayesianTuningResult(
        rxns=ec_rxn_ids_tunable, old_kcat=kcat0.copy(), groups=list(groups),
    )

    generation = 0
    while True:
        generation += 1
        schedule_idx = [
            i for i, sg in enumerate(params.schedule_generations) if generation >= sg
        ]
        n_new = params.schedule_samples[schedule_idx[-1] if schedule_idx else 0]

        if generation == 1:
            prior = build_kcat_prior(kcat0, sigma0_log)
            new_particles = np.array(
                [[prior.rvs()[c] for c in columns] for _ in range(n_new)]
            ).T
            transition: Optional[GeckoTransition] = None
        else:
            X_df = pd.DataFrame(kcat_top.T, columns=columns)
            transition = GeckoTransition(sigma0_log)
            transition.fit(X_df, weights_top)
            new_particles = np.array(
                [[transition.rvs_single()[c] for c in columns] for _ in range(n_new)]
            ).T
        new_particles = np.maximum(new_particles, np.finfo(float).tiny)
        new_rmse = np.array([score(new_particles[:, j]) for j in range(n_new)])

        if selection == "truncation":
            combined_particles = np.concatenate([new_particles, kcat_top], axis=1)
            combined_rmse = np.concatenate([new_rmse, rmse_top])
            sel = truncation_select(combined_rmse, min_keep=params.min_keep)
        else:
            combined_particles = new_particles
            combined_rmse = new_rmse
            if epsilon is None:  # generation 1: bootstrap from this batch
                epsilon = next_quantile_epsilon(combined_rmse)
            sel = quantile_epsilon_select(combined_rmse, epsilon=epsilon)
            if sel.accepted_idx.size == 0:
                sel.accepted_idx = np.array([int(np.argmin(combined_rmse))])

        n_total_this_gen = len(combined_rmse)
        new_kcat_top = combined_particles[:, sel.accepted_idx]
        new_rmse_top = combined_rmse[sel.accepted_idx]

        if regularization == "shrinkage" or transition is None:
            new_weights_top = np.full(
                len(sel.accepted_idx), 1.0 / len(sel.accepted_idx),
            )
        else:
            new_weights_top = compute_importance_weights(
                new_kcat_top, prior_logpdf,
                parents=kcat_top, parent_weights=weights_top,
                transition_logpdf=transition.component_logpdf,
            )

        if regularization == "shrinkage":
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

        if selection == "quantile_epsilon":
            epsilon = next_quantile_epsilon(
                new_rmse_top,
                weights=new_weights_top if regularization == "importance_weighting" else None,
            )

        if verbose:
            logger.info(
                "bayesian_kcat_tuning: generation %d, best RMSE = %g, "
                "accepted %d/%d",
                generation, diag.best_rmse, diag.n_accepted, diag.n_total,
            )

        kcat_top, rmse_top, weights_top = new_kcat_top, new_rmse_top, new_weights_top

        if diag.best_rmse <= params.rmse_threshold:
            result.converged = True
            break
        if generation >= params.max_generations:
            result.converged = False
            break

    result.n_generations = generation
    best_idx = int(np.argmin(rmse_top))
    result.new_kcat = kcat_top[:, best_idx].copy()
    model.ec.kcat[tunable_idx] = result.new_kcat
    apply_kcat_constraints(model, update_rxns=ec_rxn_ids_tunable)

    return result


def _make_scorer(
    model: "EcModel",
    tunable_idx: np.ndarray,
    ec_rxn_ids_tunable: list[str],
    bay_data: BayesianData,
    excarbon: dict[str, float],
    bio_rxn_id: str,
    *,
    make_anaerobic,
    change_protein_biomass,
) -> Callable[[np.ndarray], float]:
    """One persistent-model scorer, reused for every particle.

    Writes the candidate kcat vector into ``model.ec.kcat`` and
    rewrites just the tunable rows' stoichiometry via
    ``apply_kcat_constraints(update_rxns=...)`` -- no per-particle
    ``EcModel.copy()`` (prohibitively expensive; see the "Spike
    results" section of ``docs/internal/bayesian_tuning_plan.md``).
    """
    def score(kcat_vec: np.ndarray) -> float:
        model.ec.kcat[tunable_idx] = kcat_vec
        apply_kcat_constraints(model, update_rxns=ec_rxn_ids_tunable)

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
        )
        return rmse

    return score
