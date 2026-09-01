# Bayesian (ABC-SMC) kcat tuning — port plan

Status: resumed on `feat/bayesian-kcat-tuning`. Open decisions below are
now resolved (see each item) after a fuller review of MATLAB's tuning
history and a closer read of pyABC's actual generation-loop mechanics.

## What it is

MATLAB source:
`GECKO/src/geckomat/kcat_sensitivity_analysis/Bayesian/`
- `bayesianSensitivityTuning.m` (856 LOC) — the main loop
- `abc_max.m` (121 LOC) — distance function (FBA + RMSE)
- `loadBayesianData.m` (19 LOC) — loads the 3 experimental TSVs
- `getrSample.m`, `updateprior.m`, `mse.m`, `addCarbonNum.m` — helpers
- ~1370 LOC total

It is an **ABC-SMC** (Approximate Bayesian Computation, Sequential
Monte Carlo) parameter-estimation routine that searches for kcat
values letting an ecModel reproduce experimental flux / growth data.

## Algorithm

Vanilla ABC-SMC core:
1. Lognormal prior per kcat, with per-source initial sigma.
2. Each generation: sample N proposals -> apply each kcat set to the
   ecModel -> solve FBA -> compute RMSE vs experimental data -> accept
   best ~10% (`targetAccept` percentile, gated by `minKeep`/`maxKeep`).
3. Update posterior mean/sigma in log-space from accepted samples.
4. Adaptive epsilon (quantile of the RMSE distribution).
5. Stop when RMSE < `rmseThreshold` or `maxGenerations` reached.

Distance function (`abc_max.m` -> `rmsecal`):
- RMSE between simulated FBA fluxes/growth and experimental
  measurements across multiple carbon-source conditions.
- Two datasets: max growth rates, and exchange fluxes.
- Zero-flux constraints on reactions assumed inactive.
- Carbon-normalised exchange fluxes (uses `excarbon` per reaction).
- Anaerobic switch + protein-content adjustment per condition.
- Failed/infeasible solves get RMSE = 99 (high penalty).

## GECKO-specific extensions (the load-bearing custom code)

These sit on top of the standard SMC and are what make it work on real
models:
- **Per-source sticky priors**: each kcat source (`brenda`, `dlkcat`,
  `custom`) has its own initial sigma, shrinkage threshold, variance
  cap, and force-to-prior cutoff. BRENDA kcats stay near prior;
  DLKcat ones drift freely.
- **Low-rank PCA proposal kernel** (`proposeLowRankMixture`): learns
  the dominant directions of posterior variability (SVD in log-space)
  and proposes along them.
- **Explore/exploit mixture**: alpha=0.65 PCA-guided exploit kernel +
  0.35 broad isotropic explore kernel (`cExpl` inflation).
- **Sparsity enforcement**: after each update, kcats whose change is
  `< sparsityThr * sigma0` snap back to prior.
- **Variance cap per source**: caps posterior/prior variance ratio to
  prevent random walk.
- **Diagnostics**: shrinkage trace, acceptance rate, epsilon, sparsity
  count, diversity, per-source metrics, low-rank dim — mostly consumed
  by tutorial plots.
- **Optional post-run pruning**: parameters whose change is below a
  threshold revert to prior (`enablePruning`).

Hyperparameters: ~15 in `params.bayesian`
(sigma0logDefault, kcatSources, sigma0logSource, shrinkThr*,
varianceCap*, forcePriorThr*, sparsityThreshold, scheduleGenerations,
scheduleSamples, targetAccept, minKeep, maxKeep, rmseThreshold,
maxGenerations). Fixed internals: alpha, tauResidual, sigmaFloorFrac,
adaptFracEarly, cExpl, rMax.

## Experimental data files (per model adapter, in data/)

- `bayesianFluxData.tsv` — carbon sources, growth rates, env
  conditions, exchange rates
- `bayesianMaxGrowth.tsv` — max growth on various carbon sources
- `bayesianZeroExch.tsv` — exchange reactions assumed zero-flux

## Library choice: pyABC (recommended)

Compared pyABC vs ABCpy vs ELFI:

| Library | Strength | Gap for us |
|---|---|---|
| **pyABC** | custom priors, custom Transition kernels, adaptive epsilon, dynamic-scheduled parallelism (dask/MPI/redis/mp), used in systems biology | no built-in per-source shrinkage or low-rank kernel — subclass `Transition` |
| ABCpy | widest algorithm menu | static-scheduled parallelism only |
| ELFI | BOLFI GP surrogate could cut FBA calls | research-grade; high-D kcat space may not suit BOLFI |

pyABC chosen because (a) its `Transition` abstraction is the right
hook for the GECKO-specific kernel, and (b) dynamic scheduling matters:
each FBA solve is ~10-100 ms and the schedule is up to 500 samples x
50-100 generations, so wall-time dominates. MATLAB already uses
`parfor`; pyABC's dask backend matches that.

Refs:
- pyABC docs: https://pyabc.readthedocs.io/en/latest/
- pyABC JOSS paper (2022): https://arxiv.org/pdf/2203.13043
- ABC-SMC practical guide (2025): https://arxiv.org/html/2511.21587

### What pyABC absorbs vs what we still write

pyABC gives: SMC loop, adaptive epsilon (`MedianEpsilon` ~ GECKO's
quantile rule), parallel sampling, SQLite-backed posterior + per-gen
bookkeeping. ~30-40% of LOC.

We still write (~60%): the FBA+RMSE distance (`abc_max.m` port),
per-source shrinkage / sparsity / force-to-prior / variance-cap
callback, the low-rank PCA `Transition` subclass, the 3-TSV data
loader, and the extra GECKO diagnostics.

Without pyABC we'd write that same 60% plus a 200-400 LOC SMC core
with worse parallelism — so pyABC is worth the dependency.

## geckopy module shape

```
src/geckopy/kcat_sensitivity_analysis/bayesian/
    __init__.py            # re-exports; guards pyabc absence (ImportError on call, not import)
    data.py                # load bayesian{FluxData,MaxGrowth,ZeroExch}.tsv
    simulate.py             # kcat vector -> apply_kcat_constraints -> FBA -> per-condition results
    distance.py              # carbon/condition-weighted RMSE (port of abc_max.m's rmsecal half)
    priors.py                # per-kcat lognormal priors by source_group; sparsity-prior variant
    selection.py              # both selection variants: truncation_select / quantile_epsilon_select
    posterior.py              # regularization variant A: shrink/force-to-prior/sparsity-snap blend
    importance_weights.py      # regularization variant B: prior/transition-density particle weights
    transition.py              # GeckoTransition(pyabc.transition.Transition): diagonal fit/rvs/pdf
    diagnostics.py              # per-generation/per-group diagnostics, variant-agnostic
    tuning.py                  # bayesian_kcat_tuning(model, ..., selection=, regularization=)
    pruning.py                  # optional post-run pruner (deferred, see Q5)
    cli.py                      # geckopy bayesian-tune subcommand (deferred, see Q6)
tests/
    test_bayesian_data.py
    test_bayesian_simulate_distance.py
    test_bayesian_selection.py
    test_bayesian_posterior.py
    test_bayesian_importance_weights.py
    test_bayesian_transition.py
    test_bayesian_tuning.py                 # tiny synthetic EcModel, all 4 variant combinations
    test_bayesian_tuning_smoke.py           # @pytest.mark.smoke, real-scale
```

## Open decisions (resolved)

- **Q1 Library**: pyABC, confirmed. But **not** as the top-level driver:
  `ABCSMC.run()`'s epsilon scheme (`QuantileEpsilon`/`MedianEpsilon`) is
  *prospective and streaming* — epsilon for generation *t* is fixed from
  generation *t-1*'s accepted distances, then the sampler draws proposals
  until a target *accepted* count is reached (`sample_until_n_accepted`,
  unbounded simulation budget). MATLAB's `minKeep` is the opposite: a
  *fixed* batch of `scheduleSamples[gen]` simulations, then a *retrospective*
  top-fraction cut on that same batch. These don't map onto each other
  without fighting the library. Resolution: use pyABC's `Transition`/
  `Distribution` as component libraries inside a small custom generation
  loop geckopy owns, not `ABCSMC.run()`.
- **Q2 Scope**: neither "MVP, no extensions" nor "full port of every
  extension" — informed by a detailed review of MATLAB's actual tuning
  history (not just its code), several "extensions" turned out to be either
  abandoned by MATLAB itself or empirically dead/cosmetic:
  - The low-rank PCA proposal kernel (`proposeLowRankMixture`) was tried and
    *abandoned* in favor of "the SIMPLIFIED VERSION with diagonal proposals"
    that's actually committed to HEAD. Not porting it; pyABC's adaptive
    `LocalTransition`/`MultivariateNormalTransition` is the presumptive
    diagonal-proposal default instead. Also the right shape at genome scale:
    `ec.kcat` has thousands of entries vs. a few hundred particles per
    generation, so a low-rank/full covariance would be rank-deficient.
  - `targetAccept` is confirmed dead (`minKeep`'s floor always binds first);
    `varianceCapDefault`/`varianceCapSource` are confirmed cosmetic (clamp
    only a *reported* diagnostic, never feed back into proposals or point
    estimates). Not porting either — dropping `target_accept`/
    `variance_cap_*` from `BayesianParams` outright rather than shipping
    fields that look load-bearing but silently aren't.
  - `paramTighteningGen`, `shrinkThrDefTight`, `forcePriorThrDefTight` are
    dead leftovers in MATLAB itself, not present in `BayesianParams`
    already — not adding them.
  - What *is* genuinely load-bearing and kept: `minKeep` (the single most
    consequential, empirically-tuned selection lever), `sparsityThreshold`
    (fixing its Python default from 0.3, inside the confirmed-bad
    0.25–0.28 range, to the empirically-working 0.5), per-source
    sticky-prior regularization (mechanism now built as two comparable
    variants — see below).
  - **New**: rather than assuming either "MATLAB's exact mechanism" or "the
    generic ABC-SMC textbook mechanism" is right without evidence, two axes
    are built as swappable variants and compared empirically on a toy
    problem + real-scale smoke test before picking a default:
    - *Selection*: MATLAB-faithful fixed-batch top-`min_keep`-fraction
      truncation vs. pyABC-native adaptive quantile-epsilon streaming
      acceptance.
    - *Regularization*: MATLAB-faithful shrink-weight/force-to-prior blend
      vs. proper SMC importance weighting (`prior_density / transition_density`,
      standard SMC-ABC theory) — a tight per-group prior then produces
      shrinkage-to-prior automatically, no hand-tuned thresholds needed.
      Sparsity enforcement rides along this axis too: MATLAB's post-hoc
      snap-to-prior vs. a genuine sparsity-inducing prior (spike-and-slab/
      horseshoe-style).
    See `~/.claude/plans/encapsulated-kindling-sphinx.md` for the full
    comparison design (this doc doesn't duplicate it).
- **Q3 Parallel backend**: multiprocessing, not dask — reuse the
  `multiprocessing.Pool` pattern already in `src/geckopy/utilities/ec_fva.py`
  (known-working precedent for parallel LP-heavy ecModel workloads in this
  codebase, no new dependency). Revisit only if profiling shows it's the
  bottleneck.
- **Q4 Return shape**: confirmed — mutate `model.ec.kcat` in place, return
  one `BayesianTuningResult` dataclass, matching the `TunedKcatsResult`/
  `SigmaFitterResult` precedent in `sensitivity_tuning.py`/`sigma_fitter.py`.
- **Q5 Pruning**: deferred, not same-PR — validate the core loop first
  (lower risk). Will eventually have two paths mirroring the regularization
  axis: MATLAB's one-at-a-time re-simulation loop, or a free read-off-the-
  posterior check when the sparsity prior is active.
- **Q6 Tutorial**: unchanged — follow-up PR, after the core module and its
  variant comparison are validated.

Also resolved: `ec.source` doesn't carry a generic `"okp"` tag —
OpenKineticsPredictor kcats are tagged with the raw predictor method name
(e.g. `"CataPro"`). `BayesianParams.kcat_sources` (a flat literal-string
list) can't express "any OKP method" as one trust tier. Resolution: replace
it with a `source_groups` concept — a mapping from trust-tier name to either
explicit literal source strings or a rule matching any configured OKP method
name (read from `OkpParams`) as one group; unmatched sources fall into an
implicit `"unlabelled"` tier (matching MATLAB's `noKcatSource`).

## Dependency note

pyABC pulls in scipy, pandas, cloudpickle (no dask needed — see Q3). Add as
an optional extra: `pip install geckopy[bayesian]` so the core install
stays lean (mirrors the proposed WILDkCAT / wildkcat optional extra).
