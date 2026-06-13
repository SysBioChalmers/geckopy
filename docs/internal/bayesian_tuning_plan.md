# Bayesian (ABC-SMC) kcat tuning — port plan (PAUSED)

Status: investigation done, implementation not started. Paused in
favour of the OpenKineticsPredictor work. Resume from here.

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

## Proposed geckopy module shape

```
src/geckopy/kcat_sensitivity_analysis/bayesian/
    __init__.py
    data.py            # load bayesian{FluxData,MaxGrowth,ZeroExch}.tsv
    distance.py        # FBA simulation + RMSE (port of abc_max.m + rmsecal)
    transition.py      # GeckoTransition(pyabc.Transition): per-source + low-rank PCA
    posterior.py       # per-generation post-update: shrinkage / sparsity / variance cap
    tuning.py          # bayesian_kcat_tuning(model, ...) orchestrator
    pruning.py         # optional post-run pruner
    cli.py             # geckopy bayesian-tune subcommand
tests/
    test_bayesian_data.py
    test_bayesian_distance.py
    test_bayesian_transition.py
    test_bayesian_posterior.py
    test_bayesian_tuning_integration.py   # tiny ecTestGEM-style fixture
```

## Open decisions (resolve on resume)

- **Q1 Library**: pyABC (recommended) unless reconsidered.
- **Q2 Scope**: MVP (vanilla SMC, no extensions — runs but won't match
  MATLAB) vs full port (all extensions). Recommended: full port,
  staged on a feature branch.
- **Q3 Parallel backend**: dask (recommended) / multiprocessing / MPI.
- **Q4 Return shape**: MATLAB returns 7 outputs
  (ecModel, rmseTrace, kcatTrace, sigmaLogTrace, diagnostics,
  posteriorSamples, prunedModel). geckopy convention: mutate
  `model.ec.kcat` in place + return one `BayesianResult` dataclass.
- **Q5 Pruning**: port in same PR, gated behind `prune=True`.
- **Q6 Tutorial**: port the stage-4 tutorial call in a follow-up PR.

## Dependency note

pyABC pulls in scipy, pandas, cloudpickle, and optionally dask. Add as
an optional extra: `pip install geckopy[bayesian]` so the core install
stays lean (mirrors the proposed WILDkCAT / wildkcat optional extra).
