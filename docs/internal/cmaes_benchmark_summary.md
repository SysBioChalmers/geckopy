# CMA-ES kcat tuning: current state and how it was benchmarked

This session replaced ABC-SMC with CMA-ES as the package's kcat-tuning
method, then benchmarked the result directly rather than trusting the
configuration it inherited. This document covers both: what the
tuning pipeline looks like now, and what the benchmark found for each
of its hyperparameters. The full research log, including everything
tried and discarded along the way, is
`matlab_replication_results.md`/`bayesian_tuning_handover.md`; this is
the distilled version.

## Why CMA-ES, not ABC-SMC

DLKcat's own Bayesian approach to kcat tuning implements ABC-SMC
(Approximate Bayesian Computation, Sequential Monte Carlo) --
[SysBioChalmers/DLKcat/BayesianApproach](https://github.com/SysBioChalmers/DLKcat/tree/master/BayesianApproach)
(`abc_matlab_max.m`). geckopy's ABC-SMC port was tested head to head
against CMA-ES on the same objective, the same screened parameter set,
and the same experimental data (ecYeastGEM, 41 conditions), and
ABC-SMC lost on every criterion checked:

| | distance | spread across seeds |
|---|---|---|
| ABC-SMC | 0.9156 | +/- 0.0549 |
| **CMA-ES** | **0.7974** | **+/- 0.0048** |

CMA-ES fits better by 3.7 standard errors, with a spread more than ten
times tighter. The practical failure mode was worse than the fit gap
suggests: both methods leave most kcats on flat directions the data
can't pin down, but ABC-SMC's posterior sampling wanders across the
whole flat region rather than converging to one point. Across three
ABC-SMC seeds, lanosterol synthase -- the single highest-leverage kcat
in the model -- took values of 0.338, 20.3, and 0.0781 1/s: a 260-fold
spread, at distances differing by 0.009. A number that unstable cannot
be reported as a finding. ABC-SMC (`bayesian_kcat_tuning`, plus
`selection.py`, `transition.py`, `diagnostics.py`, and
`priors.build_kcat_prior`, all of which existed only to support it)
has been removed from the package entirely, along with its `pyabc`
dependency. See `cmaes_vs_abc_smc.md` for the full comparison.

## What the tuning pipeline looks like now

Three functions in `kcat_sensitivity_analysis.bayesian`, normally used
in this order:

- **`screen_kcat_leverage(model, ...)`** -- reports which kcats the
  data can actually speak to. Perturbs each kcat (isozyme copies
  sharing a prior and a source move together, per `tie_isozymes`) up
  and down by a fold change, and ranks by the resulting change in RMSE
  weighted by how much the source is trusted. No optimisation --
  useful for curation on its own, before any tuning run.
- **`select_tunable_mask(model, screen, target_impact_share=0.9)`** --
  turns a screen into the set of kcats a tuning run should search
  over: the fewest highest-ranked kcats whose combined leverage
  reaches `target_impact_share` of the total. A relative cutoff by
  design, so the same setting selects a comparable quality of
  parameter set on a model with a different leverage scale, rather
  than an absolute leverage threshold or a fixed count tuned to one
  model's structure.
- **`cmaes_kcat_tuning(model, ...)`** -- the tuning run itself. Screens
  and selects automatically if no mask is given. Configuration is the
  same `BayesianParams` used throughout: trust tiers
  (`sigma0_log_source`/`sigma0_log_default`), `max_growth_weight`,
  `prior_penalty_weight`, `tie_isozymes`, and `max_generations`/
  `rmse_threshold` as CMA-ES's own stopping conditions.
- **`tune_prior_penalty_weight(model, candidates=..., seeds=...)`** --
  a calibration helper, added this session: sweeps
  `prior_penalty_weight` at a fixed tunable set and reports fit cost
  against cross-seed reproducibility (movers, direction agreement,
  worst spread) per candidate, so the penalty doesn't have to be
  picked blind. It reports only -- there is no single best candidate
  for it to leave the model in, so it restores `model.ec.kcat` before
  returning. This is what should be used to reproduce or extend the
  lambda benchmark below rather than hand-rolling it again.

Two defaults changed this session, both now shipped:
`prior_penalty_weight` 0.0 -> **0.03**, `tie_isozymes` `False` ->
**`True`**.

## How the hyperparameters were benchmarked

`run_benchmark.py` checks the settled operating point directly: one
script, two data regimes (the full 41-condition dataset, and a
single measured glucose growth rate standing in for a data-poor
adapter), always scored on all 41 conditions regardless of what a
cell was fitted on. Free parameter count, `prior_penalty_weight`, and
`max_growth_weight` were each swept one at a time against that same
fixed yardstick. Two seeds per cell throughout -- enough to say
whether a result reproduces at all, not a rigorous spread estimate.

### Free parameter count: 112 confirmed

| free params | full-41 distance (seed 0 / 1) | kcats changed | median fold |
|---|---|---|---|
| 224 | 0.8295 / 0.9235 | 268 / 271 | 1.49 / 1.52x |
| **112** | 0.8384 / 0.8302 | 135 / 134 | 1.38 / 1.35x |
| 56 | 1.0105 / 1.0080 | 78 / 76 | 1.56 / 1.59x |
| 28 | 1.3844 / 1.3915 | 41 / 42 | 1.83 / 1.78x |
| 14 | 1.7455 / 1.7457 | 24 / 24 | 2.39 / 2.41x |

Fit degrades smoothly all the way down -- no cliff, no plateau -- and
112 sits at the top: 224 doesn't fit better (within noise of 112) but
reproduces far worse (seeds 11% apart against 1%) for double the
change set. On the single-measurement regime, parameter count barely
matters at all (full-41 distance sits in a flat 1.57-1.73 band from
112 down to 14) -- with one number to fit, the ceiling is the data,
not the search budget.

### `prior_penalty_weight`: 0.03 confirmed, 0.04 a candidate

| lambda | full-41 distance (seed 0 / 1) | fit cost vs 0.01 | direction agree | median spread | max spread |
|---|---|---|---|---|---|
| 0.01 | 0.8248 / 0.8142 | -- | 84.6% | 2.03x | 16.7x |
| **0.03** | 0.8384 / 0.8302 | +1.8% | 100% | 1.84x | 5.7x |
| 0.04 | 0.8514 / 0.9239 | +8.3% | 100% | 1.12x | 6.6x |
| 0.05 | 0.9465 / 0.9462 | +15.5% | 100% | 1.20x | 4.6x |
| 0.1 | 1.0575 / 1.1624 | +35.4% | 100% | 1.13x | 2.2x |

Direction agreement -- whether corrections beyond two-fold move the
same way across seeds -- saturates at 0.03 and never improves further.
*Magnitude* agreement keeps improving past that: median cross-seed
spread drops from 1.84x to 1.12x between 0.03 and 0.04, then goes flat
(1.12x-1.20x) all the way to 0.1, while fit cost keeps climbing
throughout (+8.3% -> +15.5% -> +35.4%). Past 0.04, the penalty is
mostly just spending fit for no further typical-case gain.

0.03 stays the shipped default: cheapest fit cost of any penalised
setting, already at 100% direction agreement. 0.04 is a real candidate
for a stronger default, but the case rests on a median computed over
15 both-moved kcats from two seeds, and 0.04's *worst*-case spread
(6.6x) is actually the highest of the four penalised settings tested
-- plausibly one outlier kcat at this sample size, not confirmed
either way. Resolve with a third seed before changing the default.

### `max_growth_weight`: 2 confirmed, with a caveat

The untuned baseline shifts with `max_growth_weight` (the same errors,
weighted differently), so compare normalised improvement, not raw
distance (untuned flux RMSE 8.796, growth RMSE 8.208):

| mgw | flux improvement | growth improvement | direction agree | median spread | max spread |
|---|---|---|---|---|---|
| 1 | 89.7% / 89.4% | 89.7% / 89.6% | 100% | 1.11x | 4.0x |
| **2** | 89.6% / 90.8% | 90.2% / 89.8% | 100% | 1.84x | 5.7x |
| 4 | 89.4% / 87.4% | 90.4% / 90.3% | 100% | 1.30x | 11.7x |

The intended effect -- growth data, spanning eight diverse carbon
sources against a flux dataset that is 30 of 33 conditions glucose,
should count for more -- is present but small (growth improvement
rises gently from 1 to 4) and mgw 4 pays for it with the widest
worst-case spread of the three. mgw 1 actually reproduces tighter than
2, but that is because it asks for no growth emphasis at all, which
the data's own asymmetry argues against. 2 stays the default: not the
winner on raw fit (a three-way tie within noise), but the point that
delivers the intended design choice at a bounded reproducibility cost.

## Recommended configuration

```toml
[bayesian]
sigma0_log_default = 0.3     # model-specific, no universal default
max_growth_weight = 2.0
prior_penalty_weight = 0.03  # shipped default; see note above on 0.04
tie_isozymes = true          # shipped default
max_generations = 212        # ~6800 CMA-ES evaluations at this model's screened size
```

`target_impact_share` (the `select_tunable_mask` cutoff that replaced
a hardcoded parameter count) is at its own default of 0.9, which
happens to select close to 112 kcats on ecYeastGEM but has not been
validated on a model with a different leverage distribution -- see
`matlab_replication_results.md`, Open items #6.

## Where the loose ends are tracked

- `matlab_replication_results.md`, Open items #6 and #7: the
  `target_impact_share` default's portability, and the `0.04` lambda
  candidate.
- `docs/bayesian_kcat_tuning.md`: the user-facing walkthrough for
  running this pipeline on a new model.
- `docs/cmaes_vs_abc_smc.md`: the ABC-SMC comparison on its own.
