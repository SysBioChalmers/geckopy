# Why this package tunes kcats with CMA-ES, not ABC-SMC

DLKcat's own Bayesian approach to kcat tuning implements ABC-SMC
(Approximate Bayesian Computation, Sequential Monte Carlo) --
[SysBioChalmers/DLKcat/BayesianApproach](https://github.com/SysBioChalmers/DLKcat/tree/master/BayesianApproach).
geckopy's `kcat_sensitivity_analysis.bayesian` module ported that
approach faithfully and tested it side by side with a plain optimiser,
CMA-ES, on the same objective, the same screened parameter set, and
the same experimental data (ecYeastGEM, 41 conditions). CMA-ES won on
every criterion that was checked, so it is what this package ships.

## What ABC-SMC does, and why that's the wrong tool here

ABC-SMC is a sampler: each generation it draws a batch of candidate
kcat vectors from a proposal distribution, scores them, and keeps the
best fraction to seed the next generation's proposals. That machinery
-- a posterior approximation, a sample schedule, truncation selection
-- earns its keep when the goal is a *distribution* over plausible
parameter values, and the model is cheap enough to sample densely
across the whole parameter space.

Once a screen has reduced the problem to the kcats the data can
actually inform (a few dozen to a few hundred, out of thousands), that
premise stops holding. There's no longer a broad space to sample
across, and each evaluation is expensive (an FBA solve per condition).
At that point the problem is a direct search for the best-fitting
vector -- exactly what an optimiser is for, and exactly what ABC-SMC
is not.

## What the comparison showed

Same model, same 112-parameter screen, same objective, three seeds
each:

| | distance | spread across seeds |
|---|---|---|
| ABC-SMC | 0.9156 | +/- 0.0549 |
| **CMA-ES** | **0.7974** | **+/- 0.0048** |

CMA-ES fits better by 3.7 standard errors, with a spread more than ten
times tighter. The gap isn't just fit quality -- it's whether an
individual tuned kcat means anything. Both methods leave most kcats on
flat directions the data can't pin down, but ABC-SMC's posterior
sampling wanders across the whole flat region: across three ABC-SMC
seeds, lanosterol synthase (the single highest-leverage kcat in the
model) took values of 0.338, 20.3, and 0.0781 1/s -- a 260-fold spread,
at distances differing by 0.009. A number that unstable cannot be
reported as a finding.

CMA-ES doesn't remove flat directions -- nothing can, they're a
property of the data, not the search method -- but a direct optimiser
converges to one point per flat direction instead of sampling across
it, and `prior_penalty_weight` (a Tikhonov-style penalty on moving away
from the prior, weighted by how much each source is trusted) picks out
which point: see `docs/bayesian_kcat_tuning.md` for how that combines
with `screen_kcat_leverage`/`select_tunable_mask` into the method this
package actually recommends.

## Details

The full comparison, including the screening methodology both methods
were run against, is in `docs/internal/bayesian_tuning_handover.md`
("Three methods, one conclusion") and
`docs/internal/matlab_replication_results.md`.
