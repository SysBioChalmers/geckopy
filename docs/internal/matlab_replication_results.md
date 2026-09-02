# Results: reproducing MATLAB's Bayesian kcat tuning in Python

Three 31-generation runs of the ported `bayesian_kcat_tuning` against
GECKO MATLAB's `bayesianSensitivityTuning`, plus the fidelity work that
made them comparable. Read `matlab_replication_handoff.md` for setup.

Headline: **the port's distance function now matches MATLAB's to
0.002%**, and on that corrected metric the port reaches 0.9093 against
MATLAB's 0.8715. The earlier apparent 38% shortfall was almost entirely
a scoring artifact, not a search one.

Second headline, and the more important one: **judged on how many kcats
move rather than on RMSE, every "best particle" result is unusable** --
MATLAB's included. Only the blend respects the principle that a tuner
should leave a kcat at its prior unless the data moves it.

## The reference target

Two MATLAB numbers were in circulation. They are not two runs of the
same thing:

* **8.6001 -> 0.8715**, `tutorials/full_ecModel/output/rmse_trace.tsv`
  with `run_metadata.txt`, committed in `8a6b7032` (2026-06-11). This
  is the reference: same branch, same model file.
* **9.5544 -> 0.9289**, from a report generated 2026-05-20 and quoted
  in the hand-off. It predates the model update in `cce743ee`
  (2026-06-09) and therefore describes a *different model*. It is
  stale; do not compare against it.

`kcatTrace[:,1]` in `bayesian_output.mat` is bit-identical to the
port's `kcat0` (all 4834 entries, rtol 1e-6), so the two
implementations demonstrably start from the same vector.

## The distance function: two missing pieces

`abc_max.m` had two behaviours the port lacked. Each one *alone* makes
agreement worse, which is why they went unnoticed and why an early
attempt to add the weighting on its own was abandoned as a dead end:

| scoring MATLAB's own prior | result | vs MATLAB's 8.60005 |
|---|---|---|
| neither | 9.6229 | +11.9% |
| `maxGrowthWeight` only | 10.0944 | +17.4% |
| anaerobic switch only | 8.5022 | -1.1% |
| **both** | **8.6002** | **+0.002%** |

* **Anaerobic switch.** `abc_max.m:81-84` calls
  `modelAdapter.makeModelAnaerobic` for every condition whose measured
  oxygen exchange is 0 -- 9 of the 33 flux conditions here, previously
  simulated with free oxygen. Ported as
  `tutorials/full_ecModel/code/anaerobic.py`.
* **`max_growth_weight`.** `abc_max.m:48` applies
  `weights = [maxGrowthWeight, 1]` to `[rmse_flux, rmse_maxGrate]`, so
  despite its name it scales the *flux* term.
  `YeastGEMAdapter.m` sets 2. Default stays 1 (a plain mean).

Verified on six independent MATLAB vectors: its prior within 0.002%
and five of its six best accepted particles within 0.02% (one outlier
at 0.75%, consistent with an alternate LP optimum).

`changeProteinBiomass` needs no port: `YeastGEMAdapter.m` defines it as
a no-op with the real call commented out.

## The three runs

All 31 generations, `seed=0`, 63 workers on one 64-core node,
`Threads=1`.

| gen | corrected | B1 adaptive | MATLAB |
|-----|-----------|-------------|--------|
| prior | 8.6002 | 8.6002 | 8.6001 |
| 4   | 2.8067 | 1.8815 | 2.9219 |
| 8   | 1.2514 | 1.1864 | 1.3790 |
| 11  | 1.2075 | 1.0204 | 0.9933 |
| 14  | 1.0778 | 1.0024 | 0.9933 |
| 22  | 0.9620 | 0.8633 | 0.8715 |
| 31  | **0.9093** | **0.8281** | **0.8715** |
| wall | 69.4 min | 71.6 min | -- |

The samples/accepted sequence reproduces MATLAB's generation for
generation in every run -- 1000/300, 800/330, ... 400/171, matching its
recorded `n_accepted_samples 171`. Schedule, pooling and truncation
selection are faithful; the remaining differences are not in the
sampler's bookkeeping.

## Parsimony: the criterion that reverses the ranking

A low RMSE bought by rewriting every kcat is not a good answer. Changes
should be few, and concentrated in the parameters we trust least.
Movement is measured as `|log(k/k0)| / sigma0_log` -- MATLAB's own
`devFromPrior` -- which weights by confidence automatically, since
`sigma0_log` is 0.10 for `custom`, 0.20 `brenda`, 0.25 `okp`, 0.30
unlabelled.

| vector | RMSE | moved >2% | mean dev | custom untouched |
|---|---|---|---|---|
| MATLAB best particle | 0.8715 | 99.4% | 8.30 sigma | 1.0% |
| port best (corrected) | 0.9093 | 99.3% | 7.55 sigma | 1.4% |
| port best (B1) | 0.8281 | 99.3% | 8.24 sigma | 0.5% |
| MATLAB blend | not scored | 11.6% | 0.64 sigma | 100.0% |
| port blend (corrected) | 4.2665 | **8.7%** | **0.48 sigma** | 100.0% |
| port blend (B1) | 3.4627 | 41.8% | 5.25 sigma | 78.7% |

**Every best particle moves the average kcat about eight prior standard
deviations, and moves them indiscriminately.** MATLAB's moves `custom`
-- the most trusted tier -- hardest of all (8.95 sigma). These vectors
have no claim to being posteriors; their RMSE only shows that 4834 free
parameters can fit 41 conditions. The port and MATLAB behave
identically here, so this is a property of the method, not of the port.

The blend is the vector that honours the principle: `custom` untouched
entirely, `brenda` barely moved, what movement there is concentrated in
`okp` and unlabelled. The port's blend is *more* parsimonious than
MATLAB's (0.48 vs 0.64 sigma) while being built the same way.

### B1 fails on this criterion

`adapt_proposal_width` steers the proposal bandwidth by the fraction of
new proposals that survive selection, instead of MATLAB's
`0.5*std_obs + 0.5*sigma0log` -- which cannot fall below half the prior
width and *grows* as the accepted set spreads (MATLAB's own trace:
width 0.213 -> 1.037 while proposal acceptance collapses 0.10 ->
0.005). The controller works as designed: it widened to its clamp
through generation 5 while proposals were landing, then annealed to
0.219 as acceptance fell, holding acceptance near its 0.25 target for
the whole second half where both baselines had stalled.

It also produces the best RMSE of any run here, 0.8281, beating MATLAB.

**And it is a regression.** Its blend moves 41.8% of kcats at 5.25
sigma against the corrected run's 8.7% at 0.48 sigma -- an eleven-fold
increase in movement -- and it disturbs `custom` (78.7% untouched,
against 100% for both baselines) and `brenda` (73.8% against 98.5%).
It buys RMSE by scattering parameters, including exactly the ones that
should be left alone.

A plausible mechanism, not yet verified: the blend's sparsity comes
from shrinkage toward the prior, gated on `devFromPrior` against
per-source thresholds. Wider early proposals push the accepted cloud
far from the prior, so the shrink and force-to-prior gates stop firing
and the blend simply follows the particles. If that is right, any
change that widens exploration will damage parsimony unless the
shrinkage is re-tuned alongside it.

## What this means for future variants

Report both axes for every new strategy or hyperparameter, never RMSE
alone. `parsimony.py` in the run scratch computes the table above from
saved `.npy` vectors, MATLAB's export included, without re-solving.
A variant that improves RMSE while increasing sigma-movement, or while
touching `custom`/`brenda`, has not improved anything.

## Open items

1. **The blend is never scored or optimised.** MATLAB computes it every
   generation, reports its statistics, then discards it and returns the
   best particle; the port faithfully does the same. Nothing in either
   search optimises the quantity that actually matters -- a sparse,
   confidence-weighted change set that also fits. The blend at 4.27
   against the particle at 0.91 is the size of that gap, and closing it
   is the real problem.
2. **Seed spread is unmeasured.** One seed per variant cannot separate
   a better sampler from a luckier one; the corrected and MATLAB traces
   cross repeatedly, each taking one large step (MATLAB at generation
   11, the port at 20) between long flat stretches. Three or four seeds
   at reduced `max_generations` would bound it cheaply.
3. **Re-tune shrinkage for wider proposals**, if B1's exploration is
   worth keeping -- test the mechanism above before assuming it.
4. **OKP** is the source whose blend movement differs most between
   implementations (70.8% untouched here against MATLAB's 60.9%),
   pointing at source classification rather than the sampler.
5. **Performance.** `model.optimize()` builds every primal, reduced
   cost and shadow price when under 40 numbers are read, costing about
   2.0 s per particle; `apply_kcat_constraints` is another 27%,
   rebuilding caches and issuing ~9700 `add_metabolites` calls per
   particle. Neither changes any result, so both are verifiable by
   reproducing a trace exactly.
