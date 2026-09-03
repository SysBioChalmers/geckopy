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

## The prior-penalty sweep

`prior_penalty_weight` adds `rmse + lambda * mean((log(k/k0)/sigma0_log)^2)`
to selection, so parsimony is optimised rather than applied afterwards.
Selection uses that objective; `rmse_trace` stays the plain RMSE and
`objective_trace` records what was minimised, so runs remain comparable
to each other and to MATLAB at any lambda.

Calibration first, since the penalty must be scaled against the RMSE:
the best particle has `mean(dev^2)` of 88, so the penalty equals the
RMSE at lambda 0.0103, the diffuse particle loses to the sparse blend
at lambda ~0.041, and above lambda ~0.087 the search never leaves the
prior at all. The sweep brackets that.

Four runs, 15 generations each, seed 0, otherwise identical:

| lambda | RMSE | mean dev | changed | knee | objective flat from |
|--------|------|----------|---------|------|---------------------|
| 0.003 | 1.0012 | 5.13 sigma | 4767 | 4237 @ 0.9752 | gen 14 |
| 0.01  | 1.0506 | 4.68 sigma | 4782 | **2434 @ 1.0096** | gen 11 |
| 0.03  | 1.3780 | 3.01 sigma | 4734 | 3779 @ 1.3458 | **gen 6** |
| 0.06  | 1.3515 | 2.99 sigma | 4737 | 4737 @ 1.3515 | **gen 6** |

RMSE here is the *returned* vector's, read off the frontier at `t=0`.
At lambda > 0 the returned vector is the best-by-objective particle,
whose plain RMSE is not `rmse_trace[-1]` -- that is the lowest RMSE
among accepted particles, which is a different particle. Reporting the
latter as the returned vector's fit understates it by up to 27%.

### What it bought, and what it did not

Movement falls monotonically with lambda, 5.13 -> 4.68 -> 3.01 -> 2.99
sigma, against 7.55 for the unpenalised run. Fit is flat to lambda 0.01
and then costs about 28%. That is a real, well-behaved trade-off curve,
and lambda 0.01 is better than the unpenalised run on *both* axes.

Two things it did not do.

**It does not produce a sparse solution.** Every lambda still moves
98%+ of all kcats. The L2 penalty compresses how *far* parameters move,
not how *many*. Genuine sparsity needs an L0-shaped mechanism -- block
or coordinate proposals that perturb a small random subset per draw --
not a smooth penalty.

**Above lambda ~0.03 the search stalls rather than economises.** Both
0.03 and 0.06 freeze their objective at generation 6 and never improve
it again, with proposal acceptance falling to 0.077 and 0.045. Their
near-identical movement (3.01 and 2.99 sigma) is not two searches
choosing similar parsimony; it is two searches that stopped. Low
movement from a frozen search must not be read as a parsimony win --
check `objective_trace` for a plateau before believing any such number.

### Recommendation

**lambda = 0.01.** It has by far the best frontier -- 2434 changed at
RMSE 1.0096, against 3716 for the unpenalised run -- its objective was
still improving to generation 11 rather than stalling, and it costs
nothing in fit. The value sits just below the point where the penalty
equals the RMSE, which is where a regulariser usually belongs.

Two caveats on the comparison. The unpenalised movement figure (7.55
sigma) comes from a 31-generation run while these are 15, so it is an
upper bound rather than a like-for-like control; a 15-generation
lambda 0 run would settle it. And each cell is one seed, so differences
under a few percent in RMSE are inside the noise measured elsewhere in
this document.

## Identifiability: most of the changes are provably meaningless

Two screens already existed from earlier work and answer the parsimony
question far better than the penalty sweep did.

**FVA reachability.** Per-condition ecFVA across all 41 conditions:
3955 of 4834 kcats belong to reactions that can carry flux in at least
one condition. **879 can never carry flux in any of them**, so no
datum in this dataset constrains their kcat. The relaxed-protein-pool
variant is a subset (3950) and adds nothing.

Every run changes them anyway: 873 of 879 in the corrected run, 876 in
the B1 run, 868 at lambda 0.01. Those edits cannot be justified by the
data under any weighting, and no amount of RMSE improvement excuses
them.

**One-at-a-time sensitivity.** Perturbing each reachable kcat by x2
either way at the prior moves the distance by more than 1e-3 for 1939
of them, more than 1e-2 for only 81.

### Reverting by importance instead of by magnitude

`parsimony.py`'s frontier reverts the *least-moved* kcats, which is a
poor proxy: in a diffuse solution a kcat can move 8 sigma and matter
not at all. Reverting by screened sensitivity instead is dramatically
better:

corrected run (full vector: 4798 changed, RMSE 0.9093, 7.55 sigma)

| keep abs(dRMSE) > | changed | RMSE | mean dev |
|-------------------|---------|------|----------|
| 1e-4  | 2736 | **0.9069** | 4.32 sigma |
| 1e-3  | 1925 | 0.9261 | 3.01 sigma |
| 3e-3  | **297** | 1.0759 | 0.46 sigma |
| 1e-2  | 80 | 1.2726 | 0.13 sigma |

lambda 0.01 run (full vector: 4782 changed, RMSE 1.0506, 3.86 sigma)

| keep abs(dRMSE) > | changed | RMSE | mean dev |
|-------------------|---------|------|----------|
| 1e-4  | 2730 | 1.0500 | 2.69 sigma |
| 1e-3  | **1920** | **1.0287** | 1.87 sigma |
| 3e-3  | 301 | 1.2340 | 0.30 sigma |

Two results stand out. Reverting everything below 1e-4 *improves* both
runs while dropping ~2000 changes -- those changes were not merely
unjustified, they were mildly harmful. And 297 changed kcats retain an
RMSE of 1.076 against the prior's 8.60, where magnitude-ranked
reverting needed 527 changes to reach only 2.97.

So the answer to "far fewer kcats should change" is not a stronger
penalty. It is that the search is free to edit thousands of parameters
the data cannot see, and nothing -- not the penalty, not the shrinkage,
not the sparsity snap -- ever consults whether a parameter is
identifiable at all.

## Restricting the search to identifiable kcats

`bayesian_kcat_tuning(tunable_mask=...)` holds excluded kcats at their
prior. The FVA mask -- reactions able to carry flux in at least one of
the 41 conditions -- has **no hyperparameters**: it is derived from the
model and the dataset, not chosen.

### It beats pruning afterwards, at every matched sparsity

At equal budget (15 generations), restricting the search against
reverting the unrestricted run by the same screen:

| changed | prune afterwards | restrict the search |
|---------|------------------|---------------------|
| ~1920 | 1.2193 | **1.0978** |
| ~2730 | 1.2015 | **1.0651** |
| ~3915 | 1.0778 | **0.9427** |

Restriction wins by 10-13% throughout. An earlier comparison appeared
to show the opposite only because it set a 31-generation pruned run
against 15-generation restricted ones; at equal budget the ordering is
consistent. Deciding before the search which parameters the data can
speak to beats tuning everything and discarding later, because the
budget is spent where it counts instead of diffusing across dimensions
that will be thrown away.

### Result

31 generations, FVA mask, no other change:

| | RMSE | changed | mean dev |
|---|------|---------|----------|
| MATLAB | 0.8715 | ~4804 | 8.30 sigma |
| port, unrestricted | 0.9093 | 4798 | 7.55 sigma |
| **port, FVA-restricted** | **0.8688** | **3921** | 7.55 sigma |

Better than MATLAB on fit and 883 fewer changed kcats, from a mask with
nothing to tune. Pruning the restricted result by sensitivity then
gives the frontier: 1920 changed at 0.9255, or 301 changed at 1.0124,
against the prior's 8.6002.

### Restriction and the penalty do different things

The FVA-restricted run has mean dev 7.55 sigma -- identical to the
unrestricted run. Restriction cut the *count* of changes without
touching their magnitude. The prior penalty did the reverse: it
compressed magnitude (7.55 -> 4.68 sigma) while leaving the count at
98% of all kcats. They are complementary rather than alternatives,
which is not how the penalty sweep framed them.

The one datapoint on combining them is discouraging -- `sens3` plus
lambda 0.01 stalled at generation 8 -- but that stacked two tight
constraints at half budget. A loose mask with a mild penalty at full
budget is a different proposition and is worth testing before either
is adopted as a default.

### Recommendation

Restrict to FVA-reachable kcats by default. It has no hyperparameter,
it improves both axes simultaneously, and it removes a class of change
-- editing kcats no condition can constrain -- that is indefensible
regardless of what the objective rewards. Report the sensitivity
frontier alongside the result so the fit/parsimony trade-off is chosen
explicitly rather than left implicit.

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
