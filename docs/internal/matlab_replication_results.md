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
deviations.** These vectors have no claim to being posteriors; their
RMSE only shows that 4834 free parameters can fit 41 conditions. The
port and MATLAB behave identically here, so this is a property of the
method, not of the port.

An earlier version of this section read the sigma figures as saying
MATLAB moves `custom` -- the most trusted tier -- hardest of all. That
was an artifact of the normalisation and is **wrong**: `sigma0_log` is
0.10 for `custom` and 0.30 for unlabelled, so equal movement in sigma
units is a three-fold *smaller* change for `custom`. In fold terms the
tiering works, and works in the right direction (see "Movement by
source" below). Report fold changes, not sigma multiples: only the
*ranking* of the sigmas across sources carries meaning, their absolute
values are a convention.

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

Combining them at full budget confirms this. `reach` plus lambda 0.01
over 31 generations gives RMSE 0.9813 with 3910 changed at 4.31 sigma,
against the mask alone at 0.8688 with 3921 changed at 7.55 sigma: the
same count, nearly half the magnitude, for 13% of the fit. It is also
the only run leaving whole groups untouched -- okp 31.8% unchanged,
unlabelled 24.8%, brenda 14.8%, where every other run is near zero.

The penalty still does not earn its place, because pruning the
lambda 0 run dominates it on the fit/count frontier at every level:
1920 changed at 0.9255 against the combination's 2158 at 1.1618, and
2728 at 0.9297 against 2551 at 0.9491. If the criterion is how many
kcats move, the penalty costs a hyperparameter and buys nothing.

### Recommendation

**Restrict to FVA-reachable kcats, prune the result by the sensitivity
screen, and leave `prior_penalty_weight` at 0.** No hyperparameter is
added: the mask is derived from the model and the conditions, and the
pruning threshold is an operating point on a reported curve rather than
something tuned per model.

| operating point | RMSE | changed |
|-----------------|------|---------|
| MATLAB | 0.8715 | ~4804 |
| best fit (mask, no pruning) | **0.8688** | 3921 |
| balanced (prune at 1e-3) | 0.9255 | 1920 |
| sparse (prune at 3e-3) | 1.0124 | 301 |

Every point beats MATLAB on parsimony; the first beats it on fit as
well. Report the whole curve, not one row: the trade-off is a choice
about the model, not a property of the algorithm.

Restriction also removes a class of change -- editing kcats no
condition can constrain -- that is indefensible regardless of what the
objective rewards, which is a stronger argument for it than the RMSE.

## Movement by source, in fold changes

Median fold change per source, most-trusted first, over *every* kcat of
the source (see "The FVA mask does not break the trust ranking after
all" below for why this reading misleads once a mask is in play, and
read the moved-only table there alongside it):

| vector | custom | brenda | okp | unlabelled | ordered? |
|--------|--------|--------|-----|------------|----------|
| MATLAB best | 2.13 | 3.98 | 6.37 | 6.36 | ~ |
| port unrestricted@31 | **1.81** | 3.69 | 5.25 | 7.96 | **yes** |
| port reach@31 | 2.21 | 3.62 | 2.87 | 5.84 | no |
| port reach+lambda@31 | 1.47 | 2.11 | 1.74 | 2.32 | no |

The unrestricted port respects the trust ranking monotonically and
moves `custom` less than MATLAB's own result does. But `custom` still
moves about two-fold at the median, with only ~9% left within 1.1x of
prior, so "trusted kcats stay near their original values" is not yet
true of any run here.

The FVA mask appears to break the ranking here -- `okp` moving less
than `brenda` -- but on the moved-only reading it does not; what it
breaks is the ranking of how many kcats each source moves. It is
source-blind, so removing 879 unreachable kcats reweights which tiers
dominate without reference to how far a source is trusted.
`identifiability_mask` is the response -- it requires a kcat's measured effect to clear a bar
proportional to its own `sigma0_log`, so a `custom` kcat needs three
times the effect of an unlabelled one before it is eligible to move at
all. One global threshold; the per-source differentiation is derived
from the sigma ranking rather than tuned.

## Few and consequential, not merely few

The goal is a small number of *impactful* changes, not a small number
of changes. `impact_share` reports the fraction of total screened
`|dRMSE|` carried by whatever a result changed:

| vector | changed | median fold | >2x | impact share |
|--------|---------|-------------|-----|--------------|
| MATLAB best | 4804 | 4.32 | 76% | 99.8% |
| MATLAB blend | 560 | 2.87 | 87% | **27.7%** |
| port reach@31 | 3921 | 5.06 | 78% | 99.4% |
| port reach+lambda@31 | 3910 | **2.58** | 61% | 99.6% |
| **reach@31 pruned at 3e-3** | **301** | **4.48** | 73% | **69.1%** |
| reach@31 pruned at 1e-3 | 1920 | 4.98 | 79% | 94.5% |

Two conclusions that reverse earlier advice in this document.

**The prior penalty is contraindicated.** It produces precisely the
unwanted shape: 3910 kcats changed at a median of 2.58 fold, only 61%
of them beyond two-fold. Many values nudged. It was already not
recommended on parsimony grounds; on impact it is worse than neutral
and should be dropped from the design rather than kept as a knob.

**The blend is worse than it looks.** 560 changes at 2.87 fold, but
carrying only 27.7% of the available impact: it is changing kcats the
data barely responds to while leaving load-bearing ones at their
priors. It has appeared defensible throughout because it scores well on
*count*. It does not survive an impact-weighted reading, which
undercuts the "blend is the principled posterior" framing in
`bayesian_tuning_plan.md`.

Sensitivity pruning gives the wanted shape directly: 301 kcats changed
at a median 4.5 fold, carrying 69% of available impact, RMSE 1.0124
against the prior's 8.6002 -- few, large, and demonstrably
load-bearing, because the screen selected them by impact.

## The trust-weighted mask

`identifiability_mask` at threshold 3e-4 leaves **1448** of 4834 kcats
eligible, against `reach`'s 3955. Because eligibility is
`|dRMSE| * sigma0_log > threshold`, a `custom` kcat must show three
times the measured effect of an unlabelled one to qualify. 31
generations, seed 0, 53.5 min:

| vector | RMSE | changed | impact | median fold (moved) |
|--------|------|---------|--------|---------------------|
| MATLAB best | 0.8715 | ~4804 | 99.8% | 4.32 |
| reach@31 | **0.8688** | 3921 | 99.4% | 5.06 |
| reach@31 pruned at 1e-3 | 0.9255 | 1920 | 94.5% | 4.98 |
| reach@31 pruned at 3e-3 | 1.0124 | **301** | 69.1% | 4.48 |
| **trust@31** | 0.8844 | **1433** | 88.5% | 3.69 |

It **beats the pruned operating point on fit and count at once**: 1433
changes at 0.8844 against 1920 at 0.9255. Its impact share is 6 points
lower for 25% fewer changes, so per change it concentrates more of the
screened effect, not less. Against the
unpruned mask it gives up 1.8% of fit for 2.7x fewer changes. The
search did not stall -- the trace was still stepping down at
generation 29 of 31 (0.9100 -> 0.8844), and with `lambda` 0 the
objective and RMSE traces coincide.

Movement by source, among the kcats each result actually moved:

| vector | custom | brenda | okp | unlabelled | ordered? |
|--------|--------|--------|-----|------------|----------|
| port unrestricted@31 | 1.83 | 3.71 | 5.32 | 7.96 | yes |
| port reach@31 | 2.52 | 4.98 | 7.67 | 11.45 | yes |
| **trust@31** | **2.70** | 3.65 | 4.51 | 7.57 | yes |

and how much of each source it leaves alone:

| vector | custom | brenda | okp | unlabelled |
|--------|--------|--------|-----|------------|
| port unrestricted@31 | 9.2% | 3.5% | 3.3% | 2.0% |
| port reach@31 | 9.7% | 17.4% | 32.1% | 26.8% |
| **trust@31** | **84.1%** | 71.0% | 72.0% | 62.2% |

This is the first run in which trusted kcats are mostly untouched: 37
of 207 `custom` kcats move, against 196 under the FVA mask alone. The
tier gradient now appears where it belongs, in *how many* kcats a
source is allowed to move -- 17.9% of `custom` against 38.6% of
unlabelled -- rather than only in how far they travel.

### The FVA mask does not break the trust ranking after all

The claim above that `reach` inverts the ranking (`okp` moving less
than `brenda`) came from a median taken over every kcat of a source.
That statistic saturates at 1.00 as soon as most of a source sits at
its prior, and it mixes count with magnitude. Among moved kcats
`reach@31` is monotone -- 2.52, 4.98, 7.67, 11.45 -- so what it breaks
is the ranking of *how many* kcats move per source (94.7% of `custom`
against 68.9% of `okp`), not how far they move. `source_movement`
reports both medians for this reason; read `median_fold_moved`
alongside `n_moved`, never `median_fold` alone.

Every mask run here passes on that reading -- unrestricted, all three
sensitivity screens, and both trust thresholds -- so no mask inverts
the trust ranking in magnitude; they differ only in the counts.

The trust mask under the same reading has an ordering that is barely a
gradient at the top (2.49 for `custom` against 3.45 for `brenda`),
which is the honest version of the bare "trust order respected:
YES" -- that flag was computed from four medians all equal to 1.00 and
meant nothing.

### What is not yet established

The 1.8% fit difference from `reach` is noise -- see "Seed spread"
below, which measured 5% spread between seeds of one configuration.
The mask neither wins nor loses on fit; it wins on count, which is the
reproducible axis.

The threshold 3e-4 is also a hyperparameter, which the FVA mask is not.
It buys the per-source differentiation, and it sets mask size the way
the pruning threshold sets change count -- an operating point on a
curve rather than a fitted value -- but the curve has one point on it
so far.

## Seed spread: what one run can and cannot show

Three seeds per mask, 15 generations, everything else fixed.

| mask | seed 0 | seed 1 | seed 2 | mean | sd | changed |
|------|--------|--------|--------|------|----|---------|
| reach | 0.9427 | 1.0324 | 1.0259 | 1.0003 | 0.0500 | 3910 +/- 5.5 |
| trust (3e-4) | 1.0353 | 0.9968 | 0.9317 | 0.9879 | 0.0524 | 1437 +/- 2.5 |

**Fit is not reproducible to better than about 5%; change count is
reproducible to 0.2%.** The two masks differ by 0.0125 +/- 0.0418 in
mean RMSE -- indistinguishable -- while differing by a factor of 2.7 in
how many kcats they move, and each mask's change count is stable across
seeds to a handful of parameters.

For two *single* runs the difference of RMSE carries sd ~0.071, so a
one-seed gap below roughly 0.14 says nothing. That threshold is larger
than most differences this document has ranked:

| comparison | gap | verdict |
|------------|-----|---------|
| trust vs reach @31 | 0.0156 (1.8%) | noise |
| port vs MATLAB @31 | 0.0378 (4.3%) | noise |
| B1 adaptive vs corrected @31 | 0.0812 (8.9%) | noise |
| prior-penalty sweep, between lambdas | <= 0.13 | noise |
| trust vs reach @15 | 0.0926 (9.8%) | noise |

The parsimony findings are unaffected, because they rest on quantities
that *are* reproducible: change counts, impact share, the per-source
movement pattern, and -- most of all -- mechanism. "This mask forbids
editing kcats no condition can constrain" is an argument about what the
search may touch, not a measured 2% of RMSE.

What this costs: the headline "the port reaches 0.9093 against MATLAB's
0.8715" is not evidence of a shortfall, and "B1 gives the best raw RMSE
of any variant" is not evidence that it searches better. Both need
seeds before they mean anything. What it buys: the trust mask no longer
has to defend a fit deficit, because at equal budget it does not have
one -- 1437 changes for the same distance as 3910.

Measured at 15 generations. Whether spread narrows by generation 31 is
unmeasured; until it is, a fit claim from one seed at any budget is an
anecdote.

## How concentrated the achievable fit is

The one-at-a-time screen answers this without any tuning run. Over the
3950 kcats it reaches, total `|dRMSE|` is 11.616, distributed with a
max/median ratio of 1972:

| impact captured | kcats needed | share of the model |
|-----------------|--------------|--------------------|
| 39.3% | top 10 | 0.2% |
| 50% | top 24 | 0.5% |
| 63.6% | top 110 | 2.3% |
| 74.0% | top 500 | 10.3% |
| 90% | top 1442 | 29.8% |
| 99% | top 2511 | 51.9% |

**Twenty-four kcats carry half the achievable improvement.** The
remaining half is spread across thousands of parameters each worth
about 1e-3 of RMSE, mutually substitutable as far as the distance is
concerned.

That is the mechanism behind the mask curve: past the first hundred
parameters the objective is nearly flat, so a search with 3910 free
kcats and one with 110 arrive at the same distance by different routes,
and the extra 3800 changes are movement the data never asked for. It
also explains why the prior penalty compressed magnitudes without
reducing counts -- there is no pressure toward any particular sparse
solution when so many are equivalent.

### The top of the impact distribution is enriched in `custom`

| source | share of top 110 by impact | share of the model |
|--------|---------------------------|--------------------|
| custom | **25.5%** | 4.3% |
| brenda | 57.3% | 66.2% |
| okp | 14.5% | 24.3% |
| unlabelled | 2.7% | 5.3% |

`custom` is six-fold over-represented among the highest-leverage kcats.
The likely reason is selection: a curated kcat exists because someone
found the model sensitive to it, so enrichment here is not evidence
that curated values are wrong.

It does sharpen what the trust mask is doing. Requiring `custom` to
clear a bar three times higher excludes, disproportionately, the
parameters with the most leverage -- so its parsimony is bought partly
by declining to touch high-impact kcats, which is a different
proposition from declining to touch kcats no condition can constrain.
Both are defensible; they are not the same argument, and the
`sigma0_log` ranking is what decides between them.

## Held-out conditions: what tuning does to data it did not fit

Conditions are held out in groups -- a whole study of flux conditions,
a whole carbon source of max-growth conditions -- because 30 of the 33
flux conditions are glucose and ten of those sit within 0.32-0.39 1/h
of each other, so a random row split leaves a near-duplicate of every
held-out condition in the training set. Three splits, each vector
scored on both halves:

| split | held out | prior | reach@31 | trust@31 |
|-------|----------|-------|----------|----------|
| A | AEM-1998 series, ethanol, acetate | **1.9292** | 3.9187 | 3.8738 |
| B | DLKcat study, galactose, maltose | 13.6301 | 1.6652 | 1.6863 |
| C | two small studies, mannose, sucrose | 5.7487 | 1.2205 | 1.2160 |

Tuning compresses every subset toward roughly 1: subsets the prior fits
badly improve by four- to eight-fold, and the one subset it fits well
degrades by two. The aggregate 8.60 -> 0.87 is therefore a mix of
genuine improvement and redistribution, and the balance depends on
which conditions are asked about. A tuned model predicts the AEM-1998
glucose series and ethanol/acetate growth *worse* than the untuned one.

`reach@31` and `trust@31` were tuned on all 41 conditions, so neither
column above is out-of-sample for them; the table shows where the
objective spends its budget, not generalisation.

### The one real held-out measurement so far

One run has been tuned on a split's training half and scored on the
half it never saw: shipped sigmas, FVA mask, 15 generations, split A.

| | train (23+6) | held out (10+2) |
|---|---|---|
| prior | 10.9548 | 1.9292 |
| tuned | 0.9734 | **4.1260** |

It fits the conditions it saw eleven times better and the conditions it
did not see twice worse than doing nothing. Split A is the least
favourable case -- it holds out the subset the prior already fits best,
so there is more to lose than to gain -- and runs on splits B and C are
queued to settle whether the result is about generalisation or about
those particular conditions.

### Widening the priors stalls the search

`sigma0_log` is both the prior width and the proposal bandwidth
(`adapt_frac_early * std_obs + (1 - adapt_frac_early) * sigma0_log`,
floored at `sigma_floor_frac * sigma0_log`). Multiplying it by 5
collapses proposal acceptance from ~0.25 to 0.02 and freezes the
objective at 3.0572 from generation 13 onward.

| sigma set | train | held out (A) | acceptance | objective |
|-----------|-------|--------------|------------|-----------|
| shipped | 0.9734 | 4.1260 | ~0.25 | still moving |
| shipped x5 | 3.0512 | 2.5452 | 0.02 | **frozen at gen 13** |

The x5 run scores better on held-out data, and that is not a
generalisation win: it barely moved, and on the one split where the
prior beats every tuned vector, barely moving is the winning play. This
is what criterion 5 exists to catch.

The measured uncertainties (BRENDA transfers across organisms at a
robust sd of 2.2 in log space, against the 0.2 asserted here) therefore
cannot be adopted as prior widths while the same number sets the
proposal bandwidth. Separating the two -- a bandwidth reference
independent of `sigma0_log` -- is the change that would let the priors
carry their real values without destroying acceptance.

## Open items

0. **Seed spread**, in flight: three seeds per mask at 15 generations.
   Everything ranked in this document rests on one seed each.
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
