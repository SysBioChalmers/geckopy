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

Measured by `screen.py`, which perturbs each reachable kcat by two-fold
up and two-fold down, keeps the larger absolute change in the distance,
and stores each perturbation's two dataset RMSEs separately so leverage
under any weighting follows without re-simulating. 3952 kcats screened;
23 minutes on 63 workers.

| leverage captured | kcats, flux x2 | kcats, max growth x2 |
|-------------------|----------------|----------------------|
| 50% | **11** | **12** |
| 83.7% / 84.3% | 110 | 110 |
| 90% | 535 | **239** |
| 99% | 3478 | 3293 |

The largest single kcat is 6967x the median under the old weighting and
12059x under the new. **Eleven or twelve kcats carry half the total
leverage in a 4834-parameter model**, and 110 carry 84%.

Under the objective that weights max growth double the distribution is
*more* concentrated, not less: 90% of the leverage sits in 239 kcats
against 535. Weighting the eight max-growth conditions more heavily
narrows the set of parameters the data can speak to.

These supersede the figures taken from `sensitivity.npy`, which read
39.3% for the top 10 and 50% for the top 24. That file's perturbation
step cannot be recovered -- no candidate step reproduces its values and
unrelated kcats share the identical entry 0.00275, so its low end looks
floored rather than measured. Its *ranking* was sound: it agrees with
the new screen on 94 of the top 110, and the two objectives agree with
each other at a rank correlation of 0.923, which is why the masks built
on it selected sensible kcats and why nothing measured on 3 September
needs revisiting on this account.

### The top of the distribution is enriched in `custom`

| source | share of top 110 by leverage | share of the model |
|--------|------------------------------|--------------------|
| custom | **22.7%** | 4.3% |
| brenda | 57.3% | 66.2% |
| okp | 17.3% | 24.3% |
| unlabelled | 2.7% | 5.3% |

`custom` is five-fold over-represented among the highest-leverage
kcats. The likely reason is selection -- a kcat gets curated because
someone found the model sensitive to it -- so this is not evidence that
curated values are wrong. It does mean the trust mask, which makes
`custom` clear a bar three times higher, buys part of its parsimony by
declining to touch high-leverage parameters. That is a weaker argument
than declining to touch parameters no condition can constrain.

## The top of the impact distribution is enriched in `custom`

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

## Held-out conditions: not part of this method's evaluation

**Decided 2026-09-04: do not split conditions into training and test
sets.** This model has 41 conditions and most adapters will have far
fewer -- eight, five, sometimes one growth rate. A validation scheme
that only works when conditions are plentiful cannot be part of how
this method is judged, and holding conditions back from a fit that is
already short of data makes the fit worse for no return.

Evaluate on all conditions, and rank variants on the quantities that
do not need a split: how many kcats a result changes, which ones, how
far, whether the trust ranking is respected, the share of screened
impact carried, and whether the objective stalled. Those are the five
criteria in the hand-over, and none of them requires holding data back.

The work already done in this direction is recorded below, both because
it produced a usable answer and because it produced a bug worth not
repeating. Nothing further is planned on it.

### What it showed before being dropped

Four runs, each fitted on part of a group split and scored on the
conditions it never saw:

| split | held out | untuned | tuned | reduction |
|-------|----------|---------|-------|-----------|
| A | 10 flux + 2 max-growth | 2.8353 | 1.5206 | 46% |
| B seed 0 | 15 flux + 2 max-growth | 13.8744 | 1.4756 | 89% |
| B seed 1 | 15 flux + 2 max-growth | 13.8744 | 1.5050 | 89% |
| C seed 0 | 6 flux + 2 max-growth | 5.8874 | 0.6099 | 90% |

Tuned kcats predict conditions they were not fitted to, in every split
tried. That is worth knowing and does not need re-establishing.

Note that split C scores better held out (0.6099) than in training
(1.0465): the two halves differ in difficulty, and the untuned model
scores 9.1434 against 5.8874 on them. A distance is only interpretable
against the same conditions' own baseline.

### The bug, so it is not repeated

`simulate_bayesian_dataset` blocks every exchange named as a
*condition* of the dataset it is handed, so one row's carbon source
cannot feed another row's simulation. The first split implementation
dropped rows, so a max-growth dataset holding only ethanol and acetate
left the six sugar uptakes open and the model grew on those instead.

| max-growth condition | vector | full 8-condition set | 2-condition subset |
|----------------------|--------|----------------------|--------------------|
| ethanol | prior | 2.1098 | 0.4379 |
| ethanol | reach@31 | **0.0683** | 11.9118 |
| acetate | prior | 5.3360 | 1.5733 |
| acetate | reach@31 | 3.7472 | 9.8942 |

The error reversed with the vector -- the prior is too constrained to
exploit the leaked sugars, a tuned vector is not -- which made it look
like a real effect and produced a confident, wrong conclusion that
tuning destroys growth prediction on non-fermentable carbon sources.
Scored correctly, tuning improves ethanol prediction thirty-fold.

Arithmetic caught it, not suspicion: per-condition values of 11.9 and
9.9 cannot sit inside `reach@31`'s aggregate of 0.8688, because
`dataset_rmse` takes a plain mean and max-growth carries weight 1 of 3.
**Reconcile per-condition detail against the reported aggregate
whenever a scoring path changes.** That check is free and would have
caught this four hours earlier.

The general lesson outlives the splits: subsetting a dataset changes
what the simulation blocks. Any future need to score a subset should
zero weights and keep every row, as `holdout.py` now does.

### Widening the prior widths: abandoned

**Decided 2026-09-04: do not widen `sigma0_log` beyond its shipped
values, coupled to the proposal or not.** The measurement that
motivated it stands -- a BRENDA kcat transferred across organisms is
off by nine-fold at one standard deviation, against the 0.20 the model
asserts -- but every attempt to act on it costs fit and returns
nothing.

All runs on the FVA mask, 15 generations, flux-weighted objective:

| widths | coupled? | distance | changed | median fold | acceptance | outcome |
|--------|----------|----------|---------|-------------|------------|---------|
| shipped 0.30/0.25/0.20/0.10 | - | **1.0003 +/- 0.050** | 3910 | 3.22 | ~0.25 | baseline, 3 seeds |
| shipped x5 | yes | 3.05, 2.70 | 3931 | 3.88 | 0.02 | stalled gen 13, gen 11 |
| measured 2.5/1.5/1.5/0.4 | yes | 3.4355 | 3935 | 6.52 | 0.02 | stalled gen 10, infeasible LPs |
| measured, proposal at shipped | no | 1.7397 | 3921 | - | 0.25 -> 0.10 | searches, still 74% worse |

Decoupling the proposal width from the prior width -- shipping
`proposal_sigma_log_default`/`_source` for it -- removes the stall and
recovers about two thirds of the deficit. It does not close it, and it
does not improve parsimony either: the wide runs change *more* kcats
than the baseline (3921-3935 against 3910) and move them further.

Two reasons the remaining gap is structural rather than a tuning
matter. The initial population is still drawn from `sigma0_log`
(`tuning.py:383`), so wide priors scatter generation 1 across orders of
magnitude and the search spends its budget recovering. And with
`prior_penalty_weight` at 0, `sigma0_log` never enters the objective at
all, so a wider prior cannot buy a better answer -- it can only change
where the sampler looks.

What the measurement is good for is reporting, not sampling. Movement
of 5.06-fold reads as eight prior sigma under the asserted widths and
0.74 sigma under the measured ones; the second is the honest figure and
belongs in how results are described. It does not belong in the search.

The decoupling stays in the package: it costs nothing when unset, and
it separates two meanings that were wrongly sharing one field.

## How many kcats does this data actually constrain?

Restricting the search by the sigma-weighted identifiability screen,
15 generations, seeds as marked:

| mask | eligible | changed | RMSE | sd | seeds |
|------|----------|---------|------|----|-------|
| none | 4834 | 4787 | 1.0778 | - | 1 |
| reach (FVA) | 3955 | 3910 | 1.0003 | 0.050 | 3 |
| trust 1e-4 | 2491 | 2458 | 0.9832 | - | 1 |
| trust 3e-4 | 1448 | 1437 | 0.9879 | 0.052 | 3 |
| trust 1e-3 | 111 | 110 | 1.0598 | 0.073 | 3 |
| **trust 3e-3** | **62** | **62** | **1.0852** | - | 2 |

Across a 77-fold range in how many kcats may move, fit varies by 0.12
against a seed spread of 0.05-0.07. Sixty-two kcats take the distance
from 8.6002 to 1.09; 3910 take it to 1.00. The gap between the tightest
and the loosest mask is roughly two standard errors -- suggestive, not
established -- so the honest statement is that **mask size is not a
lever on fit anywhere in the range tested**.

### Where the fit actually comes from

Take a tuned vector, revert all but its K most impactful changes, and
score. Ranking is by the one-at-a-time screen, which never saw either
result, so the ordering is not circular.

| top K kept | reach@31, 3921 changed | % of gain | trust 1e-3, 110 changed | % of gain |
|------------|------------------------|-----------|-------------------------|-----------|
| 0 (prior) | 8.6002 | 0% | 8.6002 | 0% |
| 1 | 4.8637 | **48.3%** | 5.0727 | 46.7% |
| 3 | 2.4485 | **79.6%** | 2.9512 | 74.8% |
| 5 | 1.9841 | 85.6% | 2.8246 | 76.4% |
| 10 | 1.8080 | 87.9% | 2.7670 | 77.2% |
| 20 | 1.4434 | 92.6% | 1.7586 | 90.6% |
| 40 | 1.1985 | 95.7% | 1.3116 | 96.5% |
| all | 0.8688 | 100% | 1.0446 | 100% |

**One kcat carries half the improvement. Three carry 80%. Forty carry
96%.** The other 3881 changes in `reach@31` buy the remaining 4.3%.

### The kcats in question

| rank | reaction | enzyme | source | kcat0 (1/s) | share |
|------|----------|--------|--------|-------------|-------|
| 1 | r_0698 | lanosterol synthase | brenda | **0.0019** | **16.9%** |
| 2 | r_0109 | acetyl-CoA carboxylase | brenda | 1.23 | 4.4% |
| 3 | r_0079 | phosphoribosylformylglycinamidine synthetase | brenda | 0.05 | 4.3% |
| 4 | r_0486_EXP_2 | (custom) | custom | 29 | 2.8% |
| 7 | r_1264 | succinate transport | okp | 12.3 | 2.0% |
| 9 | r_1239 | oxaloacetate transport | okp | 17.3 | 1.8% |
| 10-14 | r_0438_EXP_1..5 | five isozyme copies | okp | 54.3 | 5.8% |
| 19 | r_0015 | diaminopyrimidinone reductase | brenda | 0.00067 | 0.6% |

The head of the distribution is largely implausibly low BRENDA
assignments creating artificial bottlenecks -- a lanosterol synthase at
0.0019 1/s, a reductase at 0.00067 -- which every tuned vector raises
by one to two orders of magnitude.

The five isozyme copies of `r_0438` share a prior of 54.3 and receive
1.4x, 141x, 6.7x, 1.6x and 184x from `reach@31`. They are
interchangeable, so the search distributes values among them
arbitrarily. That is unidentifiability made visible, in a single row.

### What this means for the method

On this dataset the tuner's product is a few dozen kcat corrections,
dominated by a handful of bad database values. Everything downstream
follows from the objective being flat past that point: mask size does
not matter, seed spread (5%) exceeds every variant difference measured,
change counts reproduce to 0.2% while fit does not, the prior penalty
compressed magnitudes without reducing counts, and identical isozymes
receive wildly different values.

This does not make the method worthless -- it found the corrections,
and the screen that ranks them is part of the same machinery. It does
mean the result should be *reported* as those corrections, with their
provenance, rather than as a 4834-long posterior. A reviewer can check
whether lanosterol synthase is really 0.0019 1/s. Nobody can check
3921 simultaneous changes.

## Re-baselining under max growth weighted double

`max_growth_weight` was changed on 2026-09-04 to weight the max-growth
term, so 2 makes the eight max-growth conditions count double against
the 33 flux conditions. Every result before that date used the opposite
emphasis, which is 0.5 under the new convention. The untuned model
scores 8.4043 under the new objective against 8.6002 under the old.

Masks for these runs come from `screen_components.npy` evaluated at the
same weighting, so the mask selects for the objective being minimised.
Thresholds were chosen to admit the same mask sizes as the 3 September
curve.

### Re-optimising for the new objective gains nothing

| mask | fitted to old, scored old | fitted to old, scored new | fitted to new, scored new |
|------|---------------------------|---------------------------|---------------------------|
| reach, 31 gen | 0.8688 | **0.9727** | **0.9822** |
| trust, 31 gen | 0.8844 | **0.9501** | **0.9891** |

A vector fitted to the old objective scores *better* on the new
objective than one fitted to the new objective does, for both masks.
The differences -- 0.010 and 0.039 -- sit inside the 0.05 seed spread,
so the honest statement is that **65 minutes of re-optimisation
recovers nothing over simply rescoring the old result.** The weighting
is not a lever on what the tuner produces.

That is consistent with the screen: leverage under the two weightings
correlates at 0.923, so both objectives ask about substantially the
same parameters and both find them.

### The mask curve under the new objective

Fourteen runs, masks rebuilt from `screen_components.npy` at the same
weighting, thresholds chosen to admit the sizes the 3 September curve
used. Untuned distance is 8.4043.

| mask | eligible | changed | distance | sd | seeds | impact |
|------|----------|---------|----------|----|-------|--------|
| none | 4834 | 4781 | 1.0641 | - | 1 | 99.7% |
| reach | 3955 | 3924 | 1.0654 | 0.0804 | 3 | 98.6% |
| trust 3e-5 | 1461 | 1446 | 1.0828 | 0.0697 | 3 | 95.4% |
| **trust 9.36e-4** | **112** | **110** | **1.0422** | **0.0569** | 3 | 82.5% |
| trust 3.07e-3 | 63 | 62 | 1.1736 | 0.0204 | 2 | 77.2% |
| reach, 31 gen | 3955 | 3925 | 0.9822 | - | 1 | 98.8% |
| trust 3e-5, 31 gen | 1461 | 1451 | 0.9891 | - | 1 | 95.4% |

**Restricting the search to 112 of 4834 kcats gives the best mean
distance on the curve**, and the tightest seed spread with it. Removing
the restriction entirely is no better than the FVA mask, and both are
behind the 112-kcat mask by about a quarter of a standard deviation --
which is to say all three are the same number.

The curve is flat from 4781 changes down to 110 and bends at 62, where
the distance rises 0.107 above the FVA mask at roughly 1.3 standard
errors. The old objective bent in the same place, so two objectives
independently locate the useful mask size between about 60 and 110
kcats.

Sixteen extra generations buy 0.08 (reach) and 0.09 (trust 3e-5),
about one standard deviation, and spend it by pushing the same
parameters three times further: `reach` at 15 generations moves the
four sources by 1.72 / 3.05 / 3.30 / 4.55 fold, and at 31 by 2.41 /
4.70 / 6.71 / 14.50.

### Report the trust-ordering check as inconclusive on small samples

`source_movement`'s ordering test compares four medians, and at the
tight masks those medians are taken over very few kcats. At the
112-kcat mask the four sources contribute 15 / 62 / 24 / 5 moved kcats;
at the 63-kcat mask, 7 / 35 / 13 / 3. Three values cannot support a
median, and the verdicts behave accordingly: the 112-kcat mask passes
at seeds 0 and 2 and fails at seed 1, and the 63-kcat mask fails at
both seeds while its `unlabelled` median swings between 2.85 and 9.78.

The check is sound where it was designed to work -- every mask down to
1446 changed kcats passes it at every seed -- and meaningless below
roughly ten moved kcats in any source. Report it as inconclusive there
rather than as a pass or a fail, and read `n_moved` before reading the
verdict.

## A better optimiser fits better and identifies nothing

ABC-SMC draws and truncates; it is not an optimiser, and once a screen
has reduced the problem to about a hundred parameters it is the wrong
instrument. Separable CMA-ES in log-space over the same-sized parameter
set, three seeds, 6800 evaluations each, 45 minutes per seed on 63
workers.

The parameter set was chosen with tying decided first: isozyme copies
sharing a reaction, prior and source are one parameter, the screen
moved each group as a unit, and the mask admitted 112 free parameters
covering 144 kcats across **100 distinct reactions**, against 79 for
the same budget spent per copy.

| | distance | sd | seeds |
|---|----------|----|----|
| ABC-SMC, 112 kcats, 31 generations | 0.9156 | 0.0549 | 3 |
| **CMA-ES, 112 free parameters** | **0.7923** | **0.0071** | 3 |

**The optimiser wins by 0.1233 +/- 0.0320 -- 3.9 standard errors -- and
its run-to-run spread is eight times tighter.** So 0.9156 was the
sampler's floor, not the objective's, and the last step of a tuning
pipeline should be an optimiser rather than a sampler.

### And yet it pins nothing

| agreement across three seeds | of 144 kcats moved |
|------------------------------|--------------------|
| within 1.2x | **0** |
| within 2x | 9 |
| beyond 5x | **97** |

Three runs agree on the distance to 0.7% and disagree by more than
five-fold on two thirds of the kcats they change. The sampler managed
0 of 112 within 1.2x and 70 beyond 5x -- proportionally the same.

This settles a question left open earlier. The sampler's parameter
scatter might have been search noise that a better-converged method
would resolve; it is not. A method converging to 0.7% on the objective
scatters the parameters just as widely, which means **these are flat
directions in the objective, not failures of search.** Many kcat sets
fit these 41 conditions equally well and nothing in the data chooses
between them.

The practical consequence: optimising harder buys *fit*, never
*knowledge*. A tuned kcat vector is a model that predicts better, not a
set of measurements that can be published as corrections. What survives
across methods and seeds is which parameters matter -- the screen's
ranking -- and that comes without any optimisation at all.

### The comparison is not yet like-for-like

CMA-ES ran unbounded where the sampler clips proposals to
`kcat_lo`/`kcat_hi`. Its largest single changes are 3523x, 28039x and
1357x across the three seeds, which no biochemist would accept and
which the sampler was structurally prevented from producing. Some of
the 0.1233 advantage is therefore freedom rather than skill. Repeating
it with the same bounds would settle how much; it is one small change
and one run.

### The proposal floor, and how much it actually mattered

`kcat_bounds` floored proposals at 1e-2 1/s and widened only upwards,
so any prior below the floor sat outside its own window: 350 of 4834
kcats on this model, clipped up before scoring, a 62-fold increase for
the slowest. A prior must be a proposable value, and this one was not.
Fixed in `01b6ec9`; MATLAB's `proposeSimple` has the same defect.

It is worth being precise about the damage, because the first estimate
was too strong. The floor was **not binding at the optimum**:

| lanosterol synthase | value | vs prior |
|---------------------|-------|----------|
| prior | 0.00194 | - |
| lowest the old floor allowed | 0.01 | 5.2x |
| CMA-ES optimum, unbounded | 0.900 | 465x |
| CMA-ES optimum, corrected bounds | 0.338 | 174x |

The search drives this kcat two orders of magnitude past the floor, so
the floor never determined its answer. The claim that it contaminated
every magnitude reported here is withdrawn; what it did was start every
particle 5x above the prior for those 350 kcats, which is a defect in
its own right without changing where the search ended up.

The two runs disagreeing by 2.7-fold on that same parameter -- 0.900
against 0.338, at distances of 0.8005 and 0.8011 -- is the flat-
directions result again, now visible on the single most consequential
kcat in the model.

### The optimiser's advantage is not bounds-freedom

| | distance | changes beyond 100x |
|---|----------|---------------------|
| CMA-ES, unbounded | 0.8005 | 8 |
| CMA-ES, corrected bounds | 0.8011 | 5 |
| ABC-SMC | 0.9156 +/- 0.0549 | - |

Clipping the optimiser to the same window the sampler uses costs
0.0006. The ~0.12 gap to the sampler is searching, not licence.

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
