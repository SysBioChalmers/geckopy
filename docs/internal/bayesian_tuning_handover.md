# Hand-over: Bayesian kcat tuning, current state

Where the work stands, what is settled, and what to do next. Detailed
numbers are in `matlab_replication_results.md`; the method's design
history is in `bayesian_tuning_plan.md`. This file is the entry point.

Branch: `feat/bayesian-kcat-tuning` (pushed). The earlier
`feat/bayesian-matlab-fidelity` is fully merged into it and can be
deleted; `origin/feat/bayesian-kcat-tuning-agent` is stale and carries
a superseded conclusion.

## What is settled

**The port reproduces MATLAB's scoring to 0.002%.** Two pieces of
`abc_max.m` were missing, and each one *alone* makes agreement worse,
which is why they were hard to find: the anaerobic switch (9 of the 33
flux conditions have zero measured oxygen; `abc_max.m:81-84` swaps the
model for those) and `max_growth_weight = 2` (`abc_max.m:48` weights
the *flux* term, despite the name). Verified against MATLAB's own
exported prior and accepted particles.

**The sampler is faithful.** The samples/accepted sequence reproduces
MATLAB's generation for generation, down to the terminal 171/571 fixed
point. Schedule, pooling and truncation selection need no further work.

**The reference target is 8.6001 -> 0.8715**, from
`GECKO/tutorials/full_ecModel/output/rmse_trace.tsv`. The widely quoted
9.5544 -> 0.9289 predates the June 2026 model update and describes a
different model; do not compare against it.

**Most kcat changes were never justifiable.** FVA across the 41
conditions leaves 879 kcats unable to carry flux in any of them, and
every unrestricted run changed ~873 of those. Nothing in the method
asked whether a parameter was identifiable.

**Best fit to date**: FVA-restricted, 31 generations, RMSE **0.8688**
with 3921 kcats changed, against MATLAB's 0.8715 with ~4804. Better on
both axes, from a mask with no hyperparameters.

**Best result to date**: the trust-weighted mask
(`identifiability_mask` at 3e-4, 1448 eligible), 31 generations, RMSE
**0.8844** with **1433** changed and 88.5% impact share. It beats the
pruned operating point below on fit and on count at once (1920 changed
at 0.9255), while carrying 6 points less of the screened impact for
25% fewer changes, and it is the first run that leaves trusted kcats
alone -- 84% of `custom` within 1.1x of
prior, against 10% under the FVA mask. Single seed; see open problem 3.

## How to judge any new strategy

Never on RMSE alone. Report all of these, every time:

1. **Plain RMSE**, comparable to MATLAB's 0.8715. Keep it plain even
   when selection optimises something else -- `rmse_trace` does this,
   with `objective_trace` alongside. One seed cannot support a fit
   claim: spread between seeds is 5%. Report three, or report no fit
   ranking.
2. **Number of kcats changed** (>2% from prior).
3. **Median fold change per source**, most-trusted first, and whether
   the trust ranking is respected. Report `median_fold_moved` with
   `n_moved` beside it: a median over the whole source saturates at
   1.00 once a mask leaves most of it at prior, and then the ordering
   check passes vacuously. Below about ten moved kcats in any source
   the ordering verdict is noise -- at the 63-kcat mask it is computed
   over 7 / 35 / 13 / 3 and flips between seeds -- so report it as
   inconclusive there rather than as a pass or fail.

   Never report movement in sigma units. Dividing by `sigma0_log`
   removes exactly the thing being examined, and it has now inverted
   two readings: across sources, where sigma differs by tier, and
   across configurations, where `mean_dev` made runs with 5x wider
   priors look conservative (1.32 sigma against 5.36) while they were
   in fact moving kcats further (median 3.88-fold against 3.22, and
   6.52-fold for the measured widths). `mean_dev` is comparable only
   within one sigma setting. Fold changes are comparable everywhere.
4. **Impact share** -- the fraction of total screened `|dRMSE|` carried
   by the changed kcats. Distinguishes "few and consequential" from
   merely few. The blend changes 560 kcats but carries 28% of the
   available impact; 301 screen-chosen changes carry 69%.
5. **Whether the objective stalled.** A frozen `objective_trace` with
   collapsing acceptance produces low movement for the wrong reason.
   Two runs in the penalty sweep looked parsimonious because they had
   stopped searching.

6. **What it did to conditions it already fit.** An aggregate gain can
   hide a per-condition regression, so score the prior and the result
   on the same conditions and report both. Subset a dataset by zeroing
   weights, never by dropping rows: dropping a row unblocks that
   condition's carbon source for every other row.

`parsimony.py` computes 2-4 from saved `.npy` vectors, MATLAB's export
included, without re-solving. Leverage comes from
`screen_components.npy` (written by `screen.py`), which stores each
perturbation's two dataset RMSEs so any weighting is arithmetic. The
older `sensitivity.npy` is superseded: its perturbation step cannot be
reproduced, though its ranking agrees on 94 of the top 110. `holdout.py` and `holdout_baseline.py` in
the run scratch do 6, and `run_sigma.py` reports train, held-out and
all-41 distances for a run trained on one half of a split.

## Two methods, one conclusion

CMA-ES over ~112 screened parameters reaches **0.7923 +/- 0.0071**
against ABC-SMC's 0.9156 +/- 0.0549 -- better by 3.9 standard errors,
with eight times tighter run-to-run spread. The last step of a pipeline
should be an optimiser, not a sampler.

It identifies nothing the sampler did not. Across three seeds agreeing
to 0.7% on the distance, **0 of 144 changed kcats agree within 1.2x and
97 disagree beyond 5x**. These are flat directions in the objective,
not search failures, so optimising harder buys fit and never knowledge.
(CMA-ES ran unbounded and produced changes up to 28039x; repeat with
`kcat_lo`/`kcat_hi` before quoting the advantage.)

## The finding that reframes the rest

Reverting a tuned vector to its K most impactful changes: **1 kcat
gives 48% of the improvement, 3 give 80%, 40 give 96%**, and the
remaining 3881 changes buy 4.3%. Masks admitting 62 to 4834 kcats all
reach the same distance within seed noise. The screen agrees from the
other direction: **11 kcats carry half the total leverage** and 110
carry 84%.

The head of the impact distribution is a short list of implausible
database values -- lanosterol synthase at 0.0019 1/s carries 16.9% of
all achievable improvement on its own. The tuner's real product on this
dataset is a few dozen corrections, and it should be reported that way,
with provenance, rather than as a 4834-long posterior.

## The target shape

Few changes, each large, each load-bearing. Specifically **not** many
kcats nudged slightly, and **not** few changes to parameters the data
cannot see. Current best answers:

Under the objective that weights max growth double (2026-09-04), the
same curve reads: unrestricted 1.0641 with 4781 changed, FVA mask
1.0654 +/- 0.080 with 3924, trust 3e-5 1.0828 +/- 0.070 with 1446,
**trust 9.36e-4 1.0422 +/- 0.057 with 110**, trust 3.07e-3 1.1736 with
62. The 112-kcat mask has the best mean and the tightest spread.

| operating point | RMSE | changed | median fold | impact |
|-----------------|------|---------|-------------|--------|
| best fit (FVA mask) | 0.8688 | 3921 | 5.06 | 99.4% |
| **trust mask** | 0.8844 | **1433** | 3.69 | 88.5% |
| balanced (FVA, pruned 1e-3) | 0.9255 | 1920 | 4.98 | 94.5% |
| few and large (FVA, pruned 3e-3) | 1.0124 | **301** | 4.48 | 69.1% |

All four beat MATLAB on parsimony. Read the RMSE column as one draw
from a distribution 5% wide: these rows are separated by parsimony,
not by fit, and the trust mask changes the fewest kcats of any row that
keeps impact share near 90%.

## Recommended method

FVA-reachable mask, then prune the result by the sensitivity screen.
`prior_penalty_weight` stays 0. No hyperparameter is added: the mask is
derived from the model and conditions, and the pruning threshold is an
operating point on a reported curve, not something tuned per model.

The trust-weighted mask is the likely replacement and beats it on the
numbers, but it carries a threshold of its own and rests on one seed
and one threshold value. Adopt it once the seed sweep lands and a
second threshold is on the curve.

Two directions are closed and should not be reopened: **train/test
splits over conditions** (most adapters have far fewer than 41
conditions) and **widening the prior widths** (stalls when coupled,
loses when decoupled).

Two mechanisms exist but are **not** recommended:

* `prior_penalty_weight` compresses how far kcats move rather than how
  many, which is the opposite of the target shape (median fold 2.58,
  only 61% beyond two-fold). Dropped.
* `adapt_proposal_width` gives the best raw RMSE of any variant
  (0.8281) but scatters parameters badly. Off by default; do not enable
  without reporting all five criteria above.
* **Widening `sigma0_log` is closed.** Coupled to the proposal it
  stalls the sampler (acceptance 0.02, two seeds at x5, one at the
  measured widths); decoupled it searches but still lands 74% worse
  than the shipped widths while changing more kcats and moving them
  further. The measured spreads are for *reporting* movement honestly,
  not for sampling. Do not reopen this.

## Known open problems

0. **Do not split conditions into train and test.** Decided
   2026-09-04: this model has 41 conditions and most adapters will have
   far fewer, so a scheme needing spare conditions cannot be how the
   method is judged, and withholding data from an already thin fit
   costs more than it returns. Fit on everything; rank on the five
   criteria below, none of which needs a split. Four runs were done
   before this was decided and they answer the question they were
   asked -- tuned kcats predict conditions they were not fitted to, by
   46% to 90% -- so it does not need revisiting.

1. **The trust mask needs a second threshold and a seed.** 3e-4 is the
   only value run. Sweeping it traces mask size against fit the way
   the pruning threshold traces change count, and turns the
   hyperparameter into a reported curve. `custom` movement is largely
   answered by it -- 84% within 1.1x of prior, 37 of 207 moved -- but
   the 37 that move still move 2.49-fold at the median.
2. **The FVA mask's ordering failure is in counts, not magnitudes.**
   Among moved kcats it respects the trust ranking; it changes 94.7%
   of `custom` against 68.9% of `okp`, which is where being
   source-blind shows. The trust mask addresses exactly this.
3. **Seed spread is measured, and it invalidates the fit rankings.**
   sd 0.050 (reach) and 0.052 (trust) over three seeds at 15
   generations; the masks differ by 0.0125 +/- 0.0418. The port-vs-
   MATLAB gap (4.3%), the B1 advantage (8.9%) and the whole penalty
   sweep all sit inside that. Whether spread narrows at 31 generations
   is unmeasured -- that is the cheap follow-up.
4. **The blend is never scored or optimised**, in either
   implementation, and does not survive an impact-weighted reading.
   Whether to keep reporting it at all is an open question.
5. **The proposal kernel has no move that removes a change.** Sparsity
   is unreachable by search, only by masking and pruning. Adding a
   reversion move -- snap a random subset back to prior -- is the
   obvious structural fix and needs no per-model tuning.
6. **Performance.** `model.optimize()` builds every primal, reduced
   cost and shadow price when under 40 numbers are read, ~2.0 s per
   particle; `apply_kcat_constraints` is another 27%, re-deriving
   caches and issuing ~9700 `add_metabolites` calls per particle.
   Neither changes any result, so both are verifiable by reproducing a
   trace exactly.

## Running it

See `matlab_replication_handoff.md` for the cluster environment: venv
on node-local disk, Gurobi `Threads 1` via a `gurobi.env` file (setting
`Params.Threads` does not reach the pool workers), `n_proc` left unset.
A 31-generation run is ~65 min on 63 workers.

Scripts used for the runs above live in the run scratch at
`/cephyr/NOBACKUP/groups/compmeteng/geckopy-matlab-repl-64/`:
`run_restricted.py` (masks, seeds via `MR64_SEED`, all five criteria
reported), `sweep_seeds.sh`, `parsimony.py`, `folds.py`, `settle.py`.
The venv lives on node-local disk and dies with the allocation that
built it; `build_venv_20.sh` rebuilds it in ~3 min. They are working scripts, not part of the
package; the reusable parts are in
`geckopy.kcat_sensitivity_analysis.bayesian.parsimony`.
