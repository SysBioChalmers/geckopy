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
   with `objective_trace` alongside.
2. **Number of kcats changed** (>2% from prior).
3. **Median fold change per source**, most-trusted first, and whether
   the trust ranking is respected. Report `median_fold_moved` with
   `n_moved` beside it: a median over the whole source saturates at
   1.00 once a mask leaves most of it at prior, and then the ordering
   check passes vacuously. Do *not* report movement in sigma units:
   `sigma0_log` differs by source, so it divides out the tier being
   examined and has already inverted one conclusion here. Only the
   ranking of the sigmas is meaningful.
4. **Impact share** -- the fraction of total screened `|dRMSE|` carried
   by the changed kcats. Distinguishes "few and consequential" from
   merely few. The blend changes 560 kcats but carries 28% of the
   available impact; 301 screen-chosen changes carry 69%.
5. **Whether the objective stalled.** A frozen `objective_trace` with
   collapsing acceptance produces low movement for the wrong reason.
   Two runs in the penalty sweep looked parsimonious because they had
   stopped searching.

`parsimony.py` computes 2-4 from saved `.npy` vectors, MATLAB's export
included, without re-solving.

## The target shape

Few changes, each large, each load-bearing. Specifically **not** many
kcats nudged slightly, and **not** few changes to parameters the data
cannot see. Current best answers:

| operating point | RMSE | changed | median fold | impact |
|-----------------|------|---------|-------------|--------|
| best fit (FVA mask) | 0.8688 | 3921 | 5.06 | 99.4% |
| **trust mask** | 0.8844 | **1433** | 3.69 | 88.5% |
| balanced (FVA, pruned 1e-3) | 0.9255 | 1920 | 4.98 | 94.5% |
| few and large (FVA, pruned 3e-3) | 1.0124 | **301** | 4.48 | 69.1% |

All four beat MATLAB on parsimony; the first beats it on fit. The
trust mask is off the others' trade-off curve: fewer changes and a
better fit than the balanced point, at a slightly lower impact share.

## Recommended method

FVA-reachable mask, then prune the result by the sensitivity screen.
`prior_penalty_weight` stays 0. No hyperparameter is added: the mask is
derived from the model and conditions, and the pruning threshold is an
operating point on a reported curve, not something tuned per model.

The trust-weighted mask is the likely replacement and beats it on the
numbers, but it carries a threshold of its own and rests on one seed
and one threshold value. Adopt it once the seed sweep lands and a
second threshold is on the curve.

Two mechanisms exist but are **not** recommended:

* `prior_penalty_weight` compresses how far kcats move rather than how
  many, which is the opposite of the target shape (median fold 2.58,
  only 61% beyond two-fold). Dropped.
* `adapt_proposal_width` gives the best raw RMSE of any variant
  (0.8281) but scatters parameters badly. Off by default; do not enable
  without reporting all five criteria above.

## Known open problems

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
3. **Seed spread.** Three seeds per mask at 15 generations are running
   (`logs/seeds_latest.log`, ~2.5 h, `sweep_seeds.sh`). Until they
   land, every ranking here rests on one seed per variant and
   differences of a few percent -- the reach-to-trust gap included --
   are not safe.
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
