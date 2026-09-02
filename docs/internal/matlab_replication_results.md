# Results: reproducing MATLAB's Bayesian kcat tuning in Python

The run specified by `matlab_replication_handoff.md`, carried to
MATLAB's full 31 generations. Read that hand-off first for the setup;
this file records what came out and what it implies.

Headline: the port reaches **RMSE 9.6229 -> 1.2038 in 31 generations**,
an 87.5% reduction, against GECKO's own exported run's 8.6001 ->
0.8715 (89.9%). The search machinery reproduces MATLAB's exactly,
generation for generation. The remaining difference is in the distance
function, not the sampler -- and the MATLAB reference numbers are
themselves inconsistent between runs, so "the target" is not currently
a single well-defined value.

## Provenance

* Code: `feat/bayesian-kcat-tuning-agent` at `b26a9b8`.
* Model: `ecYeastGEM_preTune_GECKOderived.yml` (8001 rxns, 3893 mets,
  4834 tunable kcats, `bio_rxn = r_4041`).
* Hyperparameters: `YeastGEMAdapter.m` verbatim, as listed in the
  hand-off; `selection="truncation"`, `seed=0`, `max_generations=31`.
* Node: C3SE Vera `vera-r01-19`, 64 cores, `n_proc=63` (left unset,
  resolved by `cobra.Configuration().processes`), Gurobi `Threads=1`
  per worker.
* Wall time: **72.0 min**. For scale, the hand-off records 178.5 min
  for 16 generations on 16 cores.
* Artifacts:
  `/cephyr/NOBACKUP/groups/compmeteng/geckopy-matlab-repl-64/`
  (`mr64_{new_kcat,old_kcat,blend,groups}.npy`, `mr64_result.json`,
  `reweight.py`, `logs/`).

## RMSE trace

Two MATLAB references are quoted throughout, because they disagree
(see "The reference target is not a single number"):

* **export** -- `GECKO/tutorials/full_ecModel/output/rmse_trace.tsv`
  with `run_metadata.txt`, a complete per-generation trace.
* **hand-off** -- the numbers recorded in
  `matlab_replication_handoff.md`.

| Gen | This run | MATLAB export | MATLAB hand-off |
|-----|----------|---------------|-----------------|
| prior | 9.6229 | 8.6001 | 9.5544 |
| 2   | 4.3621 | 4.1256 | 4.7157 |
| 4   | 3.3345 | 2.9219 | 3.4267 |
| 6   | 2.0517 | 1.9415 | 2.4551 |
| 8   | 1.7317 | 1.3790 | 1.6337 |
| 10  | 1.4575 | 1.2807 | 1.2842 |
| 14  | 1.3351 | 0.9933 | 0.9730 |
| 22  | 1.2775 | 0.8715 | 0.9289 |
| 31  | **1.2038** | 0.8715 | -- |

Full trace: 4.9865, 4.3621, 3.5495, 3.3345, 2.6274, 2.0517, 1.7317,
1.7317, 1.6092, 1.4575, 1.4575, 1.4254, 1.3351 x3, 1.3033 x5,
1.2775 x2, 1.2761, 1.2038 x8. `converged=False`; the only convergence
test is `rmse_threshold=0.2`, far out of reach.

## Per-source movement

Fraction of kcats changed by more than 2%:

| Source | n | Best particle | Blend | MATLAB reports |
|--------|-----|---------------|-------|----------------|
| brenda | 3198 | 99.1% | 1.8% | 2.1% |
| custom | 207 | 97.6% | 0.0% | 0.0% |
| okp | 1175 | 99.4% | 50.5% | 15.4% |
| unlabelled | 254 | 99.6% | 30.7% | 35.0% |

This answers the open question in the hand-off's reporting item 4:
**MATLAB's published percentages describe the blend, not the best
particle.** brenda, custom and unlabelled land close (1.8 vs 2.1, 0.0
vs 0.0, 30.7 vs 35.0). Only OKP is far off, and OKP is also the group
whose classification differs most between adapters, so it is the one
to look at first. The best particle moves essentially every kcat, as
`bayesian_tuning_plan.md` predicts.

## The blend

`posterior_trace[-1].kcats` scores **4.1822**, against the best
particle's 1.2038 and the prior's 9.6229. So the vector that looks
like a defensible posterior -- sparse, shrunk toward the prior -- is
a much worse fit than the vector the method actually returns. MATLAB
reports this vector's statistics and never scores it. The tension is
real and unresolved; it is not an artifact of the port.

## The search machinery is exactly faithful

MATLAB's `diagnostics_pergeneration.tsv` gives samples/accepted per
generation as 1000/300, 800/330, 800/339, 800/341, 800/342 x4,
600/282, 600/264, 600/259, 600/257 x3, 400/197, 400/179, 400/173,
400/171... This run reproduces that sequence generation for
generation, down to the terminal 171/571 fixed point and MATLAB's
recorded `n_accepted_samples 171`.

Schedule, pooling of new proposals with the previous accepted set, and
truncation selection are therefore identical. Verified by reading as
well: `targetAccept = 10` is below `minKeep = 0.3`, so minKeep always
binds and truncation to the best 30% is equivalent to MATLAB's
percentile-then-clamp; the proposal-width formula and its constants,
the parent-spreading step, the `proposeSimple` bounds rule and the
generation-1 prior's `-0.5*sigma^2` mean correction all match.

## The reference target is not a single number

The two MATLAB references report prior RMSEs of **9.5544** and
**8.6001** -- an 11% spread on a quantity that is deterministic given
a fixed model, dataset and distance function. The port gives 9.6229.

The hand-off's claim that the port's distance is faithful to 0.7%
rested on comparing against the 9.5544 figure alone. Against the
export's 8.6001 the same comparison is 12% out. **The distance
function is therefore not verified**, and the hand-off's conclusion
that "a shortfall in the final RMSE is a search problem, not a scoring
problem" does not hold. Given the machinery result above, the reverse
is closer to the truth.

Reconciling the prior RMSE against GECKO's export is the prerequisite
for treating any final number as a target.

## `maxGrowthWeight` is unimplemented, and does not explain the gap

`YeastGEMAdapter.m:95` sets `maxGrowthWeight = 2`, and `abc_max.m:48`
combines the two dataset terms as
`(maxGrowthWeight*rmse_flux + rmse_maxgrate) / (maxGrowthWeight + 1)`
-- note the weight named for max-growth multiplies the *flux* term.
geckopy has no such parameter; `distance.py` takes a plain mean.

That is a genuine gap in the port, but applying the 2:1 weighting to
this run's vectors makes agreement with MATLAB *worse*, not better:

| vector | flux F | maxgrow M | port 1:1 | MATLAB 2:1 |
|--------|--------|-----------|----------|------------|
| prior | 11.0374 | 8.2083 | 9.6229 | 10.0944 |
| best particle | 1.4123 | 0.9953 | 1.2038 | 1.2733 |
| blend | 3.2195 | 5.1450 | 4.1822 | 3.8613 |

The 2:1 prior (10.0944) is further from both 9.5544 and 8.6001 than
the plain mean (9.6229). Implementing the parameter is defensible for
fidelity to current GECKO, but it should not be expected to close the
gap, and it is not the explanation for it.

The per-condition `bayesianRMSEweight` column is a separate mechanism
and is already implemented faithfully (`distance.py`, matching
`abc_max.m`'s multiply-then-mean and its NaN -> 99 penalty). It is
uniformly 1.0 across all 33 flux and 8 max-growth conditions in this
dataset, so it is currently inert.

## Why both implementations plateau

MATLAB's `proposalWidthTrace` *grows* monotonically over the run,
0.213 -> 1.037, while `diversityTrace` grows 1.23 -> 9.08.
`updateProposalWidth` is `0.5*std_obs + 0.5*sigma0log` floored at
`0.15*sigma0log`: because the blend is 50/50 against the fixed prior
width, the proposal width can never fall below `0.5*sigma0log` -- the
floor is unreachable -- and it inflates as the accepted set spreads.

Across 4834 parameters this means every proposal displaces the whole
vector by roughly `0.5*sigma0*sqrt(4834)` in log space. No proposal is
ever a local refinement, so late-stage progress depends on a lucky
global draw and the incumbent best survives untouched. That is the
mechanism behind this run's flat generations (16-20 and 24-31 are
identical) and behind MATLAB's own plateau from generation 22. It is
shared by both implementations and explains the shape of the curve,
not the difference in level.

Annealing the proposal width is the obvious lever, and it is a
deliberate departure from MATLAB rather than a fidelity fix.

## Open items

1. Reconcile the prior RMSE (9.6229 here, 8.6001 and 9.5544 in the two
   MATLAB references) before treating any final RMSE as a target.
2. Establish the seed-to-seed spread of the port's final RMSE. This
   run is a single seed; the same seed at `n_proc` 16 and 63 reproduces
   identically, so within-seed determinism holds, but across-seed
   variance is unmeasured and bounds how much of any gap is real.
3. Decide whether to implement `max_growth_weight` in `BayesianParams`
   for fidelity to current GECKO, knowing it moves the numbers away
   from both references.
4. Investigate OKP's blend movement (50.5% here vs MATLAB's 15.4%),
   the one per-source figure that does not match.
