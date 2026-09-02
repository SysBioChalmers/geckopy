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

## Two diagnostics on the distance function

**The objective is carried by a handful of conditions.** Per-condition
prior RMSE across the 33 flux conditions spans 0.04 to 39.01 -- three
orders of magnitude. Six conditions (39.0, 28.3, 26.0, 25.8, 22.5,
21.5) supply over half the mean; the ten best supply under 2% of it.
Because the datasets are combined as plain means of raw RMSEs, the
search is effectively fitting those six conditions, and any small
implementation difference in how *they* simulate moves the reported
score a lot. Condition 10 alone contributes 0.59 to the 9.6229 total.
This is the first place to look for the 11% prior disagreement.

**The missing anaerobic switch is a real defect but not the
explanation.** `abc_max.m:81-84` calls
`modelAdapter.makeModelAnaerobic` for every condition whose oxygen
exchange is 0; 9 of the 33 flux conditions qualify, and the port runs
all 9 aerobically. That is 27% of the flux dataset simulated under the
wrong physiology. It does not, however, show up as inflated RMSE:
those 9 average **10.50** against the aerobic conditions' **11.24**,
and excluding them *raises* the combined score from 9.6229 to 9.7229.
Three of them (0.72, 1.05, 1.25) score suspiciously well, which is
what an unconstrained oxygen supply would produce. So the hook must be
implemented for fidelity, but it cannot be assumed to close the gap.

`changeProteinBiomass` (`abc_max.m:87-88`) stays inert: `Ptot` is
`NaN` for every condition in this dataset.

## Where to go from here

Two tracks. They are not exclusive, but only the first can be measured
against MATLAB, and the second cannot be evaluated at all until the
first has fixed a trustworthy score.

### Track A -- finish the replication (fidelity)

**A1. Localise the prior-RMSE disagreement, per condition. Do this
first; it blocks everything else.** `abc_max` already builds a named
`rmseList` (`fluxData_1..33`, `maxGrowth_1..8`) and throws it away.
Instrument it to dump that list for the prior kcat vector, run it once
-- a single evaluation, seconds, no tuning loop -- and diff against the
port's vector, which is reproduced by `percond.py` in the run scratch.
Given the spread above, the disagreement will localise to a few
conditions, and each one is then a concrete simulation question
(bounds, carbon normalisation, blocked-flux handling) rather than a
diffuse 11%. Until this lands, no final RMSE from either
implementation is interpretable, and the ranking of any method change
is unreliable.

**A2. Implement the anaerobic hook.** Needed for correctness whatever
A1 finds: 9 conditions are currently scored under the wrong
physiology. It belongs on the tutorial adapter rather than in
`geckopy`, which is where its organism-specific knowledge lives.

**A3. Settle which MATLAB reference is authoritative.** Re-run
`bayesianSensitivityTuning` once with the committed
`YeastGEMAdapter.m` and export the trace, so the target is one number
with known provenance instead of two that differ by 11%.

**A4. Implement `max_growth_weight`** in `BayesianParams`, defaulting
to 1, for fidelity to current GECKO. Cheap, but note it moves the
numbers *away* from both references, so land it after A1 and A3 or it
will confound them.

### Track B -- improve the method (performance)

The port is faithful; the algorithm it faithfully reproduces is the
weak part. In descending order of expected value:

**B1. Anneal the proposal width.** This is the single largest lever.
`updateProposalWidth` is `0.5*std_obs + 0.5*sigma0log`, so the width
is pinned above `0.5*sigma0log` and, because `std_obs` grows as the
accepted set spreads, it *inflates*: MATLAB's own
`proposalWidthTrace` goes 0.213 -> 1.037 while `diversityTrace` goes
1.23 -> 9.08. The sampler diffuses rather than converges, which is why
both implementations plateau and why late generations produce
identical best particles for eight generations at a stretch. Replace
the fixed 50/50 blend with a schedule that decays toward `std_obs`, or
scale it adaptively against `proposalAcceptRate` -- MATLAB records
that signal (0.10 at generation 1, 0.14 at generation 2) and never
uses it.

**B2. Rebalance the per-condition contributions.** With per-condition
RMSE spanning 0.04 to 39.0, a plain mean spends nearly the whole
budget on six conditions. Normalising each condition -- relative
rather than absolute error, or the existing per-condition
`bayesianRMSEweight` column set to the inverse of each condition's
prior scale -- would let the remaining 35 conditions influence the
result. The column is already parsed and applied faithfully; it is
uniformly 1.0 today, so this needs only data. Note this changes what
"RMSE" means and breaks comparability with MATLAB, so it belongs to
Track B, not Track A.

**B3. Cut the dimensionality.** 4834 parameters against 400 samples
per generation, with every proposal displacing all of them at once,
gives no credit assignment: a good move in one kcat is masked by
thousands of simultaneous bad ones. The identifiability screens
already computed in the shared scratch (`sensitivity.npy`,
`reachable.npy`) can restrict tuning to parameters that can affect
the objective at all.

**B4. Reconsider what is returned.** The best particle fits far better
(1.2038) than the blend (4.1822) but moves ~99% of all kcats, which is
not a defensible posterior. Neither vector is satisfactory; see "Two
kcat vectors" in `bayesian_tuning_plan.md`.

### Recommendation

A1, then A2, then B1. A1 is a day's work at most and decides whether
there is a scoring bug to fix or a genuine 25% method shortfall to
close. A2 is required for correctness regardless. B1 is where the
performance actually is, but it is unmeasurable before A1 and it
deliberately departs from MATLAB, so it should be a named variant
rather than a change to the faithful path.

Also outstanding: the seed-to-seed spread of the port's final RMSE is
unmeasured. The same seed at `n_proc` 16 and 63 reproduces identically,
so within-seed determinism holds, but a single seed cannot say how much
of a 25% gap is noise. Three or four seeds at reduced
`max_generations` would bound it cheaply. And OKP's blend movement
(50.5% here against MATLAB's 15.4%) is the one per-source figure that
does not match, which points at source classification rather than the
sampler.
