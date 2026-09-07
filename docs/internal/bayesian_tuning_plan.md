# Bayesian (ABC-SMC) kcat tuning — port plan

Status: resumed on `feat/bayesian-kcat-tuning`. Open decisions below are
now resolved (see each item) after a fuller review of MATLAB's tuning
history and a closer read of pyABC's actual generation-loop mechanics.

## Handoff (read this first)

This work lives on **`feat/bayesian-kcat-tuning-agent`**, forked from
`feat/bayesian-kcat-tuning` at commit `fdf8d2c`. That parent branch no
longer exists locally or on `origin`, and `fdf8d2c` is in this branch's
history, so no reconciliation is outstanding: this branch is the single
line of work.

### Done (Sequencing steps 1–10, all committed, all tests green)

Everything in the "Module layout" section below is implemented in
`src/geckopy/kcat_sensitivity_analysis/bayesian/` **except** `pruning.py`
and `cli.py` (deferred, see Sequencing step 12). `BayesianParams` was
reworked (`source_groups` dict, dropped `target_accept`/`variance_cap_*`,
fixed `sparsity_threshold`'s default) and `template.py` extended to
render dict-valued fields. Run the fast suite to confirm:

```bash
pytest tests/test_bayesian_*.py tests/test_params.py -q
```

(76 passed as of the latest commit on this branch — see `git log`.)
Commits carry the history; read their messages for the reasoning behind
each design decision — in particular the "Spike results" section below
(why `simulate.py` reuses one persistent model per worker instead of
`EcModel.copy()` per particle) and the discovery, from tracing MATLAB's
actual data flow, that its shrink-weight posterior update never feeds
back into sampling or the returned model (documented in `tuning.py`'s
module docstring) — that finding is why both regularization variants
sample from the same raw accepted-particle population and differ only
in particle weighting.

**Parallelization (`n_proc`) is implemented**, per user request (an HPC
run was the target deployment). First cut reused `utilities/ec_fva.py`'s
hand-rolled `multiprocessing.Pool` pattern; asked whether cobrapy itself
already has something for "repeated FBA on slightly-perturbed model
copies, parallelised" (yes -- it's exactly what
`single_gene_deletion`/`single_reaction_deletion`/etc. do), so this was
refactored to use `cobra.util.process_pool.ProcessPool` instead: same
one-persistent-`EcModel`-copy-per-worker-process shape, but it also
handles a real Windows-specific `multiprocessing.Pool(initializer=...)`
performance issue (opencobra/cobrapy#997) for free, and does a proper
`close()+join()` on exit rather than a bare `Pool`'s `terminate()`. Note
`ec_fva.py` itself still uses the older hand-rolled pattern -- that's a
separate, pre-existing, out-of-scope inconsistency, not something this
branch touched. Sampling itself stays single-threaded in the main
process, so a run is bit-for-bit reproducible across `n_proc` values for
the same `seed` — verified directly (`test_parallel_scoring_matches_serial`,
all 4 combinations, re-verified after the `ProcessPool` refactor). One
caveat surfaced by that test and worth knowing before an unattended HPC
run: Python 3.12 emits `DeprecationWarning: This process is
multi-threaded, use of fork() may lead to deadlocks in the child` under
`n_proc>1` (POSIX's default `fork` start method, likely from Gurobi's
internal threads) — it did not manifest as an actual deadlock in any
test run here, and `ec_fva.py` already accepts the same trade-off in
production, but a long unattended multi-generation HPC run is exactly
the scenario where a rare
fork-related hang would actually bite. If that becomes a real problem,
the fix is calling `multiprocessing.set_start_method("spawn", force=True)`
once, early, before any `bayesian_kcat_tuning(..., n_proc>1)` call
(`ProcessPool` itself doesn't expose a `context=` argument to pick this
per-call -- it just uses whatever the process-wide default is). Slower
pool startup, but immune to this class of bug.

### In progress: Sequencing step 11 (real-scale smoke test)

Per the user's direction, the three `bayesian*.tsv` files were copied
**verbatim from the real GECKO MATLAB tutorial data**
(`GECKO/tutorials/full_ecModel/data/bayesian{FluxData,MaxGrowth,ZeroExch}.tsv`,
`develop4` branch, real experimental data — not synthetic) into
`tutorials/full_ecModel/data/` here, and are committed. The real model
(`ecYeastGEM.yml`) is a **gitignored build output**
(`tutorials/full_ecModel/models/ecYeastGEM*.yml`) — it is *not* committed
anywhere; regenerate it by running `tutorials/full_ecModel/protocol.py`,
or copy it from the main geckopy checkout's `tutorials/full_ecModel/models/`
if one has already been built there.

`tests/test_bayesian_tuning_smoke.py` runs all four combinations
against the real model (`@pytest.mark.smoke`, skipped when the model
file is absent — see its docstring). Real scale here is substantial —
**4834 tunable kcats**, **33 flux conditions + 8 max-growth
conditions**, so each particle costs ~41 real FBA solves.

1. Re-run just that one cheapest combination first and note wall-time:
   ```bash
   pytest tests/test_bayesian_tuning_smoke.py -m smoke -k 'truncation-shrinkage' -v -s
   ```
2. If it completes in a reasonable time (low minutes), run all 4:
   ```bash
   pytest tests/test_bayesian_tuning_smoke.py -m smoke -v -s
   ```
   Watch the two `importance_weighting` combinations in particular:
   `importance_weights.compute_importance_weights` and
   `transition.GeckoTransition.component_logpdf`/`_pdf_single` are
   **not vectorised** (plain Python loops over every tunable parameter,
   called once per particle x parent pair) — at ~4800 parameters this is
   a real, not-yet-measured scaling risk, separate from (but analogous
   to) the `EcModel.copy()` finding in "Spike results" below. If a
   combination is prohibitively slow, that is itself the finding to
   record here (with numbers) rather than something to silently work
   around — vectorising those two functions (replace the per-parameter
   Python loop with a single numpy expression over the whole parameter
   vector) would be the fix, but wasn't attempted since it hadn't been
   shown necessary yet.
3. Once all 4 combinations run, compare RMSE trace, wall-time, and
   per-source-group diagnostics (`frac_active`/`frac_near_prior`/
   `mean_deviation` — printed by the smoke test) across combinations;
   record the comparison in this doc and use it to pick a default
   `selection=`/`regularization=` for the eventual MATLAB port-back, per
   the plan's Non-goals section (this decision is explicitly deferred
   until this data exists).
4. Note: `_SMOKE_PARAMS` in the test uses a tiny schedule (8
   samples/generation, 2 generations) deliberately, to make "does it
   run at all" cheap to check. A real validation run should use a
   larger schedule once the above scaling question is answered — but
   that's a slow, expensive, non-per-commit check per the plan's
   Sequencing notes, not something to default to.

### Not started: Sequencing step 12 (deferred/optional)

`pruning.py` (both variants), a `geckopy bayesian-tune` CLI subcommand,
and an opt-in path wiring `GeckoTransition` into a real `pyabc.ABCSMC`
for massively-parallel runs. Not blocking; do this after step 11's
comparison data exists and a default variant is picked.

## Real-scale results (Sequencing step 11)

Run against `ecYeastGEM_preTune_GECKOderived.yml` — GECKO MATLAB's own
`tutorials/full_ecModel/models/ecYeastGEM_beforeBaySens.yml`, the model
MATLAB feeds to `bayesianSensitivityTuning` (8001 rxns, 3893 mets, 4834
tunable kcats, 1144 enzymes, growth 0.13 /h). Note that
`protocol.py`'s own `ecYeastGEM.yml` is **not** a valid input: it is a
Stage 3 output, already past `sensitivity_tuning` (8 kcats retagged
`sensitivityTuning`, 1 `manual`, growth lifted 0.0013 → 0.43 /h), i.e.
already carrying the hand-fix Bayesian tuning is meant to replace.

Smoke schedule throughout: 8 particles x 2 generations, `n_proc=8`,
`_SMOKE_PARAMS` mirroring `GECKO/tutorials/full_ecModel/YeastGEMAdapter.m`
(`kcatSources = {'OpenKineticsPredictor','brenda','custom'}`, sigma0
0.3/0.25/0.2/0.1, shrink 5/3/10/12, force-prior 3.5/2.5/11/13).

| selection | regularization | wall | RMSE trace |
|---|---|---|---|
| truncation | shrinkage | 170.9 s | 8.98 → 8.59 |
| truncation | importance_weighting | 173.7 s | 8.98 → 8.59 |
| quantile_epsilon | shrinkage | 172.2 s | 8.98 → 8.71 |
| quantile_epsilon | importance_weighting | 177.3 s | 8.98 → 8.71 |

Per-group `near_prior` under shrinkage orders exactly by prior tightness
(custom 0.99 > brenda 0.48 > okp 0.35 > unlabelled 0.27), confirming the
sticky-prior machinery transmits at genome scale.

### Why importance weighting was dropped

- **Degenerate at genome scale.** Measured effective sample size
  (`1 / sum(w^2)`) is **1.000** for both selection variants — weights
  `[1.0, 0, 0]` with smallest weight 2.5e-57 (truncation) and 5.2e-91
  (quantile_epsilon). With 4834 parameters the `prior/transition`
  density ratio spans tens of orders of magnitude, so one particle
  takes all the weight and the variant collapses to "keep one
  particle". Its toy-scale advantage (trusted kcats moving less, 22/24
  seeds) does not survive to real scale.
- **Quadratically expensive.** `component_logpdf` is a Python loop over
  every parameter, called once per particle x parent pair: 0.232 s per
  pair at 4834 parameters. Per generation, excluding FBA: 0.3 min at 8
  particles, 39 min at 100, ~10 h at 400, ~64 h at MATLAB's 1000.

Both point the same way, so the regularization axis was removed
entirely; particle weights are now uniform and per-source
regularization is carried by the priors alone.

### Sampling cost bug (fixed)

The generation loop drew `prior.rvs()` / `transition.rvs_single()` once
per *(particle, parameter)* pair rather than once per particle — both
return a full parameter vector, so this was a 4834-fold multiplier:
395.6 s for a single generation-1 particle, ~1.75 h/particle from
generation 2 on. It was also wrong for the transition kernel, which
samples one parent and perturbs it as a unit: per-column draws spliced a
different parent into every coordinate.

### The distance function was not deterministic (fixed)

Scoring the *same* kcat vector three times in one process returns
9.6228, 9.2706, 9.3998. After A→B→A the model is byte-identical (0
reactions differing in bounds, 0 in stoichiometry), so this is not
state leaking between particles: the LP has alternate optima and GLPK
warm-starts from the previous solve's basis, so the reported exchange
fluxes depend on solve history. MATLAB is not exposed to this because
RAVEN's `solveLP` builds the problem fresh per call.

Consequences: every generation selected on contaminated distances, so
this was a correctness problem, not just a reproducibility one; it
affected the serial path too, so `n_proc=1` was not a workaround; and
`test_parallel_scoring_matches_serial` cannot catch it because a
2-parameter toy LP has no meaningful degeneracy. It also explains why
runs of the same combination differed across processes (8.39 / 8.45 /
8.59) and across `n_proc` values (`n_proc=4` → 8.587, `n_proc=8` →
8.450, same seed, same process).

What was ruled out, in order: state leaking between particles (after
A→B→A the model is byte-identical — 0 reactions differing in bounds, 0
in stoichiometry); the solver object (`model.solver = "gurobi"` is a
no-op in cobra when the interface is unchanged, so an apparent
"rebuild per score" test had in fact rebuilt nothing); and Gurobi's
concurrent LP (pinning `Method`/`Threads` changes nothing, because a
solver resuming from an already-optimal basis returns that vertex in
zero iterations regardless of method).

Measured, same vector scored A, B, A, A on the real model:

| configuration | deterministic | s/score |
|---|---|---|
| GLPK, default | no (9.62 / 9.27 / 9.40) | ~60-90 |
| Gurobi, default | no (spread ~1e-3) | 9.3 |
| Gurobi, Method=1 Threads=1 | no (spread ~1e-4) | 6.8 |
| Gurobi, Method=0 Threads=1 | no (spread ~7e-4) | 9.6 |
| **Gurobi + `problem.reset()`, default** | **yes** | **6.4** |
| **Gurobi + `problem.reset()`, Method=1 Threads=1** | **yes** | **7.0** |

Fix: `tuning._reset_solver_basis` discards the incumbent basis before
each particle's solves, so a worker's result cannot depend on which
particle it scored previously. `reset()` alone is sufficient — no
solver parameters are imposed by library code — and it is a no-op for
interfaces whose problem object exposes no `reset`. Guarded by
`test_scoring_is_independent_of_previous_particle`, which has to be a
smoke test since the degeneracy only exists at real scale.

### Parallelism and where the time goes

Measured on the real model (4834 tunable kcats, 33+8 conditions),
Gurobi, `Threads=1` per worker, 16 particles scored per configuration:

| n_proc | total | per particle | speedup | efficiency |
|---|---|---|---|---|
| 1 | 162.4 s | 10.15 s | 1.00x | - |
| 2 | 82.3 s | 5.14 s | 1.97x | 99% |
| 4 | 42.1 s | 2.63 s | 3.86x | 97% |
| 8 | 22.6 s | 1.41 s | 7.19x | 90% |
| 16 | 14.0 s | 0.87 s | 11.62x | 73% |

Each worker holds its own `EcModel` copy at 0.67 GB resident, so memory
is not a constraint at these core counts (16 workers ~ 11 GB).

`Threads=1` is deliberate: Gurobi defaults to `Threads=0` ("use what is
available"), which oversubscribes badly once workers fill the cores.
Particle-level parallelism is the useful axis; per-solve threading only
competes with it.

Proposal sampling runs serially in the main process, so it caps the
achievable speedup. `GeckoTransition.rvs_single` was a per-parameter
scipy loop costing 1313 ms per particle -- 525 s of serial time for a
400-particle generation, against 349 s of scoring on 16 cores, i.e. the
serial part dominated. Vectorised (`lognorm(s=h, scale=p)` is
`p * exp(h * Z)`, one draw over the whole vector) it costs 1.19 ms, most
of which is building the `Parameter` dict; a 1000-particle generation's
proposal cost falls from ~22 min to 1.2 s. Equivalence checked against
the previous implementation (KS p=0.65; per-coordinate log-ratio mean 0,
sd = bandwidth) and guarded by
`test_rvs_single_draws_the_intended_lognormal_around_its_parent`.

With that removed, scoring dominates and cores pay off up to the batch
size: for a 400-particle generation, 16 cores give ~349 s against ~564 s
on 8. Efficiency is 73% at 16 cores versus 90% at 8, so 16 is the right
choice when the cores are already allocated, and 8 is the more
economical point if they are shared.

Remaining serial cost worth knowing: generation 1 draws from
`build_kcat_prior`'s pyabc `Distribution` at 82 ms per particle (82 s
for a 1000-particle first generation). Only generation 1 uses it.

### Identifiability: what the data can actually inform

The point estimate is the deliverable, and a kcat should stay at its
prior unless moving it demonstrably helps -- so the first question is
how many kcats the 41 measured conditions can inform at all.

| screen | kcats | note |
|---|---|---|
| active at the FBA optimum | 499 | misleading: an FBA vertex is maximally sparse |
| FVA reachable, current kcats | 3955 | can carry flux in >= 1 condition |
| FVA reachable, protein pool relaxed | 3950 | ceiling over *any* kcat assignment |
| `\|dRMSE\| > 1e-3` under a x2 perturbation | 1939 | measurably moves the residuals |
| `\|dRMSE\| > 1e-2` | 81 | strong influence |

Reachability must come from FVA, not from one FBA solve: a reaction
idle at the optimal vertex can still carry flux under other kcats, and
the single-solve screen undercounts by 8x. Relaxing the protein pool
changes nothing (3950 vs 3955), so the ~880 unreachable rows are
blocked by stoichiometry and the condition constraints, not by enzyme
cost -- no kcat assignment can make them relevant.

Influence is steeply skewed: nearly everything reachable responds
slightly, but the count falls two orders of magnitude between
`\|dRMSE\|` 1e-3 and 1e-2. This is the "sloppy model" spectrum, and it
matches what MATLAB's own converged run did -- it changed **338 of
4834 kcats (7%)**, which sits between those two thresholds:

| source | total | % unchanged (MATLAB) |
|---|---|---|
| custom | 207 | 100.0% |
| brenda | 3198 | 97.9% |
| OpenKineticsPredictor | 1175 | 84.6% |
| unlabelled | 254 | 65.0% |

Caveat: this is one-at-a-time sensitivity, so it understates parameters
that only matter jointly (co-limiting enzymes in a pathway).

Distance-function validation: scoring the untuned model gives RMSE
**9.6229** against MATLAB's reported prior RMSE **9.5544** -- 0.7%
apart on the same model, so the port's distance is faithful. MATLAB's
trace then runs 9.55 -> 4.72 (gen 2) -> 0.97 (gen 14) -> 0.9289,
plateauing around generation 14 of 31.

### Sampler: tempered SMC instead of hard-epsilon ABC

The ported scheme weights particles by prior/transition density, which
degenerates as dimension grows (ESS 1.000 at 4834 parameters, and ~1
even at 50 in a synthetic test). Subsetting alone does not fix that, so
the sampler was replaced rather than shrunk:

* pseudo-likelihood `exp(-d^2 / 2 tau^2)` instead of hard acceptance,
  so weights vary smoothly;
* `tau` lowered adaptively, each step chosen by bisection to hit an ESS
  target, with a forced minimum decrement so the schedule cannot stall
  when the population sits far from the data;
* systematic resampling plus random-walk Metropolis rejuvenation, so
  particles move rather than merely being reweighted.

Two failures were found on synthetic problems before spending real
compute: the epsilon schedule let the best distance degrade across
generations (no elitism, non-monotone epsilon), and the tempering
stalled at a temperature where the pseudo-likelihood was flat.

Calibration was then checked on synthetic data with known truths (24
replicates, 30 parameters): **90% intervals covered 94.4%** and 50%
intervals 56.9% -- slightly conservative, i.e. honest rather than
overconfident, with intervals 27% narrower than the prior. This is what
makes a null result on the real model interpretable: poor contraction
there is now a statement about the data, not about the sampler.

### First real-scale posterior, and why `tau` must be anchored

On the 81-kcat subset (`\|dRMSE\| > 1e-2`), 6800 evaluations in 89 min
took the best RMSE from 9.62 to 6.88, and every parameter contracted
(100% below 0.9, 65% below 0.5). But contraction was **flat at ~0.45
across a 16x range of influence**, which is not what genuine
per-parameter learning looks like.

The cause is `tau_final`: at 0.5 against a mean distance of ~7.2,
`exp(-d^2/2 tau^2)` is negligible for every particle, so the sampler
concentrates on the best-fitting corner and reports a width set by the
temperature rather than by the data. `tau` represents model discrepancy
plus measurement error and must be anchored to the *attained* fit --
the experimental files carry no error estimates (only
`bayesianRMSEweight`, all 1), so it has to be swept around what each
subset can actually achieve, not asserted below it.

The 81-kcat subset also caps the achievable fit at ~6.9 against
MATLAB's 0.93, confirming that a 1e-2 threshold is too aggressive:
parameters that matter only jointly are excluded by a one-at-a-time
screen.

### Replicating MATLAB's result

MATLAB's reference numbers for this model: 9.5544 -> 0.9289 in 31
generations, plateauing around generation 14. GECKO's own exported run
(`tutorials/full_ecModel/output/rmse_trace.tsv`) instead records
8.6001 -> 0.8715; the two disagree by 11% on the prior, which is a
deterministic quantity. See `matlab_replication_results.md`.

A first replication at MATLAB's own schedule (1000/800/600/400 samples
at generations 1/2/9/15, `min_keep` 0.3, `YeastGEMAdapter.m`'s sigma and
threshold values) reached **1.4849 in 16 generations** -- the right
shape, but roughly 1.7x behind MATLAB at matched generations. Two
MATLAB behaviours were missing from the proposal step, and adding them
moved the early generations ahead of MATLAB's curve:

| generation | before | after | MATLAB |
|---|---|---|---|
| 1 | 7.4277 | 4.9865 | - |
| 2 | 6.9314 | 4.3621 | 4.7157 |
| 3 | 5.9830 | 3.5495 | - |
| 4 | 5.2505 | 3.3345 | 3.4267 |
| 5 | 4.9741 | 2.6274 | - |
| 6 | 4.2706 | 2.0517 | 2.4551 |

The two fixes (`b4c1db3`):

* **Biological bounds on proposals**, 1e-2 to 1e4 1/s (or prior/100 to
  1e8 above 1e4), ported from `proposeSimple`. Without them a share of
  every generation's budget went on kcats no organism could have.
* **Evenly spread parents** via systematic resampling, matching
  MATLAB's "select parents with minimal duplication". Independent draws
  leave some accepted particles unused and others repeated, discarding
  population diversity.

**The early lead does not survive to convergence.** Carried to
MATLAB's full 31 generations on a 64-core allocation, the same run
ends at **1.2038** against MATLAB's 0.9289 and 0.8715
(`matlab_replication_results.md`): ahead through generation 6, level
around generation 8, behind from generation 10 onward. In relative
terms the port cuts RMSE 87.5% from its own prior against MATLAB's
89.9%, so the shortfall is smaller than the absolute numbers suggest,
but it is real and it appears late.

The fixes are still worth having -- they beat the pre-fix run at every
generation, and the pre-fix run's *final* 1.4849 is passed by
generation 12. What they do not do is close the gap.

Two results reframe what remains. The port's samples/accepted sequence
reproduces MATLAB's generation for generation, down to the terminal
171/571 fixed point, so schedule, pooling and selection are faithful
and the remaining difference is not in the sampler. And the distance
function is unverified: with the two MATLAB priors disagreeing by 11%,
the earlier 0.7% agreement was a coincidence of picking one of them.
Localising that per-condition disagreement is the prerequisite for
reading any final RMSE as a target.

### Two kcat vectors: the best particle and the blend

Each generation produces two different parameter vectors, and only one
of them is ever tested.

* **Best particle** -- `kcatTop(:, bestIdx)`, the single accepted
  candidate with the lowest RMSE. Every one of its entries was randomly
  perturbed on the way there, so nearly all differ from the prior; most
  of those differences are noise that happened to ride along in a
  vector that fit well overall. This is what MATLAB writes to
  `ecModel.ec.kcat` (`bayesianSensitivityTuning.m:395-396`) and what
  this port writes too.
* **Blend** -- computed from the whole accepted population: the
  per-parameter mean in log space, shrunk toward the prior by
  `shrinkWeight = min(devFromPrior / shrinkThr, 1)`, forced to the
  prior when the deviation is small, then snapped exactly to the prior
  when the remaining change is below `sparsityThr * sigma0log`.
  Averaging over the population cancels the perturbation noise, so a
  parameter the data does not push ends up back at its prior. This is
  the vector that embodies "leave a kcat alone unless the data moved
  it".

MATLAB stores the blend in `kcatTrace` and never scores it, never
proposes from it and never writes it to the model. Its tuning report
(`writeBayesianTuningReport.m:90`, `kcatFinal = kcatTrace(:, end)`)
therefore describes the *blend*, while the RMSE in the same report
describes the *best particle*. The two headline numbers in that report
belong to different parameter vectors: RMSE 0.9289 (best particle,
measured) alongside 93% unchanged (blend, never scored). This port
reproduces MATLAB's model faithfully -- 4788 of 4834 kcats changed --
because it reproduces the same choice of returned vector.

### Planned: ship the blend, and report what we ship (deferred)

To be explored once the port reproduces MATLAB's fit:

1. **Score the blend.** One evaluation. Either the snapping preserves
   the fit (in which case the better-behaved model has been discarded
   every generation), or it destroys it (in which case the "% unchanged"
   statistic is cosmetic).
2. **Judge on held-out conditions, not training RMSE.** 4834 parameters
   against 41 measurements means the best particle is overfitted almost
   by construction, and the blend is a regularised competitor. On
   held-out conditions the blend may well *win*, in which case there is
   no trade-off to manage. The k-fold harness exists for this.
3. **Pick the operating point from a measured curve.** Rather than a
   fixed rule such as "blend half the approved changes", sweep the
   sparsity threshold (or the MAP objective's `tau`) and choose by a
   stated criterion -- e.g. the sparsest model within 5% of the best
   achievable held-out fit. One interpretable knob instead of MATLAB's
   seven thresholds (`shrinkThr*`, `forcePriorThr*`, `sparsityThreshold`).
4. **Polish after snapping.** Once unevidenced kcats are pinned to
   their priors, re-optimise only the survivors -- the standard
   debiasing step after sparse selection. This should recover most of
   the fit lost to snapping while keeping the model sparse.
5. **Report both numbers, attached to the right vectors.** The shipped
   model's RMSE is the headline; the best particle's RMSE is quoted as
   the unconstrained reference, and the gap between them is the price
   of parsimony. Never report a fit for a vector that is not the one
   written to the model -- that is the flaw found above.

### Solver

All FBA runs on Gurobi (WLS licence via `GRB_LICENSE_FILE`).
`tests/conftest.py` sets `cobra.Configuration().solver = "gurobi"` for
the session, falling back to GLPK when gurobipy is absent, and pool
workers inherit the interface with the pickled model. Beyond
determinism, Gurobi cut the smoke run from ~172 s to ~47 s per
selection variant.

## What it is

MATLAB source:
`GECKO/src/geckomat/kcat_sensitivity_analysis/Bayesian/`
- `bayesianSensitivityTuning.m` (856 LOC) — the main loop
- `abc_max.m` (121 LOC) — distance function (FBA + RMSE)
- `loadBayesianData.m` (19 LOC) — loads the 3 experimental TSVs
- `getrSample.m`, `updateprior.m`, `mse.m`, `addCarbonNum.m` — helpers
- ~1370 LOC total

It is an **ABC-SMC** (Approximate Bayesian Computation, Sequential
Monte Carlo) parameter-estimation routine that searches for kcat
values letting an ecModel reproduce experimental flux / growth data.

## Algorithm

Vanilla ABC-SMC core:
1. Lognormal prior per kcat, with per-source initial sigma.
2. Each generation: sample N proposals -> apply each kcat set to the
   ecModel -> solve FBA -> compute RMSE vs experimental data -> accept
   best ~10% (`targetAccept` percentile, gated by `minKeep`/`maxKeep`).
3. Update posterior mean/sigma in log-space from accepted samples.
4. Adaptive epsilon (quantile of the RMSE distribution).
5. Stop when RMSE < `rmseThreshold` or `maxGenerations` reached.

Distance function (`abc_max.m` -> `rmsecal`):
- RMSE between simulated FBA fluxes/growth and experimental
  measurements across multiple carbon-source conditions.
- Two datasets: max growth rates, and exchange fluxes.
- Zero-flux constraints on reactions assumed inactive.
- Carbon-normalised exchange fluxes (uses `excarbon` per reaction).
- Anaerobic switch + protein-content adjustment per condition.
- Failed/infeasible solves get RMSE = 99 (high penalty).

## GECKO-specific extensions (the load-bearing custom code)

These sit on top of the standard SMC and are what make it work on real
models:
- **Per-source sticky priors**: each kcat source (`brenda`, `dlkcat`,
  `custom`) has its own initial sigma, shrinkage threshold, variance
  cap, and force-to-prior cutoff. BRENDA kcats stay near prior;
  DLKcat ones drift freely.
- **Low-rank PCA proposal kernel** (`proposeLowRankMixture`): learns
  the dominant directions of posterior variability (SVD in log-space)
  and proposes along them.
- **Explore/exploit mixture**: alpha=0.65 PCA-guided exploit kernel +
  0.35 broad isotropic explore kernel (`cExpl` inflation).
- **Sparsity enforcement**: after each update, kcats whose change is
  `< sparsityThr * sigma0` snap back to prior.
- **Variance cap per source**: caps posterior/prior variance ratio to
  prevent random walk.
- **Diagnostics**: shrinkage trace, acceptance rate, epsilon, sparsity
  count, diversity, per-source metrics, low-rank dim — mostly consumed
  by tutorial plots.
- **Optional post-run pruning**: parameters whose change is below a
  threshold revert to prior (`enablePruning`).

Hyperparameters: ~15 in `params.bayesian`
(sigma0logDefault, kcatSources, sigma0logSource, shrinkThr*,
varianceCap*, forcePriorThr*, sparsityThreshold, scheduleGenerations,
scheduleSamples, targetAccept, minKeep, maxKeep, rmseThreshold,
maxGenerations). Fixed internals: alpha, tauResidual, sigmaFloorFrac,
adaptFracEarly, cExpl, rMax.

## Experimental data files (per model adapter, in data/)

- `bayesianFluxData.tsv` — carbon sources, growth rates, env
  conditions, exchange rates
- `bayesianMaxGrowth.tsv` — max growth on various carbon sources
- `bayesianZeroExch.tsv` — exchange reactions assumed zero-flux

## Library choice: pyABC (recommended)

Compared pyABC vs ABCpy vs ELFI:

| Library | Strength | Gap for us |
|---|---|---|
| **pyABC** | custom priors, custom Transition kernels, adaptive epsilon, dynamic-scheduled parallelism (dask/MPI/redis/mp), used in systems biology | no built-in per-source shrinkage or low-rank kernel — subclass `Transition` |
| ABCpy | widest algorithm menu | static-scheduled parallelism only |
| ELFI | BOLFI GP surrogate could cut FBA calls | research-grade; high-D kcat space may not suit BOLFI |

pyABC chosen because (a) its `Transition` abstraction is the right
hook for the GECKO-specific kernel, and (b) dynamic scheduling matters:
each FBA solve is ~10-100 ms and the schedule is up to 500 samples x
50-100 generations, so wall-time dominates. MATLAB already uses
`parfor`; pyABC's dask backend matches that.

Refs:
- pyABC docs: https://pyabc.readthedocs.io/en/latest/
- pyABC JOSS paper (2022): https://arxiv.org/pdf/2203.13043
- ABC-SMC practical guide (2025): https://arxiv.org/html/2511.21587

### What pyABC absorbs vs what we still write

pyABC gives: SMC loop, adaptive epsilon (`MedianEpsilon` ~ GECKO's
quantile rule), parallel sampling, SQLite-backed posterior + per-gen
bookkeeping. ~30-40% of LOC.

We still write (~60%): the FBA+RMSE distance (`abc_max.m` port),
per-source shrinkage / sparsity / force-to-prior / variance-cap
callback, the low-rank PCA `Transition` subclass, the 3-TSV data
loader, and the extra GECKO diagnostics.

Without pyABC we'd write that same 60% plus a 200-400 LOC SMC core
with worse parallelism — so pyABC is worth the dependency.

## geckopy module shape

```
src/geckopy/kcat_sensitivity_analysis/bayesian/
    __init__.py            # re-exports; guards pyabc absence (ImportError on call, not import)
    data.py                # load bayesian{FluxData,MaxGrowth,ZeroExch}.tsv
    simulate.py             # kcat vector -> apply_kcat_constraints -> FBA -> per-condition results
                            # (per-worker persistent model + incremental apply/revert -- see
                            # "Spike results" below; NOT one EcModel.copy() per particle)
    distance.py              # carbon/condition-weighted RMSE (port of abc_max.m's rmsecal half)
    priors.py                # per-kcat lognormal priors by source_group
    selection.py              # selection variants: truncation_select / quantile_epsilon_select
    posterior.py              # shrink/force-to-prior/sparsity-snap blend (trace only)
    transition.py              # GeckoTransition(pyabc.transition.Transition): diagonal fit/rvs/pdf
    diagnostics.py              # per-generation/per-group diagnostics
    tuning.py                  # bayesian_kcat_tuning(model, ..., selection=, regularization=)
    pruning.py                  # optional post-run pruner (deferred, see Q5)
    cli.py                      # geckopy bayesian-tune subcommand (deferred, see Q6)
tests/
    test_bayesian_data.py
    test_bayesian_simulate_distance.py
    test_bayesian_selection.py
    test_bayesian_posterior.py
    test_bayesian_transition.py
    test_bayesian_tuning.py                 # tiny synthetic EcModel, both selection variants
    test_bayesian_tuning_smoke.py           # @pytest.mark.smoke, real-scale
```

## Spike results (Sequencing step 2)

Measured on the real-scale `ecYeastGEM.yml` tutorial model (8001 rxns, 3893
mets, `ec.n_rxns=4842`, `ec.n_enzymes=1144`; Python 3.12.3, `.venv`):

- **`pyabc>=0.13` installs cleanly** and imports its full needed surface
  (`Distribution`/`RV`, `transition.Transition`/`MultivariateNormalTransition`,
  `epsilon.QuantileEpsilon`/`MedianEpsilon`) with no conflicts against
  geckopy's existing pinned deps (cobra 0.31.1, scipy 1.17.1, numpy 2.4.4).
  Installed as `geckopy[bayesian]` per the Dependency note below.
- **`EcModel.copy()` costs ~19.3 s/call; one `model.optimize()` FBA solve
  costs ~86 ms** — a **~225x** ratio (n=20 each, warmed up). This kills the
  "`EcModel.copy()` per particle" assumption `simulate.py` was sketched
  around: a population of even 100 particles/generation would cost >32
  minutes in copying alone, before any FBA time, and would dwarf the actual
  simulation cost by two orders of magnitude.
- **`apply_kcat_constraints` is not safely reversible via `with model:`** —
  it writes stoichiometric coefficients directly via
  `rxn.add_metabolites(..., combine=False)`, which cobra's context-manager
  history does not track/revert (only bound/objective/media changes are).
  So per-particle mutation cannot rely on `with ec_model:` for cleanup
  either; reverting means explicitly re-applying the previous kcat values
  afterward.

**Resulting design for `simulate.py`**: follow `ec_fva.py`'s existing
worker-pool precedent exactly — one `EcModel.copy()` per **worker process**
(via a `Pool(initializer=...)` holding the copy in worker-local global
state), not one per particle. Each particle then: write the new kcat vector
into the worker's persistent `model.ec.kcat`, call
`apply_kcat_constraints(model, update_rxns=affected)` scoped to just the
rows that changed since the last particle solved on that worker, solve,
read results, then leave the model as-is for the next particle to overwrite
(no revert-to-baseline needed between particles — only a revert to the
*previous* particle's state is ever implied, and the next `apply_kcat_constraints`
call already achieves that by writing fresh values over the affected rows).
This is exactly the incremental-update contract `apply_kcat_constraints`'s
docstring already promises ("idempotent... running it twice yields the same
result as running it once"), so no new capability is needed in that
function — only `simulate.py`'s call pattern needs designing around
per-worker reuse instead of per-particle copying.

## Open decisions (resolved)

- **Q1 Library**: pyABC, confirmed. But **not** as the top-level driver:
  `ABCSMC.run()`'s epsilon scheme (`QuantileEpsilon`/`MedianEpsilon`) is
  *prospective and streaming* — epsilon for generation *t* is fixed from
  generation *t-1*'s accepted distances, then the sampler draws proposals
  until a target *accepted* count is reached (`sample_until_n_accepted`,
  unbounded simulation budget). MATLAB's `minKeep` is the opposite: a
  *fixed* batch of `scheduleSamples[gen]` simulations, then a *retrospective*
  top-fraction cut on that same batch. These don't map onto each other
  without fighting the library. Resolution: use pyABC's `Transition`/
  `Distribution` as component libraries inside a small custom generation
  loop geckopy owns, not `ABCSMC.run()`.
- **Q2 Scope**: neither "MVP, no extensions" nor "full port of every
  extension" — informed by a detailed review of MATLAB's actual tuning
  history (not just its code), several "extensions" turned out to be either
  abandoned by MATLAB itself or empirically dead/cosmetic:
  - The low-rank PCA proposal kernel (`proposeLowRankMixture`) was tried and
    *abandoned* in favor of "the SIMPLIFIED VERSION with diagonal proposals"
    that's actually committed to HEAD. Not porting it; pyABC's adaptive
    `LocalTransition`/`MultivariateNormalTransition` is the presumptive
    diagonal-proposal default instead. Also the right shape at genome scale:
    `ec.kcat` has thousands of entries vs. a few hundred particles per
    generation, so a low-rank/full covariance would be rank-deficient.
  - `targetAccept` is confirmed dead (`minKeep`'s floor always binds first);
    `varianceCapDefault`/`varianceCapSource` are confirmed cosmetic (clamp
    only a *reported* diagnostic, never feed back into proposals or point
    estimates). Not porting either — dropping `target_accept`/
    `variance_cap_*` from `BayesianParams` outright rather than shipping
    fields that look load-bearing but silently aren't.
  - `paramTighteningGen`, `shrinkThrDefTight`, `forcePriorThrDefTight` are
    dead leftovers in MATLAB itself, not present in `BayesianParams`
    already — not adding them.
  - What *is* genuinely load-bearing and kept: `minKeep` (the single most
    consequential, empirically-tuned selection lever), `sparsityThreshold`
    (fixing its Python default from 0.3, inside the confirmed-bad
    0.25–0.28 range, to the empirically-working 0.5), per-source
    sticky-prior regularization (mechanism now built as two comparable
    variants — see below).
  - **New**: rather than assuming either "MATLAB's exact mechanism" or "the
    generic ABC-SMC textbook mechanism" is right without evidence, two axes
    are built as swappable variants and compared empirically on a toy
    problem + real-scale smoke test before picking a default:
    - *Selection*: MATLAB-faithful fixed-batch top-`min_keep`-fraction
      truncation vs. pyABC-native adaptive quantile-epsilon streaming
      acceptance.
    - *Regularization*: MATLAB-faithful shrink-weight/force-to-prior blend
      vs. proper SMC importance weighting (`prior_density / transition_density`,
      standard SMC-ABC theory), with sparsity enforcement riding along the
      same axis (MATLAB's post-hoc snap-to-prior vs. a sparsity-inducing
      spike-and-slab prior). **Settled**: importance weighting is degenerate
      at genome scale and quadratic in particle count (see "Why importance
      weighting was dropped"), so this axis is gone — weights are uniform
      and the per-source priors carry regularization. The spike-and-slab
      prior went with it: it existed as that variant's sparsity complement
      and was never wired into the loop. Should the deferred pruning work
      want it, it is in the history.
    The comparison design is recorded in the commits that introduced
    each variant (`selection.py`/`posterior.py`) and in the results
    section above.
- **Q3 Parallel backend**: multiprocessing, not dask — reuse the
  `multiprocessing.Pool` pattern already in `src/geckopy/utilities/ec_fva.py`
  (known-working precedent for parallel LP-heavy ecModel workloads in this
  codebase, no new dependency). Revisit only if profiling shows it's the
  bottleneck.
- **Q4 Return shape**: confirmed — mutate `model.ec.kcat` in place, return
  one `BayesianTuningResult` dataclass, matching the `TunedKcatsResult`/
  `SigmaFitterResult` precedent in `sensitivity_tuning.py`/`sigma_fitter.py`.
- **Q5 Pruning**: deferred, not same-PR — validate the core loop first
  (lower risk). Will eventually have two paths mirroring the regularization
  axis: MATLAB's one-at-a-time re-simulation loop, or a free read-off-the-
  posterior check when the sparsity prior is active.
- **Q6 Tutorial**: unchanged — follow-up PR, after the core module and its
  variant comparison are validated.

Also resolved: `ec.source` doesn't carry a generic `"okp"` tag —
OpenKineticsPredictor kcats are tagged with the raw predictor method name
(e.g. `"CataPro"`). `BayesianParams.kcat_sources` (a flat literal-string
list) can't express "any OKP method" as one trust tier. Resolution: replace
it with a `source_groups` concept — a mapping from trust-tier name to either
explicit literal source strings or a rule matching any configured OKP method
name (read from `OkpParams`) as one group; unmatched sources fall into an
implicit `"unlabelled"` tier (matching MATLAB's `noKcatSource`).

## Dependency note

pyABC pulls in scipy, pandas, cloudpickle (no dask needed — see Q3). Add as
an optional extra: `pip install geckopy[bayesian]` so the core install
stays lean (mirrors the proposed WILDkCAT / wildkcat optional extra).
