# Hand-off: reproduce MATLAB's Bayesian kcat tuning in Python

Goal: run the ported `bayesian_kcat_tuning` as close to GECKO MATLAB's
`bayesianSensitivityTuning` as the port allows, and report how close the
result gets. MATLAB's reference run on this model went **RMSE 9.5544 ->
0.9289 in 31 generations**.

This is a measurement task, not a development task: run it, record the
trace, compare against the numbers below, and report deviations. Change
the algorithm only if you find a genuine divergence from MATLAB, and say
what you changed and why.

## What is already known

Results of the full 31-generation run are in
`matlab_replication_results.md`: **9.6229 -> 1.2038**, with the search
machinery reproducing MATLAB's samples/accepted sequence generation
for generation.

That run retracts what this section used to claim. Scoring the untuned
model gives **9.6229** here, which is 0.7% from the **9.5544** quoted
below but 12% from the **8.6001** in GECKO's own exported run
(`tutorials/full_ecModel/output/rmse_trace.tsv`). Those two MATLAB
figures describe the same deterministic quantity and disagree by 11%,
so the distance function is *not* established as faithful, and a
shortfall in the final RMSE is not safely attributable to the search.
Reconcile the prior RMSE before treating any final number as a target.

A first replication (before two fixes landed) reached **1.4849 in 16
generations**, trace:

```
7.4277 6.9314 5.9830 5.2505 4.9741 4.2706 3.8734 3.1296
3.0258 2.4268 2.3227 2.1786 1.6619 1.5848 1.4849 1.4849
```

MATLAB's trace for comparison: 9.5544 (prior), 4.7157 (gen 2), 3.4267
(4), 2.4551 (6), 1.6337 (8), 1.2842 (10), 0.9730 (14), 0.9289 (22+),
plateauing around generation 14 of 31.

Two MATLAB behaviours were missing from that first run and are now in
the code; a re-run with them is the first thing to redo if no result is
available yet:

* **Biological bounds on proposals** (`tuning.kcat_bounds`): 1e-2 to
  1e4 1/s, or prior/100 to 1e8 when the prior already exceeds 1e4.
  Ported from `proposeSimple`.
* **Evenly spread parents** (`GeckoTransition.rvs_batch`): systematic
  resampling so each accepted particle is perturbed a near-equal number
  of times, matching MATLAB's "select parents with minimal duplication".

## Known remaining differences from MATLAB

Report results in the light of these; do not silently "fix" them.

* `make_anaerobic` / `change_protein_biomass` have no organism-agnostic
  implementation in geckopy, so conditions needing an anaerobic switch
  or a protein-content adjustment run without one. `Ptot` is `NaN`
  throughout this dataset, so the second hook would be inert anyway.
* MATLAB's dead knobs are deliberately absent: `targetAccept`
  (`minKeep` always binds first) and `varianceCap*` (clamps a reported
  diagnostic only, never feeds back).
* MATLAB's low-rank PCA kernel and explore/exploit mixture were
  abandoned upstream; the committed MATLAB code uses the same diagonal
  proposal this port implements.

## Isolation: another run may be in progress

A replication on a 16-core allocation may be running concurrently from
`/cephyr/NOBACKUP/groups/compmeteng/geckopy-bayesian-scratch/`. Do not
write into that directory: its logs (`logs/matlab_replication*.log`) and
saved arrays (`matlab_repl_kcat.npy`, `matlab_repl_blend.npy`,
`sensitivity.npy`, `reachable*.npy`) would be overwritten and both sets
of results lost.

Use your own scratch directory and your own git worktree, e.g.

```bash
mkdir -p /cephyr/NOBACKUP/groups/compmeteng/geckopy-matlab-repl-64
git worktree add ../geckopy-matlab-repl feat/bayesian-kcat-tuning-agent
```

You may *read* the shared directory -- `sensitivity.npy` and
`reachable.npy` are reusable screen outputs -- but write only to your
own. Project storage has a 100k file-count quota shared by the group,
so keep the venv off it (see below) and do not leave large scratch
trees behind.

## Environment

Cluster: C3SE Vera. Nodes are 64-core; an allocation is whatever SLURM
gave you.

```bash
module load Python/3.12.3-GCCcore-13.3.0
python -m venv "$VENV"           # see note on placement below
"$VENV/bin/python" -m pip install -e /path/to/geckopy'[dev,bayesian,tutorial]'
"$VENV/bin/python" -m pip install gurobipy
```

* **Do not put the venv in the repo or in `$HOME`.** Home has a 60k
  *file count* quota that a venv (~34k files) will blow. Put it on
  node-local `/local/tmp.$SLURM_JOB_ID/...` (fast, no quota, wiped when
  the job ends) and rebuild it per job.
* **Gurobi** is the solver for all FBA. `~/gurobi.lic` is a WLS licence
  and `GRB_LICENSE_FILE` already points at it; `tests/conftest.py` sets
  `cobra.Configuration().solver = "gurobi"` for test runs, but a
  standalone script must set it itself:
  `cobra.Configuration().solver = "gurobi"` before loading the model.
  A LAN token-server licence is also available at
  `/apps/Arch/software/Gurobi/12.0.1-GCCcore-13.3.0/gurobi.lic` if the
  WLS endpoint ever stalls.
* Set `Threads = 1` per worker. Gurobi defaults to using every core,
  which oversubscribes badly once workers fill the allocation;
  particle-level parallelism is the useful axis. Set it with a
  `gurobi.env` file in the working directory the run is launched from:

  ```bash
  printf 'Threads 1\n' > "$RUNDIR/gurobi.env"
  ```

  Gurobi reads that file whenever an environment is created, so every
  pool worker picks it up (they inherit the parent's cwd) and logs
  `Set parameter Threads to value 1` at startup -- check for one such
  line per worker. Setting `model.solver.problem.Params.Threads = 1` on
  the parent model does *not* reach the workers: `ProcessPool` pickles
  the model and each worker rebuilds its own Gurobi problem on
  unpickling, at the default `Threads 0` (all cores).

## Cores

**Do not hard-code a core count.** `bayesian_kcat_tuning(n_proc=...)`
defaults to `cobra.Configuration().processes`, which follows the
allocation, so the simplest correct thing is to leave `n_proc` unset.
If you want it explicit:

```python
import os
n_proc = int(os.environ.get("SLURM_CPUS_ON_NODE", os.cpu_count() or 1))
```

Cap it at the per-generation sample count -- workers beyond the batch
size idle. Measured scaling on 16 cores (16 particles, `Threads=1`):
1.97x at 2, 3.86x at 4, 7.19x at 8 (90% efficiency), 11.62x at 16 (73%).
Each worker holds its own `EcModel` copy at ~0.67 GB resident, so budget
memory as `0.67 GB x n_proc`.

Rough cost: one particle is 41 FBA solves, ~0.87 s wall at 16 cores.
MATLAB's schedule is ~10k evaluations to generation 16, which took
178 min on 16 cores.

On a 64-core allocation, expect roughly 3-4x that throughput -- the
scaling above is still climbing at 16, though efficiency falls (73% at
16, so do not expect a clean 4x). Two things to check at that width:

* **Memory**: each worker holds its own `EcModel` copy at ~0.67 GB, so
  64 workers is ~43 GB resident. Fine on a 256 GB node, but confirm
  against what SLURM gave you.
* **Licence checkout**: each worker checks out a Gurobi WLS token over
  the network at pool startup (~1.1 s each, and 8 concurrent checkouts
  measured at ~1.1 s wall total, so they parallelise). 64 at once has
  not been measured here; if pool startup stalls, switch
  `GRB_LICENSE_FILE` to the LAN token server given above.
* `n_proc` still wants capping at the per-generation sample count, but
  MATLAB's schedule starts at 1000 samples so 64 workers are all used.

## Inputs

* **Branch**: `feat/bayesian-kcat-tuning-agent`. Work in your own
  worktree (`git worktree add ../geckopy-matlab-repl
  feat/bayesian-kcat-tuning-agent`) so you do not collide with other
  work on the same checkout.
* **Model**: `tutorials/full_ecModel/models/ecYeastGEM_preTune_GECKOderived.yml`
  -- a copy of GECKO's own
  `tutorials/full_ecModel/models/ecYeastGEM_beforeBaySens.yml`
  (`develop4`), i.e. the pre-Bayesian-tuning ecModel MATLAB feeds to the
  tuner. It is gitignored (`ecYeastGEM*.yml`), so copy it into your
  worktree:

  ```bash
  cp ~/Vera/github/GECKO/tutorials/full_ecModel/models/ecYeastGEM_beforeBaySens.yml \
     tutorials/full_ecModel/models/ecYeastGEM_preTune_GECKOderived.yml
  ```

  8001 rxns, 3893 mets, 4834 tunable kcats, 1144 enzymes, growth
  0.1307 /h. **Do not** use `protocol.py`'s own `ecYeastGEM.yml`: that
  is a Stage 3 output, already past `sensitivity_tuning`, i.e. already
  carrying the hand-fix this method is meant to replace.
* **Data**: `tutorials/full_ecModel/data/bayesian{FluxData,MaxGrowth,ZeroExch}.tsv`,
  committed, copied verbatim from GECKO. 33 flux conditions + 8
  max-growth conditions.

## Hyperparameters

Verbatim from `GECKO/tutorials/full_ecModel/YeastGEMAdapter.m`. Note
`kcatSources = {'OpenKineticsPredictor','brenda','custom'}` for this
model -- **not** the `BayesianParams` default of dlkcat/brenda/custom,
which would leave this model's 1175 OKP kcats plus `standard` and
`isozymes` unclassified and silently on the `*_default` tier.

```python
from geckopy.adapter.params import BayesianParams, SourceGroupRule

PARAMS = BayesianParams(
    sigma0_log_default=0.3,
    source_groups={
        "okp": SourceGroupRule(sources=["OpenKineticsPredictor"], match_okp=True),
        "brenda": SourceGroupRule(sources=["brenda"]),
        "custom": SourceGroupRule(sources=["custom"]),
    },
    sigma0_log_source={"okp": 0.25, "brenda": 0.2, "custom": 0.1},
    shrink_thr_default=5.0,
    shrink_thr_source={"okp": 3.0, "brenda": 10.0, "custom": 12.0},
    force_prior_thr_default=3.5,
    force_prior_thr_source={"okp": 2.5, "brenda": 11.0, "custom": 13.0},
    sparsity_threshold=0.5,
    schedule_generations=[1, 2, 9, 15],
    schedule_samples=[1000, 800, 600, 400],
    min_keep=0.3,
    max_keep=0.6,
    rmse_threshold=0.2,
    max_generations=31,      # MATLAB's run; it plateaus around 14
)
```

## Runner

```python
import logging, os, sys, time
import numpy as np
import cobra
cobra.Configuration().solver = "gurobi"

from geckopy import ModelAdapter
from geckopy.utilities.load_ec_model import load_ec_model
from geckopy.kcat_sensitivity_analysis.bayesian.data import load_bayesian_data
from geckopy.kcat_sensitivity_analysis.bayesian.tuning import bayesian_kcat_tuning

logging.basicConfig(level=logging.INFO, format="%(message)s", stream=sys.stdout)

tutorial = "tutorials/full_ecModel"
adapter = ModelAdapter.from_folder(tutorial)
model = load_ec_model("ecYeastGEM_preTune_GECKOderived.yml", adapter=adapter)
bay = load_bayesian_data(adapter)

t0 = time.perf_counter()
res = bayesian_kcat_tuning(
    model, adapter=adapter, params=PARAMS, bay_data=bay,
    selection="truncation",       # MATLAB-faithful
    n_proc=None,                  # follows the allocation
    seed=0, verbose=True,         # verbose logs best RMSE per generation
)
print(f"{(time.perf_counter()-t0)/60:.1f} min, {res.n_generations} generations")
print("RMSE trace:", [f"{r:.4f}" for r in res.rmse_trace])
```

Run it detached with unbuffered output and a log file; it is hours long.

## What to report

1. **The RMSE trace**, generation by generation, against MATLAB's above.
2. **Wall time**, `n_proc` used, and the node's core count.
3. **Per-source movement** of the returned model, as fractions of kcats
   changed by more than 2% -- comparable to MATLAB's report (custom
   100% unchanged, brenda 97.9%, OKP 84.6%, unlabelled 65.0%). Expect
   the port to change nearly everything; see the "best particle vs
   blend" section of `bayesian_tuning_plan.md` for why that is
   MATLAB-faithful and what is planned about it.
4. **The blend**, if you can: `res.posterior_trace[-1].kcats` is the
   shrunk, snapped point estimate MATLAB reports but never scores.
   Score it with
   `tuning._score_kcat_vector(model, idx, ids, bay, excarbon, bio_rxn, blend, make_anaerobic=None, change_protein_biomass=None)`
   and report its RMSE next to the best particle's. This is an open
   question, not a settled one.

## Pitfalls that have already cost time here

* **Scoring must start from a cold solver basis.** `_reset_solver_basis`
  handles this; do not remove it. Without it the same kcat vector scores
  differently depending on what the worker scored before, results vary
  with `n_proc`, and every generation selects on contaminated distances.
* **`model.solver = "gurobi"` is a no-op in cobra** when the interface
  is already Gurobi, so it cannot be used to force a fresh solver.
* **`pgrep -f <pattern>` matches your own shell** if the pattern text
  appears in its command line. A chained "wait for the previous job"
  loop written that way waits for itself forever. Wait on a PID
  (`while kill -0 $PID`) instead.
* **`tail -n +1 -f` replays a log from the start**, so a monitor can
  resurface an old traceback as if it were new. Use `tail -n 0 -f`.
