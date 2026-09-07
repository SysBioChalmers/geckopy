# Parallel FBA: use `cobra.util.process_pool.ProcessPool`, not raw `multiprocessing`

## Context

While parallelizing `kcat_sensitivity_analysis/bayesian/tuning.py` (repeated
FBA on kcat-perturbed copies of a model, one persistent `EcModel` per worker
process), it turned out cobrapy already solves exactly this shape of problem:
`cobra.util.process_pool.ProcessPool` is what `single_gene_deletion`/
`single_reaction_deletion`/`double_gene_deletion`/`double_reaction_deletion`
(and cobrapy's own `flux_variability_analysis`) use internally for "repeated
FBA on slightly-perturbed model copies, in parallel."

`tuning.py` was first written against a hand-rolled
`multiprocessing.get_context("fork"/"spawn")` + manual `pickle.dumps(model)`
pattern — copied from `utilities/ec_fva.py`, which predates this — then
refactored to `ProcessPool` once this was noticed. `ProcessPool` gets two
things right that the hand-rolled version doesn't:

1. **A real Windows performance fix.** Raw `multiprocessing.Pool(initializer=...,
   initargs=...)` has a slow initarg handoff path on Windows
   (opencobra/cobrapy#997); `ProcessPool` works around it by writing
   initializer + args to a temp pickle file instead, transparently.
2. **Correct pool teardown.** A bare `Pool` used as `with ctx.Pool(...) as pool:`
   calls `terminate()` on exit (abrupt). `ProcessPool.__exit__` does
   `close()` + `join()` first, then cleans up — a graceful shutdown.

## What to do

**`ec_fva.py` should be refactored to match** (see
`kcat_sensitivity_analysis/bayesian/tuning.py`, commit `8a2b960` on
`feat/bayesian-kcat-tuning-agent`, for a worked example of the same change).
Concretely, in `ec_fva.py`'s parallel branch:

```python
# Before (current):
ctx_name = "fork" if sys.platform != "win32" else "spawn"
ctx = mp.get_context(ctx_name)
pickled = pickle.dumps(ec_model)
with ctx.Pool(n_proc, initializer=_init_worker, initargs=(pickled,)) as pool:
    ...

# After:
from cobra.util import ProcessPool
with ProcessPool(n_proc, initializer=_init_worker, initargs=(ec_model,)) as pool:
    ...
```

Drop the `sys`/`multiprocessing`/`pickle` imports if nothing else in the file
needs them; `_init_worker` no longer needs to `pickle.loads(...)` its first
arg (`ProcessPool` handles serialization).

**Then look for other instances of the same pattern**, in geckopy and in
raven-toolbox (`~/github/raven-toolbox` — a sibling repo, not this one).
This handover's search was not exhaustive: `ec_fva.py` is the only
confirmed hand-rolled `multiprocessing`/parallel-FBA instance found in
geckopy's `src/`; a comment in raven-toolbox's
`raven_toolbox/localization/certify.py` referencing FVA multiprocessing
turned out to just be a caller passing `processes=1` to cobrapy's own
(already-`ProcessPool`-backed) `flux_variability_analysis` — not something
needing a fix. Don't take that as a clean bill of health for the rest of
raven-toolbox, though; it wasn't audited beyond that one hit.

**How to spot a candidate**: grep both repos for `multiprocessing`, `mp.Pool`,
`mp.get_context`, or `Pool(` — any hit that pickles a whole cobra/EcModel
object into `initargs` by hand is a candidate.

**How to verify a refactor didn't change behaviour**: assert `n_proc=1` and
`n_proc=2`+ give identical results for the same input, ideally including the
same RNG seed if the function involves any sampling. Precedent:
`tests/test_ec_fva.py::test_parallel_matches_serial` and
`tests/test_bayesian_tuning.py::test_parallel_scoring_matches_serial`.
