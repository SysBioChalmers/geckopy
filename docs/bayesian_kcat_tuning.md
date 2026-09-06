# Tuning kcats against experimental data

A first, brief protocol for using `bayesian_kcat_tuning` to fit an
ecModel's kcats against measured growth rates and/or fluxes. This
covers the shipped, one-call path end to end; it is not a full
reference (see the docstrings in
`geckopy.kcat_sensitivity_analysis.bayesian` for every option) and it
does not build an ecModel from scratch (see
[migrating_from_gecko_matlab.md](migrating_from_gecko_matlab.md) for
that).

## What this does

An ecModel's kcats mostly come from databases (BRENDA), predictors
(DLKcat, OpenKineticsPredictor), or manual curation, and they don't
reproduce measured growth on their own. Bayesian kcat tuning searches
for a kcat vector that does: it repeatedly perturbs kcats, simulates
the model under each experimental condition, and keeps perturbations
that improve the match to measured growth rates and fluxes, weighted
by how much each kcat's source is trusted. It reports which kcats
changed and by how much, not just a fitted model -- the point is
usually the list of corrections, not the number.

## Before you start

- A working ecModel and `ModelAdapter` project (an existing
  `model_adapter.toml` and `models/*.yml`).
- `pip install geckopy[bayesian]` -- the ABC-SMC sampler (`pyabc`) is
  an optional dependency, not part of the base install.
- Experimental data: at least one of a set of measured exchange fluxes
  per condition, or a set of measured maximum growth rates per carbon
  source. Both together works and is preferred -- flux data pins the
  *shape* of metabolism, growth rate pins its *rate*.

## Step 1: Prepare the experimental data

Three tab-separated files go under `<project>/data/`. All three are
optional -- a missing file is skipped, not an error -- but you need at
least one of the first two.

| file | contents | required |
|------|----------|----------|
| `bayesianFluxData.tsv` | measured exchange fluxes, one row per condition | at least one of these two |
| `bayesianMaxGrowth.tsv` | measured maximum growth rate, one row per condition, one active carbon source per row (at `-1000`, i.e. unconstrained) | |
| `bayesianZeroExch.tsv` | reaction IDs assumed to carry zero flux in every condition | no |

Both flux files share one column layout (the same parser geckopy uses
for regular flux data):

```
Condition   Ptot   grRate   glucose (r_1714)   fructose (r_1709)   ...   bayesianRMSEweight   source
glucose     NaN    0.41     -1000              NaN                ...   1                     DLKcat
fructose    NaN    0.338    NaN                -1000              ...   1                     DLKcat
```

`Condition` names the row, `grRate` is the measured growth rate,
`Ptot` is measured protein content (`NaN` if not measured), each
`<met> (<rxn>)` column is a measured or fixed flux for that exchange
reaction, `bayesianRMSEweight` scales how much that row counts toward
the RMSE (usually `1`), and `source` is free text carried through for
your own bookkeeping. `bayesianZeroExch.tsv` is simpler -- a `Rxns`
header, then one reaction ID per line.

## Step 2: Configure `[bayesian]` in `model_adapter.toml`

Every field has a default, so this section can be omitted entirely to
start; the defaults are conservative (every source trusted equally,
no penalty). A configured example, in plain terms:

```toml
[bayesian]
sigma0_log_default = 0.3   # trust for any source not listed below
max_growth_weight = 2.0    # weight the growth-rate data double against flux data
prior_penalty_weight = 0.03  # keep large corrections reproducible across seeds
tie_isozymes = true         # isozyme copies of one reaction move together

[bayesian.source_groups.okp]
sources = ["OpenKineticsPredictor"]
match_okp = true

[bayesian.source_groups.brenda]
sources = ["brenda"]

[bayesian.source_groups.custom]
sources = ["custom"]

[bayesian.sigma0_log_source]
okp = 0.25
brenda = 0.2
custom = 0.1
```

What each knob means:

- **Trust tiers** (`source_groups` + `sigma0_log_source` +
  `sigma0_log_default`). Group your kcats by how much you trust their
  origin, and give each group a standard deviation in log-space:
  smaller means more trusted, and the search needs stronger evidence
  before it moves that kcat far. `custom` (curated by hand) at 0.1 is
  three times more trusted than an unlabelled prediction at 0.3.
- **`max_growth_weight`**. The combined score is
  `(rmse_flux + w * rmse_max_growth) / (w + 1)`. At `2`, growth-rate
  error counts double against flux error -- useful when growth rate is
  measured across many carbon sources but flux is mostly one.
- **`prior_penalty_weight`**. Adds
  `w * mean(((log k - log k0) / sigma0)^2)` to the score: a charge for
  moving a kcat, scaled by how much you trust its prior. Without this
  the fit is often flat along many directions -- several very
  different kcat vectors fit equally well -- and the search returns an
  arbitrary point along that flat direction, which does not reproduce
  if you rerun it with a different seed. The default (`0.03`) is
  calibrated so that large corrections agree in direction and land
  within a handful of fold of each other across seeds, at a small cost
  in raw fit. Set it to `0` to score on fit alone (this is also what a
  run reproducing GECKO MATLAB, which has no such term, must do).
- **`tie_isozymes`**. Several ecModel reactions are one enzymatic
  reaction split across isozyme copies. When copies share a prior
  value and a source, no experimental condition can tell them apart,
  so an untied search is free to invent a difference between them that
  means nothing biologically. Tying gives such copies one shared kcat.
- **`schedule_generations` / `schedule_samples` / `min_keep` /
  `max_keep` / `rmse_threshold` / `max_generations`**. Control the
  ABC-SMC sampler itself: how many candidate kcat vectors are drawn
  per generation, what fraction survive each round, and when to stop.
  The defaults are a reasonable starting point; see the `BayesianParams`
  docstring if a run stalls or converges too slowly.

## Step 3: Decide which kcats are tunable

By default every kcat with a value greater than zero is tunable. On a
genome-scale model that is usually thousands of parameters, most of
which no condition in your dataset can actually inform -- carrying no
flux under any of them, or moving the score by an unmeasurable amount.
Restricting to a smaller, well-chosen set both keeps the search
honest (a change is only reported where the data could see it) and
searches better (CMA-ES/ABC-SMC's step size grows with the square
root of the dimension).

There is no shipped, automatic policy for building this set yet -- it
is a plain boolean array you construct and pass as `tunable_mask`.
`geckopy.kcat_sensitivity_analysis.bayesian.parsimony.identifiability_mask`
gives you the building block: keep a kcat when its one-at-a-time
effect on the score, weighted by trust, clears a threshold you choose.
Building and calibrating that threshold is model-specific work; start
without a mask, and reach for one once you can see which kcats the
unrestricted run leaves untouched or barely touches.

## Step 4: Run it

```python
from geckopy import ModelAdapter, load_ec_model, save_ec_model
from geckopy.kcat_sensitivity_analysis.bayesian import bayesian_kcat_tuning

adapter = ModelAdapter.from_folder("path/to/project")
model = load_ec_model(adapter=adapter)  # models/ecModel.yml by default

result = bayesian_kcat_tuning(model, adapter=adapter, n_proc=8, seed=0)

save_ec_model(model, "ecModel_tuned.yml", adapter=adapter)
```

`model` is mutated in place: on return, every tunable kcat carries the
final generation's best-scoring value, and the model's kcat
constraints are already applied -- nothing further is needed before
simulating or saving it. `n_proc` parallelises scoring each
generation's candidates across processes; omit it to use
`cobra.Configuration().processes`, or pass `1` to run serially.

If growth needs special handling for your organism -- forcing
anaerobic conditions, or scaling the protein pool with a
condition-specific biomass composition -- pass `make_anaerobic` and/or
`change_protein_biomass` callables; geckopy has no organism-agnostic
default for either. See `tutorials/full_ecModel/code/anaerobic.py` in
this repository for a worked example.

## Step 5: Read the result

`result` is a `BayesianTuningResult`:

- `rmse_trace` / `objective_trace` -- best plain-fit RMSE and best
  optimised-objective value, per generation (identical unless
  `prior_penalty_weight` is nonzero). A flat trace for several
  generations before the run ends means it stalled, not converged.
- `new_kcat` / `old_kcat` / `rxns` / `groups` -- the tuned and prior
  kcats, parallel arrays, with each kcat's trust-tier source group.
- `converged` -- whether `rmse_threshold` was reached, or
  `max_generations` was hit first.
- `diagnostics_trace` -- per-generation, per-source-group summary
  statistics, for watching how each trust tier moves over the run.

Before trusting a tuned kcat, check more than the final RMSE:

1. **How many kcats changed**, and by how much -- a result that moves
   nearly everything by a little has not identified anything; a result
   that changes a handful by a lot, concentrated in less-trusted
   sources, is the useful shape.
   `geckopy.kcat_sensitivity_analysis.bayesian.parsimony` has the
   tools for this: `n_changed`, `fold_change`, `source_movement`.
2. **Impact share** -- `parsimony.impact_share` reports what fraction
   of the total achievable improvement the changed kcats actually
   carry, distinguishing "few, large, and consequential" from mere
   sparsity.
3. **Reproducibility.** Rerun with a different `seed`. Large
   corrections (say, beyond two-fold) should land in the same
   direction and within a small factor of each other; if they don't,
   the search has found a flat direction rather than a real
   correction, and `prior_penalty_weight` is the first knob to reach
   for.
4. **Whether it generalises**, if you have data you didn't fit on --
   score the tuned model against it and compare to the untuned model.
   A large improvement on fitted conditions with no improvement
   elsewhere is a sign of overfitting to a thin dataset, not a fixed
   model.

## A note on the search itself

`bayesian_kcat_tuning` runs ABC-SMC (sequential Monte Carlo): it
samples and truncates, rather than optimising directly. Internal
validation on ecYeastGEM found that a separable CMA-ES search over the
same screened parameters, same objective, same penalty, reaches a
distinctly better fit (by several standard errors) with far tighter
cross-seed reproducibility -- see
`docs/internal/bayesian_tuning_handover.md`, "Three methods, one
conclusion". That CMA-ES path is not yet wrapped as a single function
in the package; today it means composing `priors`, `distance`,
`simulate` and the tying/screening helpers directly, the way the
project's own analysis scripts do. ABC-SMC via `bayesian_kcat_tuning`
is the fully-supported entry point and a reasonable place to start;
promoting the CMA-ES path to a shipped function is a natural next
step, not yet done.

## Where to go from here

- `docs/internal/bayesian_tuning_handover.md` -- the fuller set of
  findings behind the defaults above (why 0.03, why tying, what
  doesn't work).
- `docs/internal/matlab_replication_results.md` -- the underlying
  experiments, including the still-open question of how to build a
  `tunable_mask` that generalises across models (Open items #6).
