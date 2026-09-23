# What changed in the kcat-assignment pipeline since the June model, and whether it affects tuning

The ecYeastGEM model tuned in June (`ecYeastGEM_preTune_GECKOderived.yml`) was
built by a snapshot of geckopy that has since moved on. This investigates
every difference between that snapshot and mainline today, one change at a
time, and checks whether CMA-ES tuning (top-112 free parameters, 3 seeds,
budget 6800 evaluations, `max_growth_weight=2`) still lands at the same
place. The June model tunes to a distance of 0.798 ± 0.004 (mean ± sd over
3 seeds) — every number below is judged against that noise floor.

Everything here is reproduced by re-running the assignment stage (fuzzy
BRENDA matching through `apply_kcat_constraints`) from raw inputs, isolated
one change at a time, on the June ecModel's own reaction set. A byte-exact
replay of every June bug first reproduced the June kcats 100% identically,
so each change below starts from a verified-identical baseline.

## Summary

Six pipeline changes tune to the same place as June, within noise:

| Change | Tuned distance | vs June |
|---|---|---|
| source labels only | 0.798 | +0.0% |
| substring-wildcard fix | 0.796 | −0.1% |
| tie-break fix | 0.792 | −0.7% |
| min-coefficient fix | 0.791 | −0.9% |
| standard kcat (subsystem, not model-wide) | 0.797 | −0.1% |
| OKP predictions in real tiers (not one bloc) | 0.796 | −0.2% |
| all of the above together | 0.794 | −0.5% |

One change is a clear improvement:

| Change | Tuned distance | vs June |
|---|---|---|
| DLKcat instead of CatPred | 0.671 | **−15.8%** |
| EITLEM instead of CatPred | 0.769 | −3.6% |

One change breaks tuning outright, traced to a specific data-quality issue,
with two candidate fixes:

| Change | Tuned distance | vs June |
|---|---|---|
| BRENDA 2026.1 (unmitigated) | 4.19 | broken |
| BRENDA 2026.1, `top_origin_limit=4` | 0.805–0.850 | same, or +6.5% for the unmitigated combination |
| BRENDA 2026.1, floor filter | 0.731–0.781 | **better**, except paired with DLKcat |
| BRENDA 2026.1, floor + `top_origin_limit=4` | 0.701–0.764 | same or better, every predictor |

## The BRENDA 2026.1 data trap

BRENDA's kcat table has one row for EC 2.3.1.297 at a specific activity of
6×10⁻⁷ µmol/min/mg for *Mus musculus*, joined against a molecular weight of
42,000 g/mol — almost certainly a unit or transcription error in the
source data, giving a derived kcat of 4.2×10⁻⁷ 1/s. BRENDA also has a
*Saccharomyces cerevisiae* specific-activity row for the same EC (20.68
µmol/min/mg, sensible), but no yeast molecular weight for that EC, so the
loader's SA×MW join drops it — the mouse row survives on its own and wins.

At the mainline default `top_origin_limit=6`, any exact-EC BRENDA match,
including specific-activity-derived ones, outranks a prediction. The
4.2×10⁻⁷ kcat gets assigned, the resulting enzyme-usage coefficient
(MW/(kcat·3600)) is enormous, and CMA-ES cannot recover a sane model in
budget: distance 4.19–4.38 depending on aggregation, vs 0.80 for June.

This is not unique to this one EC. Across the whole BRENDA 2026.1 dump,
5,363 of 17,829 specific-activity rows (30%) have no same-organism
molecular weight and are silently dropped; 608 kcat-table rows and 1,913
SA-derived rows sit below 10⁻² 1/s, some of which are likely errors of
the same kind.

### Two candidate fixes

**`top_origin_limit=4`** (existing parameter, no code change): demotes
specific-activity-derived matches (origins 5–6) below predictions, so a
prediction is used instead of the mouse row. Fixes the break
(`top_origin_limit=6`: 4.19 → `=4`: 0.81, close to June) but does cost a
little on the unmitigated pipeline (`b26_o4`: +6.5%) and is redundant with
the fix below once that's applied.

**A relative floor on derived kcats** (`brenda_fix.py`, sandbox-only, not
yet in mainline): drop a specific-activity-derived kcat if it's more than
4 orders of magnitude below the median of the same EC family (widening to
the EC's 3-, 2-, then 1-digit prefix if fewer than 5 reference values are
available). This is the more direct fix — it removes the bad row instead
of demoting a whole tier — and gives the best results seen in this
investigation when combined with `top_origin_limit=4`:

| Predictor | floor alone (limit=6) | floor + limit=4 |
|---|---|---|
| CatPred | 0.731 (−8.3%) | — |
| EITLEM | 0.718 (−10.0%) | 0.701–0.826 |
| unmitigated pipeline | 0.781 (−2.1%) | — |
| DLKcat | **1.80, broken** | **0.72–0.76, fixed** |

An earlier version of the floor also carried a molecular-weight fallback
(use the EC family's median MW when the organism has no MW row, recovering
the yeast SA value directly instead of relying on the floor to catch the
alternative). Tested and rejected: combining it with the floor was worse
than either fix alone, and broke the DLKcat pairing (1.13, and 10/112 top
parameters overlapping with June, vs 96+/112 for clean variants) the same
way the unmitigated data does. The floor's own reference pool is built
from the union of the kcat and SA-derived tables; recovering more rows
into that pool before filtering can quietly lower a family's "typical"
value and let a bad entry back in. Removed; only the floor remains.

### Why the floor alone broke DLKcat, and why `top_origin_limit=4` doesn't

The floor removes the single obviously-bad row (EC 2.3.1.297), but not
every specific-activity-derived kcat that passes the floor is a *good*
value — some are real BRENDA measurements that are simply less reliable
than a direct kcat. At `top_origin_limit=6`, 365 reactions swap from
whatever the predictor said to one of these floor-surviving
specific-activity values — the same 365 reactions, with the same new
values, regardless of which predictor (CatPred or DLKcat) is paired with
BRENDA, since the swap only depends on BRENDA's own matches.

Tuned against CatPred, this is fine (`c3_fl`: 0.731). Tuned against
DLKcat, it isn't (`c5_fl`: 1.80). Tracing this by comparing per-condition
RMSE between the two: the failure concentrates in 7 of 33 flux conditions,
all moderate-to-high growth rate on glucose — exactly where *S. cerevisiae*
should switch on the Crabtree effect and secrete ethanol even aerobically,
because respiratory capacity becomes limiting. The DLKcat-tuned model
systematically **under-produces ethanol and over-produces growth** for the
same glucose uptake (e.g. one condition: measured ethanol 9.5 mmol/gDW/h,
simulated 0.0; growth 0.362 vs the CatPred-tuned model's 0.338, closer to
measured). No enzyme is anywhere near its usage cap in these conditions
(all under 2%), so this isn't a protein-budget bottleneck — it's a missed
metabolic switch. DLKcat's own predictions for the reactions it keeps (the
~700+ reactions not part of the 365-reaction swap) apparently combine with
the floor-admitted BRENDA values to make "too respiratory, not enough
overflow" the solution CMA-ES converges to within budget; CatPred's
landscape doesn't have that trap. `top_origin_limit=4` avoids it by
additionally demoting the floor-surviving specific-activity tier below
*any* prediction, so DLKcat's own values are used more broadly instead.

The exact reaction(s) responsible for the respiration/overflow trade-off
were not pinned down further (would mean comparing the fermentation
pathway — pyruvate decarboxylase, alcohol dehydrogenase — and the
respiratory chain between the working and broken variants for the largest
kcat swaps). The phenomenon and its trigger (floor + DLKcat, without the
origin-limit demotion) are solid; the single-reaction mechanism isn't
required to know **floor + `top_origin_limit=4` together is safe for every
predictor tested**, which is the practical recommendation below.

For CatPred specifically, combining the two mitigations makes almost no
difference versus `top_origin_limit=4` alone (`c3_fl4` vs `c3_o4`: 34 of
4,834 kcats differ; `b26_fl4` vs `b26_o4`: byte-identical). This isn't a
bug — `top_origin_limit=4` already demotes every origin-5/6 BRENDA match
below *any* prediction, and CatPred's coverage is broad enough that almost
every reaction with an origin-5/6 match also has a CatPred prediction
available to fall back to. The floor only has room to matter where a
predictor's coverage is patchier (DLKcat, EITLEM) and a reaction's only
alternative to the bad BRENDA value is standard kcat, not a prediction.

## Recommendation

- **Commit**: nothing yet. The floor fix (`brenda_fix.py`) is sandbox-only
  and needs porting into `geckopy.databases.brenda_loader` before it's a
  real candidate; it should ship together with `top_origin_limit=4`
  becoming the default, not alone, given the DLKcat interaction above.
- **Do not commit**: the molecular-weight fallback. Tested and rejected.
- **Document only**: BRENDA 2026.1's specific-activity-derived tail has a
  genuine outlier problem inherited from the source data, independent of
  this pipeline. EITLEM (already-supported as an OKP method, just a
  different `method=` argument) gives a modest, consistent improvement
  over CatPred and is worth considering as an alternative default
  predictor; DLKcat's larger improvement is real but its parameter
  selection overlaps least with June's (67–83 of 112 top parameters
  in common, lowest of any variant tested), meaning the two predictors
  are tuning meaningfully different aspects of the model, not just
  reproducing the same fit with different labels.

## Method

Reassignment driven by `deltas/pipeline.py` (`Cfg` dataclass toggles each
change), using a private copy of the matcher (`deltas/fuzzy_compat/`) for
the three matcher-bug switches and the real mainline `merge_kcats`,
`assign_standard_kcat`, `apply_kcat_constraints` for everything else.
Screening and CMA-ES: `deltas/variant_tune.py`, same procedure as the
June-reproduction runs (screen every free parameter ×2/÷2 with tied
isozyme copies moved as a unit, top-112 by leverage×sigma0, separable
CMA-ES, 3 seeds, budget 6800, `max_growth_weight=2`). Per-condition
decomposition for the DLKcat trace: `deltas/decompose.py`,
`deltas/enzyme_trace2.py`.
