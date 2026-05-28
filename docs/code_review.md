# geckopy code review — open items

Remaining items only. These are **modeling-intent / optional design decisions**;
each would change the numbers the model produces, so they are choices, not bugs.

## Default kcat aggregation (max / mean) vs median
`fuzzy_kcat_matching` now takes `aggregate="max"|"median"`, and
`fill_kcats_from_isozymes` takes `aggregate="mean"|"median"|"max"`; the
per-subsystem standard kcat now uses the median (like the global standard).
The **defaults are still `max` / `mean`** for MATLAB-GECKO reproducibility.
Open decision: whether to flip the shipped defaults to `median` (more robust to
assay outliers / engineered mutants), and whether to extend the flag to the
`brenda/parse.py` range collapse (currently always the upper bound) and to
`apply_kcat_list` (which already supports `criteria="median"`).

## get_conc_control_coeffs: single-solve replacement (validated, with solver caveat)
`limit_proteins/get_conc_control_coeffs.py` re-solves the LP per candidate protein
(a 2× finite-difference step). The `prot_<id>` metabolite **shadow price** from a
*single* solve gives the same coefficient — verified across many scenarios in
`tests/test_get_conc_control_coeffs.py` (single + two-enzyme models, `limit`
parameter, non-binding, subset, infeasible). One regime where they diverge: when
the 2× step crosses an LP breakpoint, the finite-difference reports the *averaged*
slope while the shadow price is the (more correct) *local* marginal.

**Solver caveat**: `optlang`'s **scipy** backend does *not* implement LP duals
(`constraint.dual` raises `NotImplementedError`), so the shadow-price approach is
**not** available there. It works on `glpk`, `glpk_exact`, and `gurobi` (and
presumably `cplex`). If switched to default, scipy-only environments would lose
this code path; the safe move is to make it opt-in (e.g. `method="shadow"`) or to
detect scipy and fall back to finite-difference.

(The reaction-level `reduced_cost` is *not* a substitute — it is scaled by the
number of coupled pseudometabolites, so it reads 2× the true coefficient.)

## ecFVA reports an envelope, not exact per-reaction bounds
`utilities/ec_fva.py`. For a reaction split into several isozyme variants
(and/or forward + `_REV`), each variant's extreme can occur in a *different* LP
solution, so reducing per variant and summing gives
`Σ_v max_k v_variant(optᵏ) ≥ max_k Σ_v v_variant(optᵏ)` — an **outer envelope**
wider than the true combined range (exact for single-variant reactions,
conservative otherwise). Now documented in the `ec_fva` docstring. Fix option
(changes results): read each conventional reaction's combined flux from its own
max/min solution (the diagonal) instead of nan-reducing every variant across every
column. Left as a semantics decision.

## Greedy relaxation over-relaxes (kept as-is by request)
`relax_proteomics_greedy` fully unconstrains each picked enzyme, so it yields a
looser model than `flexibilize_enz_concs` (the two are intentionally kept as
complementary methods — `docs/relaxation_methods.md`). Not changed per request;
recorded only so the trade-off stays visible.

## MW of `X` residue
`databases/mw.py` uses an unweighted mean of the 20 standard residues for `X`
(differs from MATLAB's 126.5 Da by ~7 Da). Minor; left for consistency discussion.
(Empty-sequence handling now returns NaN.)
