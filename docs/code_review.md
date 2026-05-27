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

## get_conc_control_coeffs: single-solve replacement (validated, not applied)
`limit_proteins/get_conc_control_coeffs.py` re-solves the LP per candidate protein
(a 2× finite-difference step). Tested alternative
(`tests/test_get_conc_control_coeffs.py::test_shadow_price_matches_finite_difference_coeff`):
the `prot_<id>` **metabolite shadow price** from a *single* solve equals the
finite-difference coefficient exactly. (The usage reaction's `reduced_cost` does
**not** — it is scaled by the number of coupled pseudometabolites, so the original
"reduced cost" phrasing was wrong.) Open: swap the per-protein loop for one solve +
shadow prices (≈ n× fewer LPs); this also speeds up `flexibilize_enz_concs`, which
calls it every iteration.

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
