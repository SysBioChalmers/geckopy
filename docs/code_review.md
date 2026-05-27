# geckopy code review — open items

Remaining items only. These are **modeling-intent design decisions** for the
original GECKO author rather than mechanical bugs — each changes the numbers the
model produces, so they are choices, not fixes. Tagged **(e)** (metabolic-modeling
intent) from the original review's category scheme.

## Max-kcat selection bias
`fuzzy_kcat_matching/_brenda_query.py` (max within a BRENDA level),
`brenda/parse.py` (range collapses to the upper bound), `merge_kcats` (keeps all
rows in the winning tier), and `apply_kcat_list(criteria="max")` together compose
into "max of max of upper-bound". This systematically favours the single highest
reported turnover (assay outliers, engineered mutants), under-estimating enzyme
demand. Faithful to MATLAB GECKO; consider exposing median/percentile as the
default for database aggregation.

## Arithmetic mean of kcats across isozymes / subsystems
`ec_model/pipeline/fill_kcats.py`, `gather_kcats/get_standard_kcat.py`. kcat spans
orders of magnitude; the geometric mean (or median, already used for the global
standard) is the defensible aggregate for log-distributed rate constants.

## Greedy relaxation over-relaxes
`relax_proteomics_greedy` fully unconstrains each picked enzyme, whereas
`flexibilize_enz_concs` ramps gradually and tightens back, so greedy yields a much
looser model (the two are kept as complementary methods —
`docs/relaxation_methods.md`). Optional improvement: give greedy a tighten-back /
graded step. It also ranks by `abs(shadow_price)` (a one-sided test would be more
correct) and reads duals after `slim_optimize()` (prefer `model.optimize()` + the
returned solution).

## Finite-difference control coefficients
`limit_proteins/get_conc_control_coeffs.py` re-solves the LP per candidate protein
with a 2× bound step. The reduced cost / shadow price of the usage upper bound from
a *single* solve is the analytic margin — both faster and not direction-biased.

## ecFVA reports an envelope, not per-reaction bounds
`utilities/ec_fva.py` nan-reduces each reaction's flux across *all* groups'
optimizations, so the reported min/max is the union envelope rather than
independent per-reaction FVA. Likely intended GECKO semantics, but should be
documented; the diagonal gives true per-reaction FVA.

## sigma fit is a 100-point grid over a monotone response
`kcat_sensitivity_analysis/sigma_fitter.py`. Growth is monotone in sigma, so
bisection finds the best sigma in ~7 solves instead of 100, each currently a full
`set_prot_pool_size` + cold solve.

## MW from sequence: `X` residue and empty-sequence handling
`databases/mw.py` uses an unweighted mean for `X` (differs from MATLAB's 126.5 Da,
~7 Da/X); an empty sequence returns 18 Da (water), which a caller could mistake for
a real protein. Consider NaN for empty input.
