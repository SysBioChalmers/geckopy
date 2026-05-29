# `flexibilize_enz_concs` vs `relax_proteomics_greedy`

Both functions answer the same question — *the measured proteomics constraints
make the ecModel grow too slowly (or not at all); which enzyme bounds do we
loosen, and by how much, to restore feasible growth?* — but they are **not**
duplicates. They use different selection signals, different relaxation step
sizes, and produce results of very different stringency. This note records the
differences and the recommendation (keep both; do not merge).

## Side by side

| aspect | `flexibilize_enz_concs` | `relax_proteomics_greedy` |
|---|---|---|
| Origin | GECKO MATLAB `flexibilizeEnzConcs.m` (faithful port) | legacy geckopy, Carrasco et al. 2023 |
| Which enzyme to relax | largest **control coefficient** (finite-difference: how much each measured enzyme's bound limits growth, from `get_conc_control_coeffs`) | largest **\|shadow price\|** (LP dual of the usage upper bound) |
| Solves per selection | re-solves the LP **once per candidate protein** (expensive) | a **single** `slim_optimize` exposes all duals (cheap) |
| Relaxation step | **gradual**: usage `ub ← conc·(1 + fold_change·k)`; the same enzyme can be picked repeatedly and loosened a bit more each time | **all-or-nothing**: usage `ub ← default_upper_bound` (≈ unconstrained) in one shot |
| Tighten-back | yes — a refinement pass minimises the protein pool at `exp_growth` and **drops** enzymes whose true usage is below their original conc (un-relaxes over-relaxed ones) | none — each picked enzyme stays fully unconstrained |
| Result stringency | **tight** (minimal, proteomics-faithful relaxation) | **loose** (tends to over-relax) |
| Pool fallback | yes — relaxes `prot_pool_exchange` when no single enzyme helps | no |
| Scope | measured enzymes only (`~isnan(concs)`) | any currently-constrained enzyme (`ub < default`) |
| Non-convergence | returns a partial `FlexEnzResult` + warning (global iteration cap) | **raises** `RuntimeError` |
| Result type | `FlexEnzResult` (relaxed enzymes, per-enzyme `frequence`, ratios) | `GreedyRelaxResult` (relaxed→orig-conc dict, step trace, `converged`) |

## Interpretation

- **`flexibilize_enz_concs`** is the right default when you want to stay as
  faithful to the proteomics data as possible: it loosens enzymes a little at a
  time and then tightens back everything that didn't actually need loosening,
  so the final model is the *minimal* deviation from the measurements. The cost
  is speed — the control-coefficient computation re-solves the LP for every
  measured protein on every outer iteration.

- **`relax_proteomics_greedy`** is the right choice when you want a fast answer
  and the infeasibility is dominated by one or two enzymes: a single solve
  ranks all enzymes by shadow price, and fully unconstraining the worst offender
  usually moves growth a lot. The cost is that it **over-relaxes** (no
  tighten-back), so the resulting model is less constrained than necessary, and
  it gives up (raises) rather than degrading gracefully.

## Recommendation

Keep both — they are complementary, not redundant. They could share small
helpers (eligible-enzyme selection, result shaping), but their cores differ
enough that merging into one function would obscure the two algorithms. Two
follow-ups would make them more interchangeable:

1. Give `relax_proteomics_greedy` an optional tighten-back / graded step so its
   output is not systematically looser than `flexibilize_enz_concs` (this is
   the over-relaxation point in the modeling-intent review, §3.3).
2. Align the failure convention — currently one returns a `converged=False`
   result and the other raises; pick one so callers don't special-case each.
