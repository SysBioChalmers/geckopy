# geckopy code review

Scope: full sweep of the `geckopy` package (~16 k lines across `ec_model`, `gather_kcats`,
`databases`, `limit_proteins`, `kcat_sensitivity_analysis`, `utilities`, `adapter`,
`get_enzyme_data`, `io`, `cli`). Findings are tagged by the requested categories:

- **(a)** potential bugs
- **(b)** inefficiencies
- **(c)** inconvenient functionality / awkward API
- **(d)** illogical organization
- **(e)** critique of the metabolic-modeling intent
- **(f)** other problems

Resolved items are removed as they are fixed; see `git log` for the corresponding commits.
The entire high-priority tier (former §1.1–1.8) and **all of §2** are addressed. Most of §4
and §5 are done too; what remains below is the modeling-intent tier (§3, design calls for the
original author), a couple of lower-priority guards (§4), and a couple of hygiene items (§5).

Verified **false positives** (no change): §2.3 (currency-pair stripping mutates `s_dense` live,
so cumulative state is tracked), §2.4 (`set_kcat` base-name does not strip `_REV`, so
forward/reverse stay distinct), §2.7-objective (usage fluxes already carry MW scaling, so the
unit-weighted sum minimises mass), §2.11-pubchem-500 (PubChem returns 500 for unparseable names
— caching it is intentional and tested). `pax_db_loader` non-numeric-abundance drop was judged
correct as-is (dropped by request). Remaining IDs are kept stable.

---

## 3. Modeling-intent critiques — **(e)**

### 3.1 Pervasive "max kcat" bias
`fuzzy_kcat_matching/_brenda_query.py:154-158` (max within a BRENDA level), `brenda/parse.py:132`
(range collapses to upper bound), `merge_kcats` (keeps all rows in the winning tier), and
`apply_kcat_list(criteria="max")` default compose into "max of max of upper-bound". This
systematically favours the single highest reported turnover (assay outliers, engineered
mutants), under-estimating enzyme demand. Faithful to MATLAB GECKO, but worth exposing
median/percentile as the default for database aggregation. This is a design call for you as the
original author.

### 3.2 Arithmetic mean of kcats across isozymes/subsystems
`ec_model/pipeline/fill_kcats.py:101`, `gather_kcats/get_standard_kcat.py:205-212`. kcat spans
orders of magnitude; the geometric mean (or median, already used for the global standard) is the
defensible aggregate for log-distributed rate constants.

### 3.3 Greedy relaxation over-relaxes
`limit_proteins/relax_proteomics_greedy.py:122` jumps a limiting enzyme straight to ub=1000
(effectively unconstrained), whereas `flexibilize_enz_concs` ramps gradually and tightens back.
The difference is now documented (`docs/relaxation_methods.md`, and both docstrings); the
**remaining** modeling work is to give greedy an optional tighten-back / graded step so it
isn't systematically looser. Also `:114` ranks by `abs(shadow_price)` — relaxation should help
the objective in one direction, so rank by signed contribution; and reading duals after
`slim_optimize()` (`:95,112`) is fragile, prefer `model.optimize()` + the returned solution.

### 3.4 Finite-difference control coefficients
`limit_proteins/get_conc_control_coeffs.py:107-115` re-solves the LP per candidate protein with
a 2× bound step. The reduced cost / shadow price of the usage upper-bound from a *single* solve
is the analytic margin and collapses the loop — both faster and not direction-biased.

### 3.5 ecFVA reports an envelope, not per-reaction bounds
`utilities/ec_fva.py` nan-reduces each reaction's flux across *all* groups' optimizations, so
the reported min/max is the union envelope rather than independent per-reaction FVA. Likely
intended GECKO semantics, but should be documented; the diagonal gives true per-reaction FVA.

### 3.6 sigma fit is a 100-point brute-force grid over a monotone response
`kcat_sensitivity_analysis/sigma_fitter.py:143-150`. Growth is monotone in sigma, so bisection
finds the best sigma in ~7 solves instead of 100, each currently a full
`set_prot_pool_size` + cold solve.

### 3.7 MW from sequence: `X` residue and empty-sequence handling
`databases/mw.py:58` uses an unweighted mean for `X` (differs from MATLAB's 126.5 Da, ~7 Da/X);
`:105-106` returns 18 Da (water) for an empty sequence, which a caller could mistake for a real
protein. Consider NaN for empty input.

---

## 4. Lower-priority bugs & sanity guards — **(a)(f)** (remaining)

- `kcat_sensitivity_analysis/sensitivity_tuning.py:283` — tuned kcat has no physical ceiling
  (can exceed the ~1e7 s⁻¹ diffusion limit); `:313-322` recovers the pre-tune kcat by
  string-parsing the `notes` field — store it in a numeric dict instead. (A `max_iterations`
  cap and stall tolerance were already added.)
- `limit_proteins/fill_enz_concs.py:91` — same NaN/negative-concentration concern as
  `constrain_enz_concs` (which now skips negatives); a guard here would be consistent.
- `expand.py` objective-coefficient / GPR-explosion guard — **intentionally not changed**
  (dropped by request).

Fixed this round: `constrain_enz_concs` negative-conc guard, `calculate_f_factor` no-match
warning, `add_new_rxns_to_ec` id-collision pre-check, `load_conventional_gem` json/mat dispatch,
`run_dlkcat` timeout→RuntimeError, BRENDA SA/MW sanity bounds.

---

## 5. Organization / API hygiene — **(c)(d)(f)** (remaining)

- Two relaxation algorithms (`flexibilize_enz_concs`, `relax_proteomics_greedy`) — **evaluated**
  (`docs/relaxation_methods.md`): they are complementary, not duplicates, and should stay
  separate. Remaining nicety: align their failure conventions (one returns `converged=False`,
  the other raises).
- The deprecated-alias surface (10 in `__init__.py`, plus private `_…` aliases kept in
  `fuzzy_kcat_matching/__init__.py:43-52`) is large; plan to drop on the next major and ensure
  tests import the public names.

Fixed this round: deduped `_lookup_flux` and the gene-to-index builders, standardised the
`ec_fva` solve check on `sol.status`, removed the commented-out adapter auto-discovery block.
(The `_REV` canonicalization was unified into `canonicalize_rxn_id` earlier.)

---

## Suggested order of attack (remaining)

1. Modeling-intent items (section 3) — discuss as design decisions; you are the original author.
2. Lower-priority guards (section 4) and the remaining organization items (section 5).
