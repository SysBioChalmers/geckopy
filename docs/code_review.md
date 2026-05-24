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
Addressed and pruned so far: the entire high-priority tier (former §1.1–1.8); §2.1, §2.2,
§2.5, §2.6, §2.7, §2.8, §2.9, §2.10, §2.11 (OKP/PubChem parts), §2.12, §2.13 (read path +
ec_fseof/flux-data overrides), §2.14 (validator), §2.15 (data out of install tree),
§2.16 (BRENDA protein lookup); and the import-hygiene bullets of §5. Verified **false
positives** (no change): §2.3 (currency-pair stripping mutates `s_dense` live, so cumulative
state is tracked), §2.4 (`set_kcat` base-name does not strip `_REV`, so forward/reverse stay
distinct), §2.7-objective (usage fluxes already carry MW scaling, so the unit-weighted sum
minimises mass), §2.11-pubchem-500 (PubChem returns 500 for unparseable names — caching it is
intentional and tested). Remaining IDs are kept stable.

---

## 2. Medium-priority (remaining)

### 2.13 Adapter coupling — **(c)(e)** (partly done)
`read_sbml_ec_model`/`EcModel.from_cobra` now accept `adapter=None` (inspection loads), and a
`resolve_param` helper lets functions take a specific id and fall back to the adapter;
`ec_fseof` and `apply_flux_data_constraints` use it. **Remaining:** apply the same per-param
override (mostly `bio_rxn`) to `flexibilize_enz_concs` and `sensitivity_tuning` (which already
accept the growth-rate override). The builder/data-acquisition functions legitimately need the
project folder and are lower-value to decouple.

### 2.16 Other data-loader robustness — **(a)** (remaining)
- `uniprot_loader.py:200-225` — duplicate-ID detection only collapses *consecutive* runs, so a
  genuine adjacent file duplicate is misread as a split repetition and not flagged.
- `pax_db_loader.py:113` — a non-numeric abundance (`NA`/`-`) is dropped silently; arguably
  correct (an unusable value), but a debug/info count would aid diagnosis. Low value.
- `dlkcat_ignore_lists.py:122-125`, `flux_data.py:214-231` — these TSV readers assume no header
  / tolerate ragged rows silently; a header row or a `(rxn)`-less flux column produces wrong
  constraints with no warning. (KEGG download's fewer-entries handling and the BRENDA
  protein-id lookup were fixed.)

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
(effectively unconstrained), whereas `flexibilize_enz_concs` ramps gradually. The two routines
are framed as interchangeable but produce models of very different stringency; greedy yields a
much looser model. Add a tighten-back refinement pass, or document the difference.
Also `:114` ranks by `abs(shadow_price)` — relaxation should help the objective in one
direction, so rank by signed contribution. Reading duals after `slim_optimize()` (`:95,112`) is
fragile; prefer `model.optimize()` and read from the returned solution.

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

## 4. Lower-priority bugs & sanity guards — **(a)(f)**

- `kcat_sensitivity_analysis/sensitivity_tuning.py:283` — tuned kcat has no physical ceiling
  (can exceed the ~1e7 s⁻¹ diffusion limit); `:313-322` recovers the pre-tune kcat by
  string-parsing the `notes` field — store it in a numeric dict instead. (A `max_iterations`
  cap and stall tolerance were added.)
- `ec_model/pipeline/expand.py:106-129` — expanded reactions don't copy
  `objective_coefficient`; an OR-GPR reaction in the objective is dropped from it. No guard
  against combinatorial GPR explosion.
- `limit_proteins/constrain_enz_concs.py:117` — no `conc >= 0` validation; a negative proteomics
  value becomes `ub < lb` (infeasible). Same NaN/negative concern in `fill_enz_concs.py:91`.
- `limit_proteins/calculate_f_factor.py:90-100` — returns `0.0` both when proteome mass is zero
  and when no in-model enzyme matched (a mapping bug), silently zeroing the protein budget; warn
  when nothing matched.
- `utilities/add_new_rxns_to_ec.py:161` — no pre-check that new reaction IDs (and their
  `_EXP_`/`_REV` expansions) don't collide with existing reactions; partial mutation leaves a
  half-modified model on failure.
- `utilities/load_conventional_gem.py:54-56` — no `.json`/`.mat` dispatch; a JSON GEM is fed to
  the SBML reader and fails confusingly.
- `gather_kcats/run_dlkcat.py:113-115` — `subprocess.TimeoutExpired` propagates raw instead of
  as the documented `RuntimeError`.
- `databases/brenda/parse.py:136` — kcat ceiling exists, but SA/MW have no sanity bound, so a
  bogus MW corrupts SA-derived kcats.

---

## 5. Organization / API hygiene — **(c)(d)(f)**

- `_lookup_flux` duplicated between `enzyme_usage.py:142-150` and
  `report_enzyme_usage.py:286-293`; two `_build_*_gene_to_indices` near-duplicates in
  `ec_from_database.py:211-270`. (The `_REV` canonicalization was unified into
  `canonicalize_rxn_id`.)
- Infeasibility checks are inconsistent across the package (`sol.status` vs `sol.fluxes is None`
  vs `sol.objective_value` vs nothing). The flexibilize/tuning initial-solve and ec_fseof cases
  were fixed; remaining sites (e.g. `bottlenecks`, parts of `ec_fva`) could standardise on
  `sol.status == "optimal"`.
- `adapter/adapter.py:92-116` — 25 lines of commented-out auto-discovery shipped in the module.
- Two near-duplicate relaxation algorithms (`flexibilize_enz_concs`, `relax_proteomics_greedy`)
  with different result dataclasses and failure modes (one returns `converged=False`, the other
  raises) — consider a shared core with a pluggable ranking strategy.
- The deprecated-alias surface (10 in `__init__.py`, plus private `_…` aliases kept in
  `fuzzy_kcat_matching/__init__.py:43-52`) is large; plan to drop on the next major and ensure
  tests import the public names.

---

## Suggested order of attack (remaining)

1. Finish 2.13 (bio_rxn override on flexibilize_enz_concs / sensitivity_tuning) and the 2.16
   loader-robustness bits.
2. Modeling-intent items (section 3) — discuss as design decisions; you are the original author.
3. Lower-priority guards (section 4) and the remaining organization items (section 5).
