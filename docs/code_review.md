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
The entire high-priority tier (former §1.1–1.8), plus §2.8, §2.9 and the import-hygiene
bullets of §5, have been addressed and pruned. Remaining IDs are kept stable.

---

## 2. Medium-priority (correctness / robustness)

### 2.1 `_should_overwrite` mishandles NaN — **(a)**
`gather_kcats/select_kcat_value.py:205-215`. `is_unset = current == 0` is `False` for NaN, and
`new > NaN` is `False`, so a NaN kcat is never overwritten under `overwrite=False`/`"if_higher"`
— directly contradicting the docstring ("0/NaN counted as infinitely small"). Fix:
`is_unset = (current == 0) or np.isnan(current)`.

### 2.2 DLKcat parser aborts on a single unknown substrate — **(a)(c)**
`gather_kcats/read_dlkcat_output.py:109-120`. One substrate name absent from
`model.metabolites` raises and discards a multi-thousand-row prediction file. Fix: drop unknown
rows with a warning (as already done for non-numeric kcat).

### 2.3 Currency-pair stripping uses a stale snapshot — **(a)**
`gather_kcats/write_dlkcat_input.py:327-334`. For reactions with multiple currency pairs each
pair's "are substrates left?" guard tests a column snapshot taken before any removal, so
cumulative stripping can reduce a reaction to zero substrates. Fix: test against the live
`s_dense[:, j]` after each provisional removal.

### 2.4 `set_kcat_for_reactions` base-name match also hits `_REV` — **(a)**
`ec_model/pipeline/set_kcat.py:96-99` strips only `_EXP_<n>`, not `_REV`, so a base-name request
sets the same kcat on forward and reverse arms, contradicting the documented
forward/reverse-distinct policy. Fix: exclude `_REV` unless explicitly requested, or document.

### 2.5 Global `.replace(_REV_SUFFIX, "")` can mangle IDs — **(a)**
`utilities/ec_fva.py:233`, `map_rxns_to_conv.py:153`, `get_subset_ec_model.py:229`. A reaction
whose base ID contains the substring `_REV` is corrupted. Fix: strip only a trailing suffix and
handle the `_REV_EXP_` infix explicitly; extract one shared `canonicalize_rxn_id` helper (the
logic is triplicated).

### 2.6 ec_fseof feasibility / threshold bugs — **(a)(e)**
`utilities/ec_fseof.py`: feasibility tested via `sol.fluxes is None` (`:151,167,192,384`) — use
`sol.status != "optimal"`; strict-`==0` KO/KD classification (`:227`) and strict-monotonicity
on `np.diff` (`:222-227`) should use a tolerance; the carbon-source warning at `:155` compares
a negative `lower_bound` against `cs_flux` in a way that likely fires spuriously.

### 2.7 `pfba_enzymes` minimizes molar usage, not protein mass — **(e)**
`utilities/pfba_enzymes.py:83`. All usage-reaction objective coefficients are set to `1.0`, so
it minimizes total molar enzyme usage, while the docstring claims a mass-weighted (MW)
objective. Fix: weight by `ec.mw`, or correct the docstring. Also `:84-85` proceeds after a
possibly-`None` `slim_optimize` without checking.

### 2.10 `get_subset_ec_model` leaves orphaned enzyme machinery — **(a)(e)**
`utilities/get_subset_ec_model.py:141-192`. When a gene is trimmed, the enzyme is dropped from
`ec.enzymes` but its `usage_prot_<X>` reaction and `prot_<X>` metabolite are preserved, leaving
orphaned constraints. Several ec-array trims (`:174,176,213`) silently skip on a shape mismatch
instead of raising (`ec.validate()` is now called by `make_ec_model`, but subset extraction
builds a new model and should validate too).

### 2.11 OKP / network clients lack retry/transient-error handling — **(a)**
`gather_kcats/open_kinetics_predictor/client.py` (whole file): timeouts but no handling of
`ConnectionError`/`Timeout`/transient 5xx during an up-to-1-hour poll loop. `pubchem.py:202`
caches HTTP 500 as a permanent "no SMILES" miss. `pubchem.py:31-34` still requests the
deprecated `CanonicalSMILES` property (PubChem now prefers `SMILES`).

### 2.12 Redundant large copies / per-iteration set rebuilds — **(b)**
- `io/sbml.py` (`_annotate_ec_metadata`) — `model.copy()` followed by a separate
  `deepcopy(model.ec)`; the `ec` deepcopy is unnecessary on write.
- `io/sbml.py` — the file is parsed by libsbml a second time after cobra already
  parsed it; group membership could come from `cobra_model.groups`.
- `ec_model/pipeline/protein_pool.py:88,192` — `{m.id for m in model.metabolites}` /
  `{r.id …}` rebuilt inside the enzyme loop (O(n_enzymes·n_metabolites)); hoist out.
- `apply_kcat.py:108-111` — the prot-metabolite set is rebuilt by scanning all metabolites on
  every call; derive it from `ec.enzymes`.
- `ec_fva.py:266-273` — objective reset by looping over every reaction per solve.

### 2.13 Adapter required where the type says optional — **(c)(e)**
`io/sbml.py:152-182` and `adapter/resolve.py:23-66`. `read_sbml_ec_model(adapter=None)` is typed
optional but raises if `None`; `resolve_adapter` makes a full `ModelAdapter` mandatory across
30+ functions even when only one path/param is needed, forcing a full project folder for a quick
read/inspect. Consider letting functions take the specific param and fall back to the adapter.

### 2.14 `BayesianParams` parallel-list coupling unvalidated — **(a)**
`adapter/params.py:56-73`. `kcat_sources` and the five `*_source` lists are positionally
coupled with no validator enforcing equal length; an edited TOML with a mismatched list yields a
silent index error downstream. Add a `@model_validator`. Most `ModelParameters` physical fields
(`sigma`, `f`, `p_tot` ∈ [0,1], `gr_exp` > 0) also lack `Field` range constraints.

### 2.15 Writable data defaults inside the install tree — **(a)(f)**
`adapter/adapter.py:170-177` (`get_brenda_db_folder` → `<pkg>/data/brenda`) and `cli.py:66-71`
(`brenda-refresh` output). Writing into site-packages fails on read-only/zip installs. Default
to a user cache dir (e.g. `platformdirs`) or the project `path/data`, consistent with
`get_phyl_dist_path`.

### 2.16 Other data-loader robustness — **(a)**
- `uniprot_loader.py:200-225` — duplicate-ID detection only collapses *consecutive* runs, so a
  genuine adjacent file duplicate is misread as a split repetition and not flagged.
- `pax_db_loader.py:113` — no NaN-token handling; `NA`/`-` cells raise and silently drop a
  measured protein.
- `brenda/parse.py:151-158` — verify `proteins[]` PID key type matches `proteins_map` keys
  (JSON keys are strings); a type mismatch silently drops all organisms for a measurement.
- `dlkcat_ignore_lists.py:122-125`, `flux_data.py:214-231` — several TSV
  readers assume no header / tolerate ragged rows silently; a header row or `(rxn)`-less flux
  column produces wrong constraints with no warning. (KEGG download's fewer-entries handling
  was fixed.)

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
`utilities/ec_fva.py:151-206,286-294` nan-reduces each reaction's flux across *all* groups'
optimizations, so the reported min/max is the union envelope rather than independent per-reaction
FVA. Likely intended GECKO semantics, but should be documented; the diagonal gives true
per-reaction FVA.

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
  `report_enzyme_usage.py:286-293`; the `_REV` canonicalization triplicated (see 2.5); two
  `_build_*_gene_to_indices` near-duplicates in `ec_from_database.py:211-270`.
- Infeasibility checks are inconsistent across the package (`sol.status` vs `sol.fluxes is None`
  vs `sol.objective_value` vs nothing) — standardize on `sol.status == "optimal"`. (The
  flexibilize/tuning initial-solve cases were fixed; the others remain — ties to 2.6, 3.3.)
- `map_rxns_to_conv.py:25` and `select_kcat_value.py:8` still have a misplaced constant import /
  unused `Optional` import respectively (the rest were hoisted).
- `adapter/adapter.py:92-116` — 25 lines of commented-out auto-discovery shipped in the module.
- Two near-duplicate relaxation algorithms (`flexibilize_enz_concs`, `relax_proteomics_greedy`)
  with different result dataclasses and failure modes (one returns `converged=False`, the other
  raises) — consider a shared core with a pluggable ranking strategy.
- The deprecated-alias surface (10 in `__init__.py`, plus private `_…` aliases kept in
  `fuzzy_kcat_matching/__init__.py:43-52`) is large; plan to drop on the next major and ensure
  tests import the public names.

---

## Suggested order of attack (remaining)

1. The kcat-data bugs (2.1–2.4) — they change the numbers the model produces.
2. The `_REV` canonicalization + ec_fseof feasibility bugs (2.5, 2.6).
3. Network/loader robustness (2.11, 2.15, 2.16).
4. Inefficiencies and API/validation (2.12, 2.13, 2.14).
5. Modeling-intent items (section 3) — discuss as design decisions; you are the original author.
