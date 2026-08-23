# Changelog

All notable changes to **geckopy** are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project uses
[PEP 440](https://peps.python.org/pep-0440/) pre-release versioning.

## [Unreleased]

Parity pass against the MATLAB GECKO Toolbox: geckopy's core functions now
reproduce every value GECKO's own unit-test suite pins down.

### Added

- **`tests/test_gecko_matlab_parity.py` — a 1:1 port of GECKO's
  `test/unit_tests/geckoCoreFunctionTests.m`** (47 tests, MATLAB test cases
  `tc0001`–`tc0013`), run against the same `ecTestGEM` fixture MATLAB uses.
  Every expected value is copied verbatim from the MATLAB sources rather than
  re-derived from geckopy's output, so the file is the executable definition of
  "geckopy builds the same ecModel as GECKO". It needs no MATLAB to run and
  covers `makeEcModel` (full and light), `applyComplexData`, `setProtPoolSize`,
  `getECfromGEM`, `getECfromDatabase`, save/load round-trip,
  `fuzzyKcatMatching`, `writeDLKcatInput`, `mergeDLKcatAndFuzzyKcats`,
  `selectKcatValue`, `applyKcatConstraints`, `getKcatAcrossIsozymes`,
  `applyCustomKcats`, `findMetSmiles` and the proteomics-integration chain.

### Fixed

- **`make_ec_model` sorts identifiers for full models**, as MATLAB
  `makeEcModel` does (`sortIdentifiers`, line 203, full models only, before the
  protein pseudoreactions are appended). geckopy kept the input GEM's reaction
  order, so a full ecModel came out with the same reactions in a different
  order from MATLAB's — visible in the saved YAML and in every downstream table
  that follows `ec.rxns` order.
- **`write_dlkcat_input` emits rows in reaction order**, matching MATLAB's
  column-major `find(clearedRedS < 0)`. numpy's `where` is row-major, so
  geckopy grouped the DLKcat input by substrate instead of by reaction.
- **Reverse (`_REV`) reactions keep their EC code**, via a fix in
  raven-toolbox's `convert_to_irreversible` (it now copies annotations,
  subsystem and notes onto the reverse reaction, as MATLAB's `convertToIrrev`
  copies `eccodes` / `rxnMiriams` / `subSystems` / `rxnNotes`). Without it
  `fill_eccodes_from_gem` returned `''` for every `_REV` reaction, fuzzy BRENDA
  matching found no kcat, and reverse directions were left unconstrained.
  **Requires raven-toolbox at or after that fix.**

### Notes

- `merge_kcats` is n-ary and concatenates surviving rows list by list, so the
  equivalent of MATLAB's `mergeDLKcatAndFuzzyKcats(dlkcat, fuzzy)` — which
  always emits fuzzy rows first regardless of argument order — is
  `merge_kcats(fuzzy, dlkcat, ...)`. The deprecated
  `merge_dlkcat_and_fuzzy_kcats` alias is `merge_kcats`, so it follows the
  geckopy convention, not MATLAB's signature. Row order has no numerical
  effect: `apply_kcat_list` aggregates per reaction.

## [0.2.1] — 2026-07-16

Bugfix release. `make_ec_model` no longer mutates the GEM it is given.

### Fixed

- **`make_ec_model` leaves the input GEM unchanged.** Stages 1–5
  (`remove_pseudoreaction_gprs`, `invert_backwards_only_reactions`,
  `convert_to_irreversible`, `expand_model`) ran in place on the caller's
  model, leaving it with `_REV` / `_EXP_` reaction splits, `prot_`
  pseudometabolites, and stripped pseudoreaction GPRs. The pipeline now works
  on an internal `model.copy()`, matching MATLAB GECKO's `makeEcModel` value
  semantics. (#32)
- `test_loads_sbml_gem` on Windows: the test wrote a backslash path into the
  adapter TOML, which TOML reads as an escape sequence.

### Documentation

- README: disambiguation note stating that this geckopy is a from-scratch port
  of the MATLAB GECKO Toolbox, unrelated to and sharing no code with the
  separate PyPI `geckopy` package described in Carrasco Muriel et al. (2023).

### Internal

- **CI: leaner test matrix** — full Python range (3.11–3.13) on Linux plus a
  single Python (3.12) on macOS and Windows, instead of the full 3×3
  OS × Python grid. Same OS coverage at roughly half the runner cost. (#31)
- Removed unused test imports flagged by ruff (F401).

## [0.2.0] — 2026-06-14

Dependency rename + CI hardening. geckopy now targets **raven-toolbox**
(the renamed raven-python) and runs its test matrix on every pull request.

### Changed

- **Depend on `raven-toolbox`** (renamed from `raven-python`): all imports are
  now `raven_toolbox.*`, and the `pyproject.toml` git-URL dependency and project
  URL point at `SysBioChalmers/raven-toolbox`. (#29)
- Use the public `raven_toolbox.manipulation.expand.gpr_to_dnf` instead of the
  former private `_gpr_to_dnf`. (#26)
- **CI: run the test matrix on every pull request**, not only PRs into
  `main` / `master`.
- **CI: bump GitHub Actions to Node 24 runtimes** — `checkout@v5`,
  `setup-python@v6`, `upload-artifact@v7`, `download-artifact@v8`. (#27)

### Documentation

- Rewrote `migrating_from_gecko_matlab.md` as a complete reference (#25); moved
  internal planning notes into `docs/internal/` (#28).

## [0.1.0a3] — 2026-05-30

Build-config hotfix. v0.1.0a2 was tagged but its release.yml build
failed: hatchling refuses PEP 508 direct-URL dependencies
(`raven-toolbox @ git+https://...`) by default. The v0.1.0a2 tag and
release were deleted; v0.1.0a3 is the published artifact carrying the
0.1.0a2 content described below. No code changes vs. v0.1.0a2.

### Fixed

- `pyproject.toml`: `[tool.hatch.metadata] allow-direct-references =
  true` so hatchling accepts the `raven-toolbox` git URL dependency.
  PyPI publish stays disabled (geckopy's name is taken, raven-toolbox
  isn't on PyPI yet); the direct reference is acceptable for the
  git-install path the alpha series uses.

## [0.1.0a2] — 2026-05-30

Second public alpha. Layering pass: the bits of geckopy that aren't
GECKO-specific moved into raven-toolbox, and `gecko-light` ecModels are now
buildable end-to-end. Repo housekeeping for a discoverable first GitHub
release (LICENSE, planning docs out of the user-facing tree).

### Added

- **gecko-light end-to-end.** `make_ec_model(gecko_light=True)` now produces
  a working light ecModel: cobra reactions stay singular; per-isozyme
  coupling lives in `ec.rxns` as duplicate rows with a `###_` counter
  prefix; only the shared `prot_pool` constrains enzyme usage. `apply_kcat_constraints`
  picks the lowest-cost (smallest `MW_sum/kcat`) isozyme per cobra reaction
  and writes one `prot_pool` stoichiometric coefficient. `set_kcat_for_reactions`
  recognises both base names and the `###_` prefix; the `Enzyme` proxy's
  read paths work on light models. New `tutorials/light_ecModel/protocol.py`
  mirrors MATLAB GECKO's light tutorial; opt-in Human-GEM smoke test
  (`pytest -m smoke`) builds an ecModel from the unmodified
  Human-GEM YAML in ~2.5 min. See [`docs/gecko_light_status.md`](docs/gecko_light_status.md). ([#20])

### Changed

- **YAML I/O delegated to raven-toolbox.** `save_ec_model` / `load_ec_model`
  are now ~80-LOC wrappers around `raven_toolbox.io.read_yaml_model` /
  `write_yaml_model`. raven-toolbox owns the typed `EcData` (now re-exported
  as `geckopy.EcData`), the `ec-rxns` / `ec-enzymes` / `gecko_light` YAML
  schema, and the three legacy GECKO normalisations (top-level `smiles` →
  `annotation`, reverse-direction `usage_prot_*` flip, bare-`-` document
  root). Geckopy keeps file-extension dispatch, adapter-aware path
  resolution, and provenance/diagnostics. SBML ecModel I/O removed
  (`geckopy.io.sbml` deleted; `cobra.io.read_sbml_model` is still used for
  loading the conventional starting GEM). See [`docs/raven_integration.md`](docs/raven_integration.md). ([#19])

- **`ec_fseof` re-aligned with `raven_toolbox.analysis.fseof`.** Thin
  wrapper over raven's regression-based FSEOF: drops `usage_prot_*` from
  the scan + targets, resolves `bio_rxn` from the adapter, and emits an
  optional carbon-source consistency warning. Selection moves from MATLAB
  GECKO's strict-monotonicity + top-25%-by-slope to raven's
  `|correlation| ≥ threshold` with pFBA per step. Result type is raven's
  `FSEOFResult` (replaces `EcFseofResult`); action labels are
  `amplify` / `knockdown` / `knockout` (replaces `OE` / `KD` / `KO`); the
  per-gene `essentiality` column drops (use
  `cobra.flux_analysis.single_gene_deletion`). ([#21])

- **`raven-toolbox` pin** moved from `@main` to `@develop`. raven-toolbox's
  `main` is empty in the current iteration; all releaseable work lives on
  `develop`.

### Fixed

- **`get_conc_control_coeffs` solver-state guard.** Both `_shadow_price_coeffs`
  and `_finite_difference_coeffs` previously accepted any solve whose
  `objective_value` was non-`None` and non-NaN, including infeasible
  glpk-via-optlang solves (which return `status="infeasible"` but a
  non-NaN objective from the binding-constraint rhs). New `_solution_is_optimal`
  helper requires `sol.status == "optimal"` at every solve point.

### Removed

- `geckopy.io.sbml.{read,write}_sbml_ec_model` and the `geckopy.io` package
  (SBML ecModel I/O dropped; YAML is the supported on-disk format).
- `geckopy.EcFseofResult` (use `raven_toolbox.analysis.fseof.FSEOFResult`).
- Top-level `requirements.txt` (a `pip freeze` snapshot that duplicated
  `pyproject.toml`'s declared dependencies). pyproject is the source of
  truth.

### Internal

- Internal planning docs (`brenda_refresh_plan`, `openkineticspredictor_plan`,
  `porting_plan`, `raven_inventory`, `code_review`) moved from `docs/` to
  `docs/internal/` so the user-facing docs tree only carries reference
  documentation.
- `LICENSE` file added (MIT was already declared in `pyproject.toml` and
  classifier; the text file was missing).

## [0.1.0a1] — 2026-05-30

First public alpha. Every MATLAB GECKO 3.2.5 function used in the standard
ecModel build is ported; the yeast-GEM tutorial runs end-to-end.

### Foundations

- **Project skeleton** — `EcModel` / `EcData` / `Enzyme` core classes built on
  cobrapy; `ModelParameters` / `ModelAdapter` scaffolding (TOML-driven config,
  template generator, `geckopy init` CLI). ([#1])
- **`make_ec_model` pipeline** — 9-stage MATLAB-port pipeline that builds an
  ecModel from a conventional GEM + adapter. Fixes a unit bug where MW values
  were treated as kDa instead of Da. ([#2])
- **Kcat constraints + protein-complex / MW** — `apply_kcat_constraints`,
  `set_kcat_for_reactions`, `apply_custom_kcats`, isozyme averaging,
  Complex Portal loader, `find_met_smiles`, `calculate_mw`. ([#3])

### Kcat sources

- **BRENDA + DLKcat pipelines** — EC-code resolution helpers, `fuzzy_kcat_matching`
  MATLAB-port against BRENDA, DLKcat wrappers (`write_dlkcat_input`, `run_dlkcat`,
  `read_dlkcat_output`), standard-kcat fallback. ([#4])
- **Wide BRENDA snapshot** — `kcat.tsv` / `sa.tsv` / `mw.tsv` ship both `max`
  and `median` per `(ec, substrate, organism)` triple (plus `n` =
  measurement count). `BrendaData.kcat_for("max" | "median")` picks the
  view; three adapter fields (`kcat_aggregate_brenda` /
  `kcat_aggregate_candidates` / `kcat_aggregate_isozymes`) let projects
  flip aggregation defaults globally. See
  [`docs/kcat_aggregation.md`](docs/kcat_aggregation.md). ([#18])
- **Generalised kcat merge** — `merge_kcats` records per-row provenance and
  accepts any combination of BRENDA / DLKcat / OpenKineticsPredictor /
  manual sources. ([#13])
- **Data-source downloaders** — `geckopy brenda-refresh` CLI (BRENDA bulk
  JSON), KEGG as an alternative protein/EC source (with a warning when
  fallback uses a bare KEGG gene id), OpenKineticsPredictor REST-API
  submit/fetch wrappers, YAML I/O aligned with cobrapy + GECKO keys. ([#12])

### Proteomics integration & analysis

- **Protein concentrations + sigma + sensitivity** — `calculate_f_factor` +
  `load_pax_db`, `fill_enz_concs`, `constrain_enz_concs`,
  `get_conc_control_coeffs`, `flexibilize_enz_concs`, `constrain_flux_data`;
  sigma fitting and `sensitivity_tuning`. ([#5])
- **Analysis utilities** — `load_flux_data`, `load_prot_data`,
  `enzyme_usage` + `report_enzyme_usage`, `map_rxns_to_conv`, `ec_fva`,
  `ec_fseof`, `add_new_rxns_to_ec`, `get_subset_ec_model`. ([#6])
- **Persistence** — `load_ec_model` + `save_ec_model` (YAML + SBML
  round-trip with `ec_*` metadata), `merge_dlkcat_and_fuzzy_kcats`,
  `load_conventional_gem`. ([#7])
- **Enzyme accessor + analysis utilities + SBML I/O** —
  `Enzyme` / `EnzymeView` accessor, `pfba_enzymes`,
  `get_enzyme_bottlenecks`, parallel `ec_fva` via multiprocessing,
  `relax_proteomics_greedy`, SBML reader/writer with MW_KCAT encoding,
  flatter top-level imports. ([#9])

### Refactors

- **Dedup, renames, splits** — de-duped ec-layer constants across 12
  files, renamed 7 functions for verb-style consistency (with deprecated
  aliases), flattened top-level API, `resolve_adapter` helper,
  `fuzzy_kcat_matching` split into a 4-file subpackage, SBML name
  aliases, auto-flip of legacy reverse-direction protein reactions,
  `ec.kcat` sentinel switched from `NaN` to `0`. ([#11])

### Solver / numerics

- **Bisect sigma, shadow-price control coeffs, diagonal ecFVA, median
  kcat** — `fit_sigma(method="bisect")` exploiting monotonicity;
  `calculate_mw` returns NaN on empty / all-skipped sequence;
  `get_conc_control_coeffs` switches to LP shadow-price with scipy
  auto-fallback; `ec_fva` reports the exact per-reaction range
  (diagonal); `median` option for kcat aggregation. ([#17])

### Code review hardening

- **§1 silent-data-loss guards** — `genes`-list length validation,
  numpy-scalar coercion, kcat / relaxation loop non-termination guards,
  enzyme-provenance round-trip across SBML, UniProt-accession dedup,
  Complex Portal sub-complex classification, downloader hardening,
  NaN-aware kcat overwrite, non-fatal DLKcat unknown-substrate handling,
  `ec_fseof` / `pfba_enzymes` feasibility checks, OKP GET retries,
  current PubChem SMILES property, perf hoists, BayesianParams length
  validation, robust BRENDA protein lookup. ([#14])
- **§2 adapter decoupling** — adapter-free read / inspect of ecModels;
  per-parameter overrides decouple analysis functions from the adapter;
  `bio_rxn` overrides + loader robustness; writable BRENDA data kept out
  of the install tree. ([#15])
- **§4 / §5 cleanup** — sanity guards + DLKcat ignore-list header
  tolerance; dedup helpers + standardised infeasibility check;
  dead-code removal; tuned-kcat capping + negative-concentration guards;
  dropped private back-compat aliases in `fuzzy_kcat_matching`. ([#16])

### Docs, tutorial, packaging

- **Tutorial + packaging** — full_ecModel notebook tutorial walking
  through Stages 0–5 (build → proteomics → simulation),
  `ec_fva` tqdm progress bar, ruff baseline, `.gitattributes` LF
  enforcement, PyPI-ready packaging metadata, GitHub Actions for
  tests / lint / release. ([#8])
- **Docs sweep** — README rewrite for non-specialist readability,
  friendlier docstrings on the public-API modules, Carrasco et al. (2023)
  citation, third-party library gotchas (SBML / multiprocessing / ruamel),
  RAVEN function inventory + ravenpy scope, gecko-light status. ([#10])

### Dependencies

- **raven-toolbox is now a hard dependency** — geckopy delegates
  `expand_model` and `convert_to_irreversible` to raven-toolbox (the
  functions originated in geckopy and were adopted upstream as their
  canonical home). raven-toolbox is not yet on PyPI; install both via
  the git URLs in the README. See
  [`docs/raven_integration.md`](docs/raven_integration.md) for the
  current delegation and the planned future migrations.

[0.2.1]: https://github.com/SysBioChalmers/geckopy/releases/tag/v0.2.1
[0.2.0]: https://github.com/SysBioChalmers/geckopy/releases/tag/v0.2.0
[0.1.0a3]: https://github.com/SysBioChalmers/geckopy/releases/tag/v0.1.0a3
[0.1.0a2]: https://github.com/SysBioChalmers/geckopy/releases/tag/v0.1.0a2
[0.1.0a1]: https://github.com/SysBioChalmers/geckopy/releases/tag/v0.1.0a1
[#1]: https://github.com/SysBioChalmers/geckopy/pull/1
[#2]: https://github.com/SysBioChalmers/geckopy/pull/2
[#3]: https://github.com/SysBioChalmers/geckopy/pull/3
[#4]: https://github.com/SysBioChalmers/geckopy/pull/4
[#5]: https://github.com/SysBioChalmers/geckopy/pull/5
[#6]: https://github.com/SysBioChalmers/geckopy/pull/6
[#7]: https://github.com/SysBioChalmers/geckopy/pull/7
[#8]: https://github.com/SysBioChalmers/geckopy/pull/8
[#9]: https://github.com/SysBioChalmers/geckopy/pull/9
[#10]: https://github.com/SysBioChalmers/geckopy/pull/10
[#11]: https://github.com/SysBioChalmers/geckopy/pull/11
[#12]: https://github.com/SysBioChalmers/geckopy/pull/12
[#13]: https://github.com/SysBioChalmers/geckopy/pull/13
[#14]: https://github.com/SysBioChalmers/geckopy/pull/14
[#15]: https://github.com/SysBioChalmers/geckopy/pull/15
[#16]: https://github.com/SysBioChalmers/geckopy/pull/16
[#17]: https://github.com/SysBioChalmers/geckopy/pull/17
[#18]: https://github.com/SysBioChalmers/geckopy/pull/18
[#19]: https://github.com/SysBioChalmers/geckopy/pull/19
[#20]: https://github.com/SysBioChalmers/geckopy/pull/20
[#21]: https://github.com/SysBioChalmers/geckopy/pull/21
[#31]: https://github.com/SysBioChalmers/geckopy/pull/31
[#32]: https://github.com/SysBioChalmers/geckopy/pull/32
