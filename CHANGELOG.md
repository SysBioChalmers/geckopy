# Changelog

All notable changes to **geckopy** are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project uses
[PEP 440](https://peps.python.org/pep-0440/) pre-release versioning.

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

- **raven-python is now a hard dependency** — geckopy delegates
  `expand_model` and `convert_to_irreversible` to raven-python (the
  functions originated in geckopy and were adopted upstream as their
  canonical home). raven-python is not yet on PyPI; install both via
  the git URLs in the README. See
  [`docs/raven_integration.md`](docs/raven_integration.md) for the
  current delegation and the planned future migrations.

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
