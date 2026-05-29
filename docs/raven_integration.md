# raven-python integration

geckopy depends on [raven-python](https://github.com/SysBioChalmers/raven-python),
the Python port of the RAVEN Toolbox. raven-python owns the generic
model-manipulation primitives and the YAML I/O backbone; geckopy keeps
the enzyme-constraint layer that's specific to it.

## What geckopy delegates to raven-python today

| geckopy module / API | re-exports from raven-python |
|---|---|
| `geckopy.ec_model.pipeline.expand.expand_model` | `raven_python.manipulation.expand.expand_model` |
| `geckopy.ec_model.pipeline.preprocess.convert_to_irreversible` | `raven_python.manipulation.irreversible.convert_to_irreversible` |
| `geckopy.utilities.add_new_rxns_to_ec` (internal `_gpr_to_dnf`) | `raven_python.manipulation.expand._gpr_to_dnf` |
| `geckopy.utilities.save_ec_model` | `raven_python.io.yaml.write_yaml_model` (cobra-side serialisation + opaque round-trip of GECKO top-level keys) |
| `geckopy.utilities.load_ec_model` | `raven_python.io.yaml.model_from_yaml_data` (cobra-side parsing + capture of unknown top-level keys onto `model.notes['_yaml_sections']`) |

The first three are 1:1 re-exports — same function names, same
signatures. geckopy keeps the modules as the documented entry points so
existing user code keeps working; over time callers can `import` from
raven-python directly if they prefer.

The YAML I/O delegation is layered:

- **`save_ec_model`** builds the `ec-rxns` / `ec-enzymes` / `gecko_light`
  sections from `EcData` (with the GECKO conventions: omit empty
  `source` / `notes` / `eccodes`, treat `kcat == 0` as "no kcat
  assigned", omit NaN `mw` / `concs` and empty `sequence`), stashes
  them on `model.notes['_yaml_sections']` plus `metaData` on
  `model.notes['metaData']`, calls `raven_python.io.yaml.write_yaml_model`,
  then restores the caller's `model.notes` (the mutation is transient).
  Numerical coercion of the GECKO sections happens inside geckopy
  because raven-python's writer only coerces the cobra-shaped portion;
  numpy / ruamel scalars inside `_yaml_sections` would otherwise trip
  the safe-dumper.
- **`load_ec_model`** reads the raw YAML with its own `_read_yaml`
  (handles legacy bare-`-` sequences of single-key mappings), applies
  `_normalize_legacy_layout` (lifts `id`/`name`/`version` out of
  `metaData`, moves per-metabolite top-level `smiles` into `annotation`),
  hands the cleaned dict to `model_from_yaml_data`, then reads the
  GECKO sections back off `cobra_model.notes['_yaml_sections']` to build
  the typed `EcData`. Geckopy also keeps `_flip_legacy_prot_direction`
  (a geckopy-specific post-pass for older MATLAB ecModels written with
  reverse-sign protein reactions).

On-disk format is the same RAVEN/cobrapy YAML in both directions — files
written by geckopy load in raven-python (as a plain `cobra.Model` with
the GECKO sections stashed opaquely on `.notes['_yaml_sections']`), and
files written by raven-python / MATLAB RAVEN / MATLAB GECKO load in
geckopy as full `EcModel`s.

## Planned future migrations

| geckopy module / API | candidate raven-python module | why deferred |
|---|---|---|
| `geckopy.utilities.ec_fseof` | `raven_python.analysis.fseof` (base FSEOF) | raven-python's `fseof` is a redesign: regression-based target selection (vs strict monotonicity), pFBA per step (vs FBA), correlation threshold (vs top-25 %-by-slope). geckopy's `ec_fseof` mirrors GECKO MATLAB's `ecFSEOF.m` exactly. Migrating means either changing geckopy's documented behaviour or duplicating raven-python's primitives — both bigger than a re-export. Tracked for a post-alpha PR. |

## What geckopy keeps in-tree (not raven-python territory)

- The `EcModel` / `EcData` / `Enzyme` dataclasses and the
  `make_ec_model` 9-stage pipeline (GECKO-specific).
- BRENDA loaders + fuzzy kcat matching + DLKcat wrappers
  (`geckopy.gather_kcats`, `geckopy.databases.brenda`).
- Per-gene KEGG REST-API downloader (`geckopy.databases.kegg_download`).
  raven-python has a different KEGG concern: bulk-file dumps for de-novo
  reconstruction (`raven_python.reconstruction.kegg.download_kegg_dump`).
- Phylogenetic-distance matrix loading (`geckopy.databases.phyl_dist`).
  raven-python has a different concern: KEGG taxonomy parsing
  (`raven_python.reconstruction.kegg.taxonomy`).
- SBML round-trip with `ec_*` metadata encoded in notes
  (`geckopy.io.sbml`).
- Everything under `geckopy.limit_proteins`,
  `geckopy.kcat_sensitivity_analysis`,
  `geckopy.ec_model.pipeline.protein_pool`, etc. — pure
  enzyme-constraint mechanics, no general-GEM analog in RAVEN.
