# raven-python integration

geckopy depends on [raven-python](https://github.com/SysBioChalmers/raven-python),
the Python port of the RAVEN Toolbox. raven-python owns the generic
model-manipulation primitives geckopy used to ship in-tree; geckopy keeps the
enzyme-constraint layer that's specific to it.

## What geckopy delegates to raven-python today

| geckopy module / API | re-exports from raven-python |
|---|---|
| `geckopy.ec_model.pipeline.expand.expand_model` | `raven_python.manipulation.expand.expand_model` |
| `geckopy.ec_model.pipeline.preprocess.convert_to_irreversible` | `raven_python.manipulation.irreversible.convert_to_irreversible` |
| `geckopy.utilities.add_new_rxns_to_ec` (internal `_gpr_to_dnf`) | `raven_python.manipulation.expand._gpr_to_dnf` |

These are 1:1 re-exports — same function names, same signatures. The geckopy
modules stay as the documented entry points so existing user code keeps
working; over time, callers can `import` directly from raven-python if they
prefer.

## Planned future migrations

| geckopy module / API | candidate raven-python module | why deferred |
|---|---|---|
| `geckopy.utilities.ec_fseof` | `raven_python.analysis.fseof` (base FSEOF) | raven-python's `fseof` is a redesign: regression-based target selection (vs strict monotonicity), pFBA per step (vs FBA), correlation threshold (vs top-25%-by-slope). geckopy's `ec_fseof` mirrors GECKO MATLAB's `ecFSEOF.m` exactly. Migrating means either changing geckopy's documented behaviour or duplicating raven-python's primitives — both bigger than a re-export. Tracked for a post-alpha PR. |
| `geckopy.io` (no YAML reader/writer of its own — uses cobrapy) | `raven_python.io.yaml` (cobrapy YAML + `ec_*` enzyme-constrained fields) | This is a *new feature* geckopy could adopt, not a migration. Once geckopy starts persisting `ec` data in YAML (rather than via SBML notes), it can use raven-python's reader/writer directly. |

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
