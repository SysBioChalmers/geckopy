# geckopy

[![tests](https://github.com/edkerk/geckopy/actions/workflows/test.yml/badge.svg)](https://github.com/edkerk/geckopy/actions/workflows/test.yml)

Enzyme-constrained genome-scale metabolic modeling in Python.

geckopy is a Python port of the [GECKO Toolbox](https://github.com/SysBioChalmers/GECKO)
(MATLAB), built on [cobrapy](https://github.com/opencobra/cobrapy). It reconstructs
*ec*Models that account for the protein cost of catalysing each enzymatic
reaction, integrating kcat values from BRENDA, DLKcat predictions, and
custom curations.

> **Status: alpha.** The library is feature-complete with respect to the
> MATLAB protocol it ports (GECKO v3.2.5), and runs the full yeast-GEM
> tutorial end-to-end. Not yet on PyPI; install directly from GitHub.

## Install

```bash
pip install git+https://github.com/edkerk/geckopy.git
```

Optional extras:

```bash
pip install "geckopy[tutorial] @ git+https://github.com/edkerk/geckopy.git"  # matplotlib for tutorial plots
pip install "geckopy[dev]     @ git+https://github.com/edkerk/geckopy.git"  # pytest + ruff
```

Python 3.11 or newer.

## Quick start

```python
from geckopy import ModelAdapter, make_ec_model
from geckopy.databases import load_uniprot_tsv
from geckopy.utilities import load_conventional_gem, save_ec_model

# Load an adapter (organism parameters live in model_adapter.toml).
adapter = ModelAdapter.from_folder("my_project")

# Load the starting GEM and the UniProt cache.
model = load_conventional_gem(adapter)
uniprot = load_uniprot_tsv(adapter.params.path / "data" / "uniprot.tsv")

# Build the ecModel (adds protein pseudometabolites, pool exchange,
# usage reactions, and populates model.ec).
ec_model = make_ec_model(model, adapter, uniprot_db=uniprot)

# Save in geckopy's canonical YAML format.
save_ec_model(ec_model, "ecModel.yml", adapter=adapter)
```

For a full end-to-end protocol (yeast-GEM, BRENDA + DLKcat kcat sources,
proteomics integration, Crabtree-effect simulation), see
[`tutorials/full_ecModel/protocol.py`](tutorials/full_ecModel/protocol.py).

## YAML format

geckopy reads and writes ecModels in a canonical YAML format that is a
strict superset of cobrapy's dict schema, with two GECKO-specific
top-level keys (`ec-rxns`, `ec-enzymes`). See
[`docs/yaml_format.md`](docs/yaml_format.md) for the spec and the
migration table from the legacy RAVEN format. Models written by
cobrapy-aware tools (escher, memote, ...) silently ignore the GECKO
extensions.

## Relationship to MATLAB GECKO

Function-by-function port of the MATLAB toolbox; algorithmic fidelity
takes priority over Pythonic idiom. Each ported source file carries a
`Ported from GECKO MATLAB: <path>` docstring header, and intentional
divergences are tagged with `MATLAB-COMPAT:` comments. The cumulative
list of MATLAB-side improvements surfaced during the port lives in
[`docs/future_improvements.md`](docs/future_improvements.md).

## License

MIT. See [`LICENSE`](LICENSE) (if present) or the `license` field of
[`pyproject.toml`](pyproject.toml).
