# geckopy

[![tests](https://github.com/edkerk/geckopy/actions/workflows/test.yml/badge.svg)](https://github.com/edkerk/geckopy/actions/workflows/test.yml)

Enzyme-constrained genome-scale metabolic modelling in Python.

## What this is

Genome-scale metabolic models (GEMs) describe an organism's chemistry as a
big network of reactions. A standard GEM lets each reaction run as fast as
its substrates allow. In reality, enzymes catalyse those reactions, and
each enzyme has a finite capacity (its kcat) and a finite supply (the
cell's total protein budget). An *ec*Model (enzyme-constrained model)
adds those limits to the network.

geckopy is a Python port of the [GECKO Toolbox](https://github.com/SysBioChalmers/GECKO)
(MATLAB). It builds on [cobrapy](https://github.com/opencobra/cobrapy) for
the standard modelling layer, then layers GECKO's enzyme-constraint machinery
on top. It fetches kcat values from BRENDA, predicts them with DLKcat,
applies custom curations, and integrates proteomics measurements.

> **Status: alpha.** All the MATLAB GECKO 3.2.5 functions used in the
> standard ecModel build are ported. The yeast-GEM tutorial runs
> end-to-end. Not yet on PyPI; install from GitHub for now.

## Install

```bash
pip install git+https://github.com/edkerk/geckopy.git
```

Optional extras:

```bash
# Adds matplotlib (needed to render plots from the tutorial)
pip install "geckopy[tutorial] @ git+https://github.com/edkerk/geckopy.git"

# Adds pytest + ruff (only needed if you're contributing)
pip install "geckopy[dev] @ git+https://github.com/edkerk/geckopy.git"
```

Requires Python 3.11 or newer.

## Quick start

A minimal end-to-end build, assuming you have a project folder set up the
way the tutorial expects (`model_adapter.toml` + `models/` + `data/`):

```python
from geckopy import (
    ModelAdapter,
    load_conventional_gem,
    load_uniprot_tsv,
    make_ec_model,
    save_ec_model,
)

# 1. The adapter holds organism-specific parameters (taxonomy id,
#    biomass reaction, average enzyme saturation, etc.). They live in
#    model_adapter.toml so they're easy to swap per organism.
adapter = ModelAdapter.from_folder("my_project")

# 2. Load the starting GEM (the regular metabolic model you want to
#    extend) and a cached UniProt query (for enzyme MW and sequence).
model = load_conventional_gem(adapter)
uniprot = load_uniprot_tsv(adapter.params.path / "data" / "uniprot.tsv")

# 3. Build the ecModel. This adds protein pseudo-metabolites, a shared
#    protein pool, and per-enzyme usage reactions on top of the GEM.
ec_model = make_ec_model(model, adapter, uniprot_db=uniprot)

# 4. Save in geckopy's YAML format. SBML (.xml) is also supported.
save_ec_model(ec_model, "ecModel.yml", adapter=adapter)
```

Every public function (loaders, kcat fetchers, FBA/FVA helpers,
SBML I/O, ...) is reachable as `from geckopy import X` -- the
subpackages exist for code organisation but you don't have to use
them. See `geckopy.__all__` for the full list (~85 names).

The full protocol (which adds kcat curation from BRENDA and DLKcat,
proteomics integration, and Crabtree-effect simulation) is in
[`tutorials/full_ecModel/protocol.py`](tutorials/full_ecModel/protocol.py).
That script reproduces the MATLAB GECKO Nature Protocols tutorial in
Python, step by step, and is the easiest way to see geckopy in action.

## On-disk formats

geckopy reads and writes ecModels in two formats:

- **YAML** (recommended): a canonical, human-readable format that is a
  strict superset of cobrapy's YAML schema, with two extra top-level
  keys for the ec data (`ec-rxns`, `ec-enzymes`). Standard cobrapy
  tools (escher, memote, ...) can read the cobra portion and silently
  ignore the GECKO-specific extensions.
- **SBML** (.xml): for interoperability with other constraint-based
  modelling tools. The ec data is encoded as a `Protein` SBML group
  with MW carried in species notes.

See [`docs/yaml_format.md`](docs/yaml_format.md) for the YAML
specification and how it differs from the legacy MATLAB / RAVEN format.

## Relationship to MATLAB GECKO

geckopy is a function-by-function port: algorithmic fidelity comes
first, Pythonic idiom second. Each ported source file carries a
`Ported from GECKO MATLAB: <path>` docstring header. Where the Python
version intentionally differs (different direction conventions, fixed
MATLAB bugs, etc.), the divergence is tagged with a `MATLAB-COMPAT:`
comment in source. The cumulative list of MATLAB-side improvements
the port surfaced lives in
[`docs/future_improvements.md`](docs/future_improvements.md).

## License

MIT. See [`LICENSE`](LICENSE), or the `license` field of
[`pyproject.toml`](pyproject.toml).
