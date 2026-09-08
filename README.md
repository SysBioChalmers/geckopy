# geckopy - under active development

[![tests](https://github.com/SysBioChalmers/geckopy/actions/workflows/test.yml/badge.svg)](https://github.com/SysBioChalmers/geckopy/actions/workflows/test.yml)

Enzyme-constrained genome-scale metabolic modelling in Python.

> **Disambiguation.** This is a new, from-scratch port of the
> [GECKO Toolbox](https://github.com/SysBioChalmers/GECKO) (MATLAB). It is unrelated to
> the earlier and separate `geckopy` Python package
> ([on PyPI as release 2.0.2 and earlier](https://pypi.org/project/geckopy/2.0.2/))
> described in Carrasco Muriel, Long & Sonnenschein, *Microbiology Spectrum*
> 11(6):e01705-23 (2023), [doi:10.1128/spectrum.01705-23](https://doi.org/10.1128/spectrum.01705-23).

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

> **Status: beta (`4.0.0b1`).** All the MATLAB GECKO 3.2.5 functions used in
> the standard ecModel build are ported, plus CMA-ES-based kcat tuning
> against experimental data. The yeast-GEM tutorial runs end-to-end.
> geckopy's version tracks the GECKO project's own toolbox-generation
> numbering (MATLAB GECKO 1–3, geckopy = 4.x), not a semantic-versioning
> API-stability claim.

## Install

```bash
pip install --pre geckopy
```

geckopy is currently published as a pre-release, so `--pre` (or an exact
pin, e.g. `geckopy==4.0.0b1`) is required until the `4.0.0` final release —
plain `pip install geckopy` won't pick it up before then.
[raven-toolbox](https://github.com/SysBioChalmers/raven-toolbox) (the Python
port of the RAVEN Toolbox) is pulled in automatically.

Optional extras:

```bash
# Adds matplotlib (needed to render plots from the tutorial)
pip install --pre "geckopy[tutorial]"

# Adds cma (needed for CMA-ES kcat tuning, geckopy.kcat_tuning.evotune)
pip install --pre "geckopy[evotune]"

# Adds pytest + ruff (only needed if you're contributing)
pip install --pre "geckopy[dev]"
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

# 4. Save in geckopy's YAML format.
save_ec_model(ec_model, "ecModel.yml", adapter=adapter)
```

Every public function (loaders, kcat fetchers, FBA/FVA helpers,
YAML I/O, ...) is reachable as `from geckopy import X` -- the
subpackages exist for code organisation but you don't have to use
them. See `geckopy.__all__` for the full list (~80 names).

The full protocol (which adds kcat curation from BRENDA and DLKcat,
proteomics integration, and Crabtree-effect simulation) is in
[`tutorials/full_ecModel/protocol.py`](tutorials/full_ecModel/protocol.py).
That script reproduces the MATLAB GECKO Nature Protocols tutorial in
Python, step by step, and is the easiest way to see geckopy in action.

## On-disk format

geckopy reads and writes ecModels as **YAML**: a canonical,
human-readable format that is a strict superset of cobrapy's YAML
schema, with three extra top-level keys for the ec data (`ec-rxns`,
`ec-enzymes`, `gecko_light`). Standard cobrapy tools (escher, memote,
...) can read the cobra portion and silently ignore the GECKO
extensions. The same format is read and written by MATLAB GECKO /
RAVEN, so ecModels exchange between the two toolboxes with no
translator. The schema and the parsing live in raven-toolbox.

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
