# yeast-GEM example adapter

An example ecModel project folder for _Saccharomyces cerevisiae_, demonstrating
the standard geckopy folder layout.

## Folder layout

Every geckopy ecModel project follows the same convention:

```
yeast-GEM/
├── model_adapter.toml    # parameters (required)
├── adapter.py            # custom adapter subclass (optional)
├── data/                 # organism-specific input data
│   └── uniprotConversion.tsv   # optional gene ID mapping table
├── models/               # input and output SBML/YAML model files
│   └── yeast-GEM.xml     # the conventional GEM this ecModel is built from
└── output/               # generated files, logs, and intermediate results
```

Only `model_adapter.toml` and `models/` are strictly required. The `data/` and
`output/` folders are created when first needed.

## Using this adapter

```python
from geckopy import ModelAdapter

adapter = ModelAdapter.from_folder("examples/yeast-GEM")
print(adapter.params.org_name)
# Saccharomyces cerevisiae
```

If you need custom behavior (for example, a non-standard way to identify
spontaneous reactions), create an `adapter.py` file in this folder with a
`ModelAdapter` subclass and import it explicitly:

```python
from examples.yeast_GEM.adapter import YeastAdapter

adapter = YeastAdapter.from_folder("examples/yeast-GEM")
```

## Parameters

See `model_adapter.toml` for the full list of parameters. The schema is
defined in `geckopy.adapter.params.ModelParameters` and validated on load,
so typos and unknown fields will raise a clear error.

## Note

This folder currently contains only the adapter configuration. The actual
`models/yeast-GEM.xml` file is not included in the repository and must be
provided separately before running the full ecModel construction pipeline.
