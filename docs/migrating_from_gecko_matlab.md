# Migrating from GECKO 4 (MATLAB) to geckopy

geckopy is the Python port of the GECKO toolbox. It builds the same
enzyme-constrained models, follows the same algorithms, and writes the
same on-disk format as **GECKO 4** (MATLAB). This guide is for someone
who knows GECKO 4 in MATLAB and wants to do the same work in Python.

If you are coming from **GECKO 3**, first read
[GECKO3_to_GECKO4.md](https://github.com/SysBioChalmers/GECKO) in the
GECKO repository: GECKO 4 changed a few things (notably the
forward-direction protein reactions) that geckopy also follows.

The port's goal is **behavioural and on-disk fidelity**: same reactions,
metabolites, coefficients, kcats, and YAML/SBML output. Differences come
from being idiomatic Python and from using cobrapy instead of RAVEN.

---

## 1. At a glance

| | GECKO 4 (MATLAB) | geckopy |
|---|---|---|
| Model toolbox | RAVEN + COBRA Toolbox | cobrapy |
| Model object | RAVEN struct + `ec` field | `EcModel(cobra.Model)` with an `.ec` dataclass |
| Adapter | `ModelAdapter` classdef + `ModelAdapterManager` | `ModelAdapter` + `model_adapter.toml` (pydantic-validated) |
| Default adapter | global `ModelAdapterManager.getDefault()` | none; pass `adapter=` or set `model.adapter` |
| Return style | multiple outputs `[a,b,c] = f(...)` | single return; the model is mutated in place |
| Naming | `camelCase` | `snake_case` |
| Missing data | often a silent default | usually a raised exception |
| HTTP (BRENDA/KEGG/UniProt/OKP) | `webread`/`webwrite` | `requests` |

The two exchange ecModels directly: a model saved by either tool loads
in the other (see §6).

---

## 2. The model object

GECKO MATLAB uses a RAVEN model struct with parallel cell arrays plus a
`model.ec` struct and side cells (`model.metSmiles`, `model.eccodes`).

geckopy's `EcModel` subclasses `cobra.Model` and adds one `.ec`
attribute (an `EcData` dataclass: kcat array, enzyme list, MW array,
sequences, and the sparse reaction–enzyme coupling matrix). Per-reaction
and per-metabolite metadata lives where cobrapy keeps it:

| MATLAB | geckopy |
|---|---|
| `model.metSmiles{i}` | `metabolite.annotation['smiles']` |
| `model.eccodes{i}` | `reaction.annotation['ec-code']` (a `list[str]`) |
| `model.ec.rxnEnzMat` | `model.ec.rxn_enz_mat` (scipy sparse) |
| `model.ec.kcat`, `.mw`, `.sequence`, ... | `model.ec.kcat`, `.mw`, `.sequence`, ... |
| grRules (strings) | cobrapy GPR objects (parsed AST) |

Loaders return dataclasses / DataFrames rather than structs: e.g.
`BrendaData`, `ProtData`, and `fuzzy_kcat_matching` returns a
`pandas.DataFrame`.

---

## 3. Adapters and configuration

In MATLAB you write an adapter `classdef` and register a default with
`ModelAdapterManager`. In geckopy you create an adapter folder with a
`model_adapter.toml`, loaded via `ModelAdapter.from_folder(path)`. The
TOML is parsed into pydantic models (`ModelParameters` with nested
`kegg`, `uniprot`, `okp`, `bayesian`, ...) and validated strictly, so
config typos are caught early.

There is **no global default adapter**. Functions read `model.adapter`
or take an explicit `adapter=` argument. Where a MATLAB function takes a
`modelAdapter` only to derive a path, the geckopy version takes the path
(or pre-loaded data) directly.

**Secrets** (e.g. an OpenKineticsPredictor API key) are never stored in
the adapter. Provide them via a function argument, an environment
variable, or a git-ignored file in the project `data/` folder.

---

## 4. Calling conventions

### Mutate in place, return one thing

MATLAB functions return several outputs including a copy of the model.
geckopy mutates the model in place and returns at most one value,
logging what MATLAB returned as extra outputs.

| MATLAB | geckopy |
|---|---|
| `[model,noUniprot] = makeEcModel(...)` | `make_ec_model(...)` mutates; logs unmatched genes, annotates `rxn.notes['geckopy_warning']` |
| `[model,rxnUpdated,notMatch] = applyCustomKcats(...)` | `apply_custom_kcats(...)` returns `None`; logs the rows |
| `[model,...] = getStandardKcat(...)` | `assign_standard_kcat(...)` mutates; diagnostics via logging |
| `[minFlux,maxFlux] = ecFVA(...)` | `ec_fva(...)` returns a DataFrame indexed by reaction id |
| `[rxns,kcat,name,gpr,idx] = getReactionsFromEnzyme(...)` | `get_reactions_from_enzyme(...)` returns one DataFrame |

### Errors instead of silent defaults

Where MATLAB silently substitutes a default, geckopy tends to raise so
problems surface early — e.g. missing proteome data raises
`FileNotFoundError` instead of returning an f-factor of `0.5`; an
unknown protein id raises `ValueError` instead of returning empty.

### Names and indexing

`camelCase` → `snake_case`. A couple of deliberate renames, e.g.
`getKcatAcrossIsozymes` → `fill_kcats_from_isozymes` (kept as a
deprecated alias). MATLAB is 1-indexed; geckopy is 0-indexed, which
shows up in arguments like `data_col` (default `0`) and `condition`.

### gecko-light

Not implemented yet; light-model code paths raise `NotImplementedError`.
See [gecko_light_status.md](gecko_light_status.md).

---

## 5. A minimal workflow comparison

MATLAB (sketch):

```matlab
adapter = ModelAdapterManager.setDefault('MyAdapter.m');
model   = makeEcModel(loadConventionalGEM(), false, adapter);
kcatList = fuzzyKcatMatching(model);
model   = selectKcatValue(model, kcatList);
model   = applyKcatConstraints(model);
saveEcModel(model, 'ecModel.yml');
```

geckopy (sketch):

```python
import geckopy
from geckopy import (ModelAdapter, make_ec_model, fuzzy_kcat_matching,
                     apply_kcat_list, apply_kcat_constraints, save_ec_model,
                     load_brenda_data, load_phyl_dist)

adapter = ModelAdapter.from_folder("my_model")
model   = make_ec_model(conv_gem, adapter)         # mutates / returns EcModel
brenda  = load_brenda_data(...)
kcats   = fuzzy_kcat_matching(model, brenda, phyl_dist)   # returns a DataFrame
apply_kcat_list(model, kcats)                       # mutates in place
apply_kcat_constraints(model)
save_ec_model(model, "ecModel.yml", adapter=adapter)
```

Note the geckopy steps take pre-loaded data (BRENDA, phylogenetic
distance) explicitly rather than re-loading inside each function.

---

## 6. File exchange

ecModels move between the two tools without conversion:

- **Saving** in either tool writes the same canonical, cobrapy-style
  YAML, so a geckopy `.yml` loads in GECKO 4 and vice versa.
- **Loading** is backward compatible: geckopy reads older MATLAB
  ecModels too, and it auto-flips legacy reverse-direction protein
  reactions on load — exactly as GECKO 4's `loadEcModel` does. So a
  GECKO 3 model opened in geckopy and saved comes out as a current
  model.

The BRENDA database files use a refreshed schema. geckopy ships
`kcat.tsv` / `sa.tsv` / `mw.tsv` with both max and median rows per
(ec, substrate, organism) triple, while MATLAB GECKO still ships
only max in `max_KCAT.txt` etc. See
[kcat_aggregation.md](kcat_aggregation.md) for the rationale.

---

## 7. What's the same, what's extra, what's missing

**Same algorithms, aligned results.** The port found and fixed several
GECKO 3 bugs (EC-code assignment, sigma fitting, standard-kcat
subsystem means, duplicate EC codes, ...); GECKO 4 carries the same
fixes, so the two toolboxes agree. The full per-function divergence list
lives in the in-source `MATLAB-COMPAT:` comments
(`grep -rn "MATLAB-COMPAT:" src/geckopy/`) and in
[future_improvements.md](future_improvements.md).

**Available in both** (added during the port): KEGG as a fallback
protein/EC source, and OpenKineticsPredictor over its REST API
(`submit_open_kinetics_predictor` / `fetch_open_kinetics_predictor`).

**Not yet in geckopy:** the light ecModel build
([gecko_light_status.md](gecko_light_status.md)) and Bayesian (ABC-SMC)
kcat tuning ([bayesian_tuning_plan.md](bayesian_tuning_plan.md)).

---

## 8. Where to look when something differs

- In the geckopy source, every intentional divergence carries a
  `MATLAB-COMPAT:` comment, and each ported function's docstring opens
  with `Ported from GECKO MATLAB: <path>.` so you can find the exact
  MATLAB original.
- MATLAB-side bugs and rough edges the port found are tracked in
  [future_improvements.md](future_improvements.md).
