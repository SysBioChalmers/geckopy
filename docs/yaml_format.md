# geckopy ecModel YAML format

This document defines the YAML format that geckopy uses to save and
load ecModels.

## Why a custom format

Standard cobrapy already has its own YAML format for ordinary
metabolic models. But an ecModel carries extra information that
cobrapy doesn't know about: a list of kcat values, a list of
enzymes with molecular weights and protein concentrations, and
the matrix that links reactions to the enzymes that catalyse them.

We had two options:

1. Match the legacy MATLAB / RAVEN YAML format used by GECKO MATLAB
   at the time, so files exchange perfectly between MATLAB and
   Python.
2. Match cobrapy's YAML format, then add the ec-specific data as
   extra top-level keys.

We picked option 2. cobrapy already serialises metabolic models in
a clean schema; reusing it means any cobrapy-aware tool (Escher,
Memote, etc.) can load the cobra portion of a geckopy file and
just ignore the GECKO extras. MATLAB GECKO and RAVEN (via
`writeYAMLmodel` / `readYAMLmodel`) now emit and accept the same
format, so ecModels exchange directly between the two toolboxes
with no translator. RAVEN's reader still loads the older
`---` / `!!omap` files as well, so existing distributions keep
working.

## Design goals

1. **Compatible with cobrapy.** The cobra-shaped portion of the
   file is exactly what `cobra.io.dict.model_from_dict` produces;
   loading it doesn't need a custom parser.
2. **Round-trippable with MATLAB GECKO.** RAVEN's `writeYAMLmodel`
   and `readYAMLmodel` emit and accept this same format, so
   ecModels move between the two toolboxes without conversion.
3. **Human-readable.** Plain YAML mappings, no `!!omap` tags. Empty
   strings and NaN values can be omitted; the reader fills in the
   defaults.

## Top-level structure

The whole file is one YAML mapping (a top-level dictionary). It
contains the following keys:

| Key | Required? | Owner | What it holds |
|---|---|---|---|
| `id` | yes | cobra | Model id |
| `name` | optional | cobra | Human-readable name |
| `version` | optional | cobra | cobra schema version (omit if unused) |
| `compartments` | optional | cobra | Mapping of compartment id -> name |
| `metabolites` | yes | cobra | List of metabolite entries (see below) |
| `reactions` | yes | cobra | List of reaction entries (see below) |
| `genes` | yes | cobra | List of gene entries (see below) |
| `ec-rxns` | yes (for ecModel) | geckopy | Per-reaction ec data (see below) |
| `ec-enzymes` | yes (for ecModel) | geckopy | Per-enzyme ec data (see below) |
| `gecko_light` | optional | geckopy | Boolean flag; defaults to `false` |
| `metaData` | optional | geckopy | Free-form provenance: version, date, author, taxonomy, note |

The "Owner" column says which library the field belongs to. Tools
that only understand the cobra side (cobrapy, escher, memote, ...)
silently ignore the geckopy-specific keys (`ec-rxns`, `ec-enzymes`,
`gecko_light`, `metaData`), so a geckopy YAML loads cleanly as a
plain cobra model in those tools — you just lose the ec layer.

## Cobra-shaped section

This is exactly the schema produced by
`cobra.io.dict.model_to_dict`. A few conventions worth noting:

- `compartments` is a flat mapping (`{c: cytosol, e: extracellular}`),
  not a list.
- `metabolites`, `reactions`, `genes` are lists of plain mappings
  (one per metabolite/reaction/gene).
- A reaction's `metabolites` field is a flat mapping
  `met_id -> stoichiometric_coefficient`.
- `annotation` is a mapping `key -> [list, of, strings]`. RAVEN's
  per-metabolite `smiles` field lives at
  `annotation: {smiles: ["..."]}`.

Example metabolite:

```yaml
- id: s_0001
  name: "(1->3)-beta-D-glucan"
  compartment: ce
  formula: C6H10O5
  charge: 0
  annotation:
    bigg.metabolite: ["13BDglcn"]
    chebi: ["CHEBI:37671"]
    kegg.compound: ["C00965"]
    metanetx.chemical: ["MNXM6492"]
    sbo: ["SBO:0000247"]
    smiles: ["C(C1C(C(C(C(O1)O)O)O)O)O"]
```

Example reaction:

```yaml
- id: r_0001
  name: "(R)-lactate:ferricytochrome-c 2-oxidoreductase"
  metabolites:
    s_0027: -1.0
    s_0556: 1.0
  lower_bound: 0
  upper_bound: 1000
  gene_reaction_rule: "(YDL174C and YML054C) or (YEL024W and YML054C)"
  subsystem: "Oxidative phosphorylation"
  annotation:
    ec-code: ["1.1.2.4"]
    kegg.reaction: ["R00196"]
```

## `ec-rxns`: per-reaction ec data

A list, one entry per catalysed reaction. When the build pipeline
splits a reaction across isozymes (`_EXP_<N>` suffix) or into
forward/reverse pairs (`_REV` suffix), each variant gets its own
entry here. Fields:

| Field | Type | Required? | Notes |
|---|---|---|---|
| `id` | string | yes | Matches a model reaction id; may carry `_EXP_<N>` and/or `_REV` suffixes |
| `kcat` | number | yes | Turnover number in s^-1; write NaN as `.nan` |
| `source` | string | optional | Where the kcat came from (`"brenda"`, `"dlkcat"`, `"manual"`, ...). Default `""` |
| `notes` | string | optional | Free-form note. Default `""` |
| `eccodes` | string \| list | optional | One EC code as a string (`"1.1.2.4"`), or a list when several apply. Default `""` |
| `enzymes` | mapping `enzyme_id -> stoich` | yes | Which enzymes catalyse this reaction, and the subunit count of each. Every key must appear in `ec-enzymes` |

Example:

```yaml
- id: r_0001_EXP_1
  kcat: 1500.0
  source: "brenda"
  eccodes: ["1.1.2.4", "1.1.99.40"]
  enzymes:
    P00045: 1
    P32891: 1
```

## `ec-enzymes`: per-enzyme ec data

A list, one entry per unique enzyme. Fields:

| Field | Type | Required? | Notes |
|---|---|---|---|
| `genes` | string | yes | Gene id; must match a `genes[].id` from the cobra section |
| `enzymes` | string | yes | Enzyme id (usually a UniProt accession) |
| `mw` | number | optional | Molecular weight in Da. Write NaN as `.nan`. Default NaN |
| `sequence` | string | optional | Amino acid sequence. Default `""` |
| `concs` | number | optional | Measured concentration in mmol/gDW (from proteomics). Write NaN as `.nan`. Default NaN |

The `standard` pseudo-gene (added by `get_standard_kcat` to cover
reactions without a real gene assignment) appears as both `genes:
"standard"` and `enzymes: "standard"`.

Example:

```yaml
- genes: "Q0045"
  enzymes: "P00401"
  mw: 58798.0
  sequence: "MVQRWLYSTNAKDIAVLY..."
```

## Sparse defaults

To keep files compact, writers can omit any field that takes the
documented default (empty string, NaN). Readers fill the defaults
back in when loading. The exception is numeric NaN that you want
to keep explicit (e.g. "this kcat is deliberately unknown"): write
it as `.nan` so the reader sees it.

## Legacy MATLAB / RAVEN format

Older RAVEN releases wrote YAML in a different on-disk layout
(outer `---` / `!!omap` wrapper, `!!omap` tags throughout,
`id`/`name`/`geckoLight` nested under `metaData`, scalar
annotation values, ...). Files in that older shape still load via
RAVEN's `readYAMLmodel` (it auto-detects the layout); geckopy's
own loader only accepts the cobrapy-style format, so if you need
to bring an older RAVEN/GECKO YAML into geckopy, re-save it once
from MATLAB (`writeYAMLmodel` now produces the cobrapy layout) or
run a one-off conversion script.

For reference, the cosmetic differences between the two layouts:

| Aspect | Legacy MATLAB / RAVEN | Current (cobrapy / this spec) |
|---|---|---|
| Top-level shape | YAML sequence of single-key mappings | YAML mapping |
| Ordered-map tags | `!!omap` everywhere | none |
| `compartments` | sequence of single-key mappings | flat mapping |
| `metabolites[]`, `reactions[]`, `genes[]` | each entry is `!!omap` | each entry is a plain mapping |
| `reactions[].metabolites` | `!!omap` list | flat mapping `met_id -> stoich` |
| `ec-rxns[].enzymes` | `!!omap` list | flat mapping `enzyme_id -> stoich` |
| Top-level `id` / `name` | nested under `metaData` | top-level |
| `metaData` provenance fields | inside `metaData` | inside `metaData` (unchanged) |
| `metaData.geckoLight` | inside `metaData`, string `"true"`/`"false"` | top-level `gecko_light: <bool>` |
| `smiles` per metabolite | top-level metabolite key | nested under `annotation` as `{smiles: ["..."]}` |
| Annotation values | scalar strings | list of strings |
| Reaction `ec-code` annotation | not exposed (only `ec-rxns[].eccodes`) | exposed at `reactions[].annotation.ec-code` if known |

The two layouts hold the same information; the difference is
purely in how it is shaped on disk.
