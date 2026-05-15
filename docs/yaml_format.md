# geckopy ecModel YAML format

This document specifies the canonical on-disk YAML format for geckopy
ecModels. The format is also intended to be the target for MATLAB
GECKO once `writeYAMLmodel`/`readYAMLmodel` are updated; until then,
the MATLAB-side legacy format (RAVEN / `!!omap`) and the geckopy-side
canonical format are different. See the migration table at the bottom
of this document.

## Design goals

1. **Round-trippable with MATLAB GECKO** so that ecModels can be
   exchanged between the two implementations without lossy
   conversion. (Requires the MATLAB-side updates listed in
   [future_improvements.md](future_improvements.md).)
2. **Compatible with cobra-py's `cobra.io.dict.model_from_dict`** so
   the cobra-shaped portion of the file can be loaded by any
   cobra-py-based tool (escher, memote, etc.) with no additional
   readers, with the GECKO-specific extensions being silently ignored
   by those tools.
3. **Human-readable.** Plain YAML mappings everywhere; no `!!omap`
   tags; sparse-by-omission allowed for empty / `NaN` fields.

## Top-level structure

The document is a single YAML mapping with these keys:

| Key | Required? | Owner | Description |
|---|---|---|---|
| `id` | yes | cobra | Model id |
| `name` | optional | cobra | Human-readable name |
| `version` | optional | cobra | cobra schema version (omit if unused) |
| `compartments` | optional | cobra | Mapping of compartment id -> name |
| `metabolites` | yes | cobra | List of metabolite mappings (see below) |
| `reactions` | yes | cobra | List of reaction mappings (see below) |
| `genes` | yes | cobra | List of gene mappings (see below) |
| `ec-rxns` | yes (for ecModel) | geckopy | Per-rxn ec data (see below) |
| `ec-enzymes` | yes (for ecModel) | geckopy | Per-enzyme ec data (see below) |
| `gecko_light` | optional | geckopy | Boolean flag; defaults to `false` |
| `metaData` | optional | geckopy | Free-form provenance: version, date, author, taxonomy, note |

Tools that consume the cobra-shaped portion of the file (cobra-py,
escher, memote, ...) silently ignore the geckopy-specific top-level
keys (`ec-rxns`, `ec-enzymes`, `gecko_light`, `metaData`).

## Cobra-shaped section

Identical to the schema produced by `cobra.io.dict.model_to_dict`.
Notable conventions:

- `compartments` is a flat mapping, not a sequence.
- `metabolites`, `reactions`, `genes` are sequences of plain
  mappings.
- `reactions[].metabolites` is a flat mapping `met_id -> stoich`.
- `annotation` is a mapping `key -> [list, of, strings]`. For the
  SMILES convention used by RAVEN, the SMILES string lives at
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

## GECKO ec-rxns section

`ec-rxns` is a sequence of mappings, one per catalysed (or
isozyme-expanded, or reversibility-split) reaction. Each entry:

| Field | Type | Required? | Notes |
|---|---|---|---|
| `id` | string | yes | Matches a model reaction id; may carry `_EXP_<N>` and/or `_REV` suffixes |
| `kcat` | number | yes | s^-1; serialise NaN as `.nan` |
| `source` | string | optional | Provenance tag (`"brenda"`, `"dlkcat"`, `"manual"`, ...). Default `""` |
| `notes` | string | optional | Free-form note. Default `""` |
| `eccodes` | string \| list of strings | optional | Single EC code as string (e.g. `"1.1.2.4"`), or a list when multiple ECs apply. Default `""` |
| `enzymes` | mapping `enzyme_id -> stoich` | yes | Encodes the rxn_enz_mat row. Each `enzyme_id` must appear in `ec-enzymes` |

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

## GECKO ec-enzymes section

`ec-enzymes` is a sequence of mappings, one per unique enzyme. Each
entry:

| Field | Type | Required? | Notes |
|---|---|---|---|
| `genes` | string | yes | Gene id (matches a `genes[].id` from the cobra section) |
| `enzymes` | string | yes | Enzyme id (typically a UniProt accession) |
| `mw` | number | optional | Da; serialise NaN as `.nan`. Default NaN |
| `sequence` | string | optional | Amino acid sequence. Default `""` |
| `concs` | number | optional | Measured concentration in mmol/gDW; serialise NaN as `.nan`. Default NaN |

The `standard` pseudo-gene (added by `get_standard_kcat` to cover
reactions without enzyme assignment) appears as both `genes:
"standard"` and `enzymes: "standard"`.

Example:

```yaml
- genes: "Q0045"
  enzymes: "P00401"
  mw: 58798.0
  sequence: "MVQRWLYSTNAKDIAVLY..."
```

## Sparse / NaN handling

Writers may omit any field whose value is the documented default
(empty string, NaN). Readers must fill back in the defaults when
loading. For numerical fields where NaN must be present (e.g. a
deliberately missing kcat), serialise as `.nan` (the YAML 1.1
standard NaN representation).

## Migration table: legacy MATLAB / RAVEN -> canonical

| Aspect | Legacy MATLAB / RAVEN | Canonical (this spec) | MATLAB-side action |
|---|---|---|---|
| Top-level shape | YAML sequence of single-key mappings | YAML mapping | Drop outer sequence wrapper |
| Ordered-map tags | `!!omap` everywhere | none (relies on Python 3.7+ / YAML 1.2 ordered mappings) | Stop emitting `!!omap` tags |
| `compartments` | sequence of single-key mappings | flat mapping | Flatten |
| `metabolites[]`, `reactions[]`, `genes[]` | each entry is `!!omap` | each entry is a plain mapping | Emit plain mappings |
| `reactions[].metabolites` | `!!omap` list | flat mapping `met_id -> stoich` | Flatten |
| `ec-rxns[].enzymes` | `!!omap` list | flat mapping `enzyme_id -> stoich` | Flatten |
| Top-level `id` / `name` | nested under `metaData` | top-level | Lift `id` and `name` out of `metaData` |
| `metaData` provenance fields (`version`, `date`, `givenName`, ...) | inside `metaData` | inside `metaData` (unchanged) | No change; cobra-py ignores |
| `metaData.geckoLight` | inside `metaData`, string `"true"`/`"false"` | top-level `gecko_light: <bool>` | Move out of `metaData`, switch to native boolean |
| `smiles` per metabolite | top-level metabolite key | nested under `annotation` as `{smiles: ["..."]}` | Move into annotation |
| Annotation values | scalar strings | list of strings | Wrap each value in a single-element list |
| Reaction `ec-code` annotation | not exposed (only `ec-rxns[].eccodes`) | exposed at `reactions[].annotation.ec-code` if known | Optional: emit per-rxn `ec-code` annotation |

The conversion is purely cosmetic; no GECKO data is lost or
reinterpreted. A one-off conversion script can rewrite legacy
`ecYeastGEM.yml`-style files into the canonical format.
