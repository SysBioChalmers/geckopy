# Future improvements

A running list of ideas the port surfaced that are worth doing but
weren't blocking. Three audiences:

- **geckopy users / contributors** — sections "Conventions / data
  model", "Third-party library gotchas", "API simplification",
  "get_enzyme_data subsystem".
- **MATLAB GECKO maintainers** — section "MATLAB GECKO changes".
  Each item is a bug or rough edge in MATLAB GECKO that the port
  spotted while comparing the two implementations side by side.
  Every item is also flagged with a `MATLAB-COMPAT:` comment in
  the geckopy source, so `grep -rn "MATLAB-COMPAT:" src/` gives
  the full list with file references.
- **anyone hitting a strange dependency bug** — section
  "Third-party library gotchas".

## Conventions / data model

Two cross-cutting decisions where geckopy and MATLAB GECKO ended
up different. Not blocking, but worth aligning eventually.

- **Missing-kcat sentinel.** geckopy uses `NaN` in `ec.kcat` to
  mark "no value yet". MATLAB GECKO uses `0`. Any I/O layer that
  crosses between the two has to translate (`0` -> `NaN` on read,
  `NaN` -> `0` on write). The cleaner long-term fix is to switch
  MATLAB GECKO to `NaN` too, after which no translation is needed.

- **Custom subunit stoichiometry file format.** The `stoicho`
  column in `customKcats.tsv` was meant to let users override
  subunit counts manually. The MATLAB code never actually reads
  it; geckopy also ignores it. Two options if we want to support
  it properly: a separate `customComplexes.tsv` file, or a new
  "user_overrides" section in the ComplexPortal JSON. Deferred
  until someone needs it.

## Third-party library gotchas

A scratchpad of brittle behaviours we ran into in geckopy's hard
dependencies. Each is worked around in source today; the entries
exist so the next person hitting the same wall can find the fix
fast, and so we have a list of things to retest when we upgrade.

- **`libsbml.Species.appendNotes` silently fails on cobra-written
  documents.** Returns status `-5` (`LIBSBML_INVALID_OBJECT`) with
  no exception, so the MW we tried to write into the species notes
  was lost. Tried both `<html>` and `<notes>` wrappers; same result.
  Workaround: set MW via `cobra.Metabolite.notes["mw"]` *before*
  calling `cobra.io.write_sbml_model`; cobra serialises that dict
  into SBML `<notes>` reliably. See
  `geckopy.io.sbml._annotate_mw`. Tested with libsbml 5.21.1.

- **`multiprocessing.Pool(context="spawn")` deadlocks on some WSL
  kernels.** Even a trivial 2-worker pool with a `lambda x: x*2`
  task hangs forever — no error, no timeout. The original plan for
  parallel `ec_fva` called for `spawn` because it's the only context
  available on Windows. Workaround: pick the context based on the
  platform — `fork` on POSIX (which is also faster, no module
  re-import), `spawn` on Windows. See `geckopy.utilities.ec_fva.ec_fva`.

- **`ruamel.yaml`'s `ScalarInt` / `ScalarFloat` leak out of cobra
  YAML loads.** After loading a YAML model with cobra, attributes
  like `Metabolite.charge` are ruamel scalar types (subclasses of
  int/float that carry YAML formatting metadata). The ruamel
  safe-dumper has no representer for them, so writing the model
  back out fails with `RepresenterError: cannot represent an
  object: 0`. Workaround: walk the assembled document and coerce
  everything to plain Python primitives just before writing. See
  `geckopy.utilities.save_ec_model._to_native`. This bites any
  YAML-to-YAML round-trip of a cobra-loaded model.

## API simplification

Two API rough edges that geckopy inherited from MATLAB GECKO.
Cleaning them up requires a coordinated rename or behaviour change
in both packages, so deferred for now.

- **Mode B in `apply_custom_kcats`.** The "proteins only, no rxns"
  mode in `applyCustomKcats.m` is unclear: the docstring and the
  code disagree about whether the full-match gate should fire.
  geckopy follows the code (since it's the authoritative
  behaviour); MATLAB should clean up the docstring or drop the
  mode entirely.

- **Function naming.** `get_kcat_across_isozymes` is a literal
  port of the MATLAB name, but the verb `get` is misleading — the
  function modifies `ec.kcat` in place. A clearer name would be
  `fill_kcats_from_isozymes`. Deferred to a coordinated rename
  across MATLAB and Python.

## MATLAB GECKO changes

Aimed at MATLAB GECKO maintainers. Each item is something the
geckopy port found while comparing the two implementations: bugs,
dead code, unit-comment slip-ups, or unclear behaviour. The
in-source `MATLAB-COMPAT:` comments in geckopy carry the same
notes alongside the Python implementation.

Notable items:

- Switch `prot_pool_exchange` and `usage_prot_*` reactions to the
  forward direction (positive flux for protein production).
- Sort `ec.genes` alphabetically in stage 7 of `makeEcModel`.
- Implement gene-cell splitting in the UniProt loader.
- Apply the `stoicho` column from `customKcats.tsv`, or drop it from
  the schema.
- Change source-string convention for `setKcatForReactions` to `'manual'`.
- Make `getReactionsFromEnzyme` case-sensitive.
- Forbid length-N kcat lists for un-suffixed `rxn_ids` in
  `setKcatForReactions` (strict matching rule).
- Delete `getECstring` from MATLAB GECKO. After the `findECInDB`
  refactor (geckopy works with raw `;`-joined EC tokens throughout)
  the function has no remaining callers. It also has three latent
  bugs not worth fixing in place: an accumulator footgun (callers
  must add a trailing space to `EC_set` before calling, otherwise
  tokens collide), an empty-input quirk (returns `"EC"` instead of
  `""` because `strsplit("", " ")` returns `{''}`), and no
  validation (would happily emit `ECEC1.1.1.1` for already-prefixed
  input or `ECnotanec` for junk).
- Fix the validation regex in `getECfromGEM`. The current pattern
  `(\d\.(\w|-)+\.(\w|-)+\.(\w|-)+)(;\w+\.(\w|-)+\.(\w|-)+\.(\w|-)+)*(.*)`
  substituted with `$3` returns the level-3 character of the first EC
  for any valid input (e.g. `"1.2.3.4"` -> `"3"`), and the subsequent
  `~cellfun(@isempty, ...)` then flags every non-empty result as
  invalid. As written, every non-empty EC string is silently
  discarded. Replace with a straightforward
  `^TOKEN(;TOKEN)*$` validation where `TOKEN` is the canonical
  four-level dotted EC pattern with `-` allowed in any level.
- Dedupe the `intersection` helper output in `findECInDB`. The
  helper currently produces duplicates when multiple subunits share
  the same EC (e.g. `"1.1.1.1;1.1.1.1"` for a two-subunit complex
  whose subunits both map to EC 1.1.1.1). Apply a final
  `compare_wild` pass on the intersection result, mirroring what
  geckopy does.
- Simplify `loadBRENDAdata` signature: drop the `modelAdapter` arg
  and `ModelAdapterManager.getDefault()` fallback; take the folder
  path directly. The adapter resolution belongs at the call site.
- Drop the unused 5th column from the BRENDA dump file format.
  `loadBRENDAdata` parses it as `%q` and never reads it; observed
  dumps consistently set it to `*`.
- Fix the stale `[1/hr]` comment on `SAcell{3}` in `loadBRENDAdata`.
  With the post-refactor scaling factors (SA * 1/60, MW * 1/1000)
  the unit is `[1/s]`, not `[1/hr]`.
- Resolve the `'add'` action in `getECfromDatabase`: either implement
  the commented-out `addMultipleMatches` branch (and drop the
  apologetic `% I don't understand the purpose of this, let's skip
  it for now` comment) or remove the option from the documented API.
  geckopy currently only supports `'display'` and `'ignore'`.
- Subset by `ecRxns` upfront in `getECfromDatabase` rather than
  iterating every reaction and discarding the unwanted results at
  the end. The function's own inline comment already calls this out
  ("Probably faster to subset with ecRxns in the beginning of the
  script, but this was at the moment simpler to implement").
  geckopy already does this.
- Fix the search-order/output-ranking inconsistency in
  `fuzzyKcatMatching`. The search order inside `mainMatch` tries
  org-SA (output rank 5) BEFORE any-organism-no-substrate-kcat
  (output rank 4). When both are available for a reaction, MATLAB
  returns the org-SA result and reports its origin as 5, even though
  the docstring's origin ranking implies origin 4 should win. Either
  swap the search order so any-no-subs-kcat is tried before org-SA,
  or update the docstring to match the search order. geckopy
  replicates the current MATLAB behavior with a MATLAB-COMPAT note.
- Cap the iterative-EC-escalation loop in `fuzzyKcatMatching` at 4
  wildcards. The current `while ~success` loop has no termination
  condition for tokens that never match; combined with
  `EC{k} = [EC{k}(1:dot_pos(4-wild_num)) ...]`, when `wild_num`
  exceeds 4 the indexing becomes invalid (``dot_pos(0)``) and MATLAB
  errors out. geckopy returns "no match" cleanly at that point.
- Remove the dead-code `if forceWClvl == 1` block in
  `fuzzyKcatMatching`. After the preceding `while forceWClvl > 0`
  loop, `forceWClvl` is always 0, so the `if` never fires. Either
  delete it or capture the original `forceWClvl` value before the
  loop and check that.
- Drop the unused `ids` field from the loaded `phylDistStruct` in
  `KEGG_struct`. It is parsed but never read.
- Fix `sigmaFitter` to actually apply the optimal sigma to the
  returned model. The current implementation tries 100 sigma values
  in a loop, picks the best, then returns the model with the LAST
  trial (sigma=1.0) applied, not the best. The docstring claims the
  model is adapted to the optimal sigma. Add a final
  `model = setProtPoolSize(model, Ptot, f, sigma, modelAdapter);`
  after the loop.
- Fix the wildcard branch of `findMaxValue`. The current code does
  `EC_cell{i} = EC_cell{i}(strfind(EC_cell{i},'-')-1:end);` which
  keeps the suffix from before the first `-` to the end (yielding
  e.g. `".-"` for `"EC1.1.1.-"`), then re-prepends `"EC"` to produce
  `"EC.-"` and uses that as a substring-search key. No real EC code
  starts with `"EC.-"`, so the wildcard branch matches nothing in
  practice. Replace with a prefix match on everything before the
  first `-` (e.g. `"EC1.1.1."` for `"EC1.1.1.-"`).
- Delete `updateProtPool` from MATLAB GECKO. The function has been
  obsolete since GECKO 3.2.0 (all enzymes, measured and unmeasured,
  draw from the protein pool); its sole runtime behaviour on a
  3.2.0+ model is to raise an error pointing users to
  `setProtPoolSize`. geckopy does not port it; the recommended
  replacement is `set_prot_pool_size` in `protein_pool.py`.
- Make `calculateFfactor` raise (or return NaN with a warning)
  when no proteome data is provided. Silently returning 0.5 hides
  the absence of real data; downstream calculations using this f
  factor look fine until they don't. geckopy's split design
  (`load_pax_db` raises FileNotFoundError, `calculate_f_factor`
  requires pre-loaded data) makes the missing-data case explicit.
- Fix the `all(kcatSubSystemIdx)` check in `getStandardKcat`. The
  intent is "does the reaction's subsystem appear in our subsystem
  list?", but `all` of a length-N boolean vector is true only when
  every entry is true (i.e. only when the model has a single unique
  subsystem). For any model with more than one subsystem, the
  per-subsystem mean kcat is effectively dead code and the function
  always falls back to the global standard kcat. Use `any(...)` (or
  the equivalent membership check) instead. geckopy uses the
  intended semantics.
- Fix `saveEcModel`'s default filename. The MATLAB default is
  `'ecModel'` (no extension); combined with the dispatch on
  `filename(end-3:end)` falling through to `writeYAMLmodel` for
  unknown extensions, the file is written as
  `<path>/models/ecModel` (no extension). Either default to
  `'ecModel.yml'` or auto-append `.yml` when no extension is
  given. geckopy defaults to `'ecModel.yml'`.
- Drop the `e-005` -> `e-05 ` post-processing in `saveEcModel`'s
  SBML branch. The Windows/Mac stoichiometric-coefficient
  formatting workaround predates current libSBML, which formats
  consistently across platforms. The backup-file dance is also a
  footgun: if the function errors mid-way the user is left with a
  dangling `backup.xml`. Either rely on libSBML's current output
  directly, or fix it via a libSBML option rather than text
  rewriting.
- Update `writeYAMLmodel` / `readYAMLmodel` (RAVEN) to emit and
  accept the canonical geckopy YAML schema specified in
  [yaml_format.md](yaml_format.md). The conversion is purely
  cosmetic (drop the outer sequence wrapper; drop all `!!omap`
  tags; flatten `compartments`, per-met / per-rxn / per-gene
  entries, `reactions[].metabolites`, and `ec-rxns[].enzymes`;
  lift `id` and `name` out of `metaData` to top-level; move
  `geckoLight` out of `metaData` into a top-level boolean
  `gecko_light`; move per-met `smiles` into `annotation:
  {smiles: [...]}`; coerce annotation values to single-element
  lists). Once both implementations agree on the canonical
  format, ecModels become directly exchangeable between MATLAB
  GECKO and geckopy with no translator. As an interim measure,
  ship a one-off conversion script (`legacy_to_canonical.py` or
  similar in `docs/` or `scripts/`) for rewriting existing
  RAVEN-format YAMLs (e.g. the published `ecYeastGEM.yml`).
- Make `readDLKcatOutput`'s substrate-name match against
  `model.metNames` case-insensitive. `ismember` is case-sensitive by
  default, which can fail spuriously when an SBML loader normalizes
  case differently from how the DLKcat input was generated. geckopy
  uses case-insensitive matching here.
- Fix the `rxnsToClear(ecRxns) = false` block in `writeDLKcatInput`.
  The `rxnsToClear` array has length `numel(ecRxns)` (the count of
  selected reactions), but `ecRxns` (after `find()`) holds indices
  into the larger `model.ec.rxns` space. For any non-trivial
  selection, indices will exceed the array length and MATLAB will
  error out. The intended behaviour ("clear unselected reactions in
  the subset matrix") is already achieved by the preceding
  `reducedS(:, origRxnIdxs)` slice; the `rxnsToClear` block is
  dead/buggy and should be removed.

## get_enzyme_data subsystem

- **KEGG as an alternative protein-sequence source for
  `make_ec_model`.** Today `populate_enzyme_data` only consults
  UniProt. KEGG returns similar information (gene id, EC number,
  MW, sequence) and is often more complete for non-model
  organisms. The KEGG loader and downloader are already ported
  (in `geckopy.get_enzyme_data`), but `make_ec_model` does not
  use them yet. The work needed: add a KEGG fallback to stage 7
  of `make_ec_model`, controlled via adapter parameters.
