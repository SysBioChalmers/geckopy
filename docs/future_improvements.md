# Future improvements

This file tracks ideas for improvements to geckopy and to MATLAB GECKO
that are deferred until after the initial port is complete.

## Conventions / data model

- **Missing-kcat sentinel.** geckopy uses `NaN` in `ec.kcat` to mark
  missing values; MATLAB GECKO uses `0`. The model I/O layer must
  translate at the boundary (`0` -> `NaN` on read, `NaN` -> `0` on write).
  A simpler long-term fix would be to switch MATLAB GECKO to `NaN` as
  well, eliminating the translation layer.

- **Custom subunit stoichiometry file format.** The `stoicho` column in
  `customKcats.tsv` is currently dead code in MATLAB and ignored in
  geckopy. Two design options for proper support: a separate
  `customComplexes.tsv` file, or extending the ComplexPortal JSON
  format with a "user_overrides" section. Decision deferred until
  after the basic port is complete.

## API simplification

- **Mode B in `apply_custom_kcats`.** The "proteins only, no rxns"
  mode in `applyCustomKcats.m` has unclear semantics (the docstring
  and code disagree about whether to apply the full-match gate).
  Geckopy follows the implementation; MATLAB should clean up the
  docstring or drop the mode.

- **Function naming.** `get_kcat_across_isozymes` is a literal port of
  the MATLAB name but the verb "get" is misleading; the function
  modifies `ec.kcat` in place. A cleaner name would be
  `fill_kcats_from_isozymes`. Deferred to a coordinated rename across
  both implementations.

## MATLAB GECKO changes

A list of all places where MATLAB GECKO behavior should be brought into
line with geckopy. Each is also tagged with a `MATLAB-COMPAT:` comment
in the geckopy source code. Run `grep -rn "MATLAB-COMPAT:" src/` for
the full list.

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

- **KEGG as an alternative protein-sequence source for `make_ec_model`.**
  Currently `populate_enzyme_data` consults UniProt only. KEGG returns
  similar information (gene, EC, MW, sequence) and is sometimes more
  complete for non-model organisms. The KEGG loader and downloader are
  ported as part of `get_enzyme_data/`, but `make_ec_model` does not yet
  use them. Future work: add a KEGG fallback path in stage 7 of
  `make_ec_model`, controlled via adapter parameters.
