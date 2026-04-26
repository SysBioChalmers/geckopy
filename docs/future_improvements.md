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
- Switch `getECstring` from accumulator (`EC_set` in/out) to a pure
  function returning just the formatted EC string for one raw input.
  The accumulator is a footgun: callers must remember to add a
  trailing space to `EC_set` before calling, otherwise tokens collide.
- Fix `getECstring` empty-input behaviour. It currently returns the
  bare string `"EC"` for empty input because `strsplit("", " ")`
  returns `{''}` and the loop runs once. Should return `""`.
- Validate tokens in `getECstring`. The current implementation will
  blindly produce `ECEC1.1.1.1` for already-prefixed input and
  `ECnotanec` for junk. Strip a leading `EC` (case-insensitive) before
  re-prefixing and warn-and-skip tokens that fail the EC-shape regex.
- Fix the validation regex in `getECfromGEM`. The current pattern
  `(\d\.(\w|-)+\.(\w|-)+\.(\w|-)+)(;\w+\.(\w|-)+\.(\w|-)+\.(\w|-)+)*(.*)`
  substituted with `$3` returns the level-3 character of the first EC
  for any valid input (e.g. `"1.2.3.4"` -> `"3"`), and the subsequent
  `~cellfun(@isempty, ...)` then flags every non-empty result as
  invalid. As written, every non-empty EC string is silently
  discarded. Replace with a straightforward
  `^TOKEN(;TOKEN)*$` validation where `TOKEN` is the canonical
  four-level dotted EC pattern with `-` allowed in any level.

## get_enzyme_data subsystem

- **KEGG as an alternative protein-sequence source for `make_ec_model`.**
  Currently `populate_enzyme_data` consults UniProt only. KEGG returns
  similar information (gene, EC, MW, sequence) and is sometimes more
  complete for non-model organisms. The KEGG loader and downloader are
  ported as part of `get_enzyme_data/`, but `make_ec_model` does not yet
  use them. Future work: add a KEGG fallback path in stage 7 of
  `make_ec_model`, controlled via adapter parameters.
