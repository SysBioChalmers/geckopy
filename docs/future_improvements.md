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
