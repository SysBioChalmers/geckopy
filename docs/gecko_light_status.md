# gecko-light status

Assessment of how much of GECKO MATLAB's `light_ecModel` protocol
is implementable in geckopy today, and what would need to be added
to support it end-to-end.

## TL;DR

**gecko-light is not implemented**. A live probe on `examples/ecTestGEM`:

```python
make_ec_model(model, adapter, gecko_light=True)
# NotImplementedError: gecko_light mode is not yet implemented in geckopy.
```

Several downstream functions also raise on light models. The
`EcData.gecko_light` flag exists and a handful of helpers already
respect it, so the scaffolding is partially in place -- but the
core build pipeline (`make_ec_model`) and the kcat application path
(`apply_kcat_constraints`) need new code paths.

## Conceptual difference from the full formulation

| | Full | Light |
|---|---|---|
| Isozyme expansion | reactions split per isozyme with `_EXP_<N>` suffix | reactions stay singular |
| `ec.rxns` | one entry per (rxn, isozyme) pair, ids match cobra rxns | duplicate entries per isozyme, ids carry a `###_` counter prefix |
| Per-enzyme pseudometabolite `prot_<id>` | yes | no |
| Per-enzyme usage reaction `usage_prot_<id>` | yes | no |
| Protein pool | yes | yes (single global cost coefficient applied directly to each catalysed cobra reaction's S-matrix) |

The light formulation is meant for very large models (HumanGEM,
Recon3D) where the full per-enzyme expansion creates too many
reactions / metabolites and the LP becomes slow.

## What's already light-aware in geckopy

| File | What it does for light |
|---|---|
| `src/geckopy/ec_model/ec_data.py` | `gecko_light: bool` field with docstring describing the duplicate-rxns convention |
| `src/geckopy/ec_model/ec_model.py` | accepts `gecko_light=` kwarg, propagates to `EcData` |
| `src/geckopy/ec_model/pipeline/expand.py` | module docstring notes "geckoLight skips this stage" |
| `src/geckopy/gather_kcats/fuzzy_kcat_matching.py` | strips the `###_` prefix when resolving cobra reactions for substrate / organism matching |
| `src/geckopy/gather_kcats/get_standard_kcat.py` | `_ec_rxn_to_cobra_id` helper strips the prefix; subsystem-mean computation handles both forms |

## What raises `NotImplementedError` today

| Location | Behaviour |
|---|---|
| `src/geckopy/ec_model/make_ec_model.py:89` | refuses to build a light model at all |
| `src/geckopy/ec_model/pipeline/apply_kcat.py:84` | refuses to apply constraints to a light model |
| `src/geckopy/ec_model/pipeline/fill_kcats.py:74` | `get_kcat_across_isozymes` is not applicable to light (no isozyme split) |
| `src/geckopy/ec_model/enzyme.py:105, 248` | `Enzyme.concentration` and `Kcats.__setitem__` setters reject light models |
| `src/geckopy/utilities/pfba_enzymes.py:46` | `pfba_enzymes` needs usage reactions, which light models don't have |

The `Enzyme` proxy itself is full-only by design (it forwards to
`prot_<id>` / `usage_prot_<id>`); the light formulation has neither,
so most of the proxy is inapplicable.

## What would need to be implemented

To match the MATLAB `light_ecModel/protocol.m` end-to-end (STEPs
8-31, plus simulation), the following code paths need new logic.

1. **`make_ec_model(gecko_light=True)`** — skip `expand_model`,
   set up `ec.rxns` with duplicate entries (one per isozyme) and a
   `###_` counter prefix, allocate ec data accordingly. Skip
   `add_protein_pseudometabolites` and `add_protein_usage_reactions`.
   Still call `add_protein_pool_pseudometabolite` and
   `add_protein_pool_exchange_reaction` (light uses a single pool).
2. **`apply_kcat_constraints` light branch** — write
   `MW_sum / (kcat * 3600)` directly as a `prot_pool` coefficient on
   each catalysed cobra reaction's S-matrix (no per-enzyme prot mets).
   Handle multi-isozyme reactions by picking the lowest-cost isozyme
   (MATLAB's `selectKcatValue` for light does this aggregation).
3. **`selectKcatValue` parity** — verify or extend the existing
   implementation handles the light ec.rxns shape (duplicate entries
   per isozyme) correctly. Likely needs an aggregation step that
   picks one kcat per cobra reaction from the per-isozyme entries.
4. **Tutorial fixture** — a `tutorials/light_ecModel/protocol.py`
   mirroring `tutorials/full_ecModel/protocol.py` but using the
   light path and a smaller model (likely a trimmed HumanGEM or the
   existing yeast-GEM with `gecko_light=True`).
5. **`Enzyme` proxy fallback for light** — relax the light-mode
   guards on read-only properties (concentration getter, mw getter,
   gene getter) so users can still query enzyme metadata on a light
   model. Setters stay blocked since there's no per-enzyme constraint
   to update.

## Estimated scope

The work is comparable to one of the original `make_ec_model`
pipeline stages: roughly the size of Ports 1, 5, and 8 combined.

- `make_ec_model` light branch: ~150 LOC + tests (~10 cases)
- `apply_kcat_constraints` light branch: ~80 LOC + tests (~5 cases)
- `selectKcatValue` light-aware aggregation: ~30 LOC + tests (~3 cases)
- `Enzyme` proxy light-mode relaxation: ~20 LOC + tests (~3 cases)
- Tutorial port + smoke test: ~200 LOC + a few hours of running

Total: ~1-2 days of focused work.

## Recommendation

Not blocking for current users (the full model works end-to-end
including the yeast-GEM tutorial). Implement when there is a
concrete user pulling for a HumanGEM ecModel via geckopy, or when
the maintainer wants to close the MATLAB parity gap.

If implementing, **start with `make_ec_model`'s light branch**: it
forces the data-layer shape to be correct, and the downstream
functions can be filled in next as their tests are added.
