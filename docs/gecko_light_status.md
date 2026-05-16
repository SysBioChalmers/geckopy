# gecko-light status

GECKO supports two ways of building an ecModel: **full** and **light**.
This document explains the difference, the state of light support in
geckopy (short version: not implemented yet), and what would need to
be added.

## What full vs light means

When a metabolic reaction is catalysed by an enzyme, GECKO needs to
encode the protein cost of that catalysis. There are two ways:

- **Full ecModel.** For each enzyme, geckopy adds a dedicated
  pseudo-metabolite (`prot_<uniprot>`) and a dedicated usage reaction
  (`usage_prot_<uniprot>`). Reactions that are catalysed by multiple
  isozymes are split into one copy per isozyme. This gives you maximum
  bookkeeping flexibility: you can ask the solver "how much of enzyme
  P00350 is being used right now?" directly. The cost is a big LP — a
  yeast-GEM full ecModel has ~8000 reactions; a HumanGEM full one
  blows up to tens of thousands.
- **Light ecModel.** Skip the per-enzyme bookkeeping. There's only a
  single shared protein pool. Each catalysed reaction gets one extra
  stoichiometric coefficient `MW / (kcat * 3600)` against that pool.
  Isozymes are not split; instead the build picks the lowest-cost
  enzyme for each reaction. The LP stays roughly the same size as the
  starting GEM. The trade-off is that you can no longer answer
  per-enzyme questions in solver outputs.

Light is the right choice for genome-scale models too big for the
full layout to solve in reasonable time (HumanGEM, Recon3D). For
smaller models (yeast, E. coli), full is usually fine.

## TL;DR for geckopy

**Light is not implemented.** A live probe on the small test
fixture:

```python
make_ec_model(model, adapter, gecko_light=True)
# NotImplementedError: gecko_light mode is not yet implemented in geckopy.
```

Several downstream functions also refuse to run on a light model.
The data layer is partially in place (the `EcData.gecko_light`
flag exists and a few helpers respect it), but the main build
path (`make_ec_model`) and the kcat application path
(`apply_kcat_constraints`) need new code.

## How the data shapes differ

|  | Full | Light |
|---|---|---|
| Isozyme expansion | reactions split, one copy per isozyme with a `_EXP_<N>` suffix | reactions stay singular |
| `ec.rxns` | one entry per (reaction, isozyme) pair; ids match the cobra reactions exactly | duplicate entries per isozyme; ids carry a `###_` counter prefix |
| Per-enzyme pseudo-metabolite `prot_<id>` | yes | no |
| Per-enzyme usage reaction `usage_prot_<id>` | yes | no |
| Shared protein pool | yes | yes (single global cost applied directly to each catalysed reaction's stoichiometry) |

## What's already light-aware in geckopy

Some of the scaffolding was added during the port even though the
build path itself isn't implemented yet:

| File | What it does for light |
|---|---|
| `src/geckopy/ec_model/ec_data.py` | `gecko_light: bool` field, with a docstring describing the duplicate-rxn convention |
| `src/geckopy/ec_model/ec_model.py` | accepts `gecko_light=` kwarg and propagates it into the `EcData` |
| `src/geckopy/ec_model/pipeline/expand.py` | docstring notes "geckoLight skips this stage" |
| `src/geckopy/gather_kcats/fuzzy_kcat_matching.py` | strips the `###_` prefix when looking up cobra reactions for substrate / organism matching |
| `src/geckopy/gather_kcats/get_standard_kcat.py` | `_ec_rxn_to_cobra_id` helper strips the prefix; the subsystem-mean computation handles both forms |

## What raises `NotImplementedError` today

| Location | What it refuses to do |
|---|---|
| `src/geckopy/ec_model/make_ec_model.py:89` | build a light model at all |
| `src/geckopy/ec_model/pipeline/apply_kcat.py:84` | apply kcat constraints to a light model |
| `src/geckopy/ec_model/pipeline/fill_kcats.py:74` | run `get_kcat_across_isozymes` (not applicable to light; light doesn't split isozymes) |
| `src/geckopy/ec_model/enzyme.py:105, 248` | `Enzyme.concentration` and `Kcats.__setitem__` setters reject light models |
| `src/geckopy/utilities/pfba_enzymes.py:46` | `pfba_enzymes` needs usage reactions, which light models don't have |

The `Enzyme` proxy is full-only by design — it forwards to
`prot_<id>` / `usage_prot_<id>`, neither of which a light model has.
So most of the proxy is inapplicable, not just the setters.

## What would need to be implemented

To run the MATLAB `light_ecModel/protocol.m` tutorial end-to-end
in Python (steps 8-31 plus simulation), these code paths need new
logic:

1. **`make_ec_model(gecko_light=True)` body.** Skip `expand_model`;
   set up `ec.rxns` with duplicate entries (one per isozyme) and a
   `###_` counter prefix; skip `add_protein_pseudometabolites` and
   `add_protein_usage_reactions`. Still call
   `add_protein_pool_pseudometabolite` and
   `add_protein_pool_exchange_reaction` — light still uses a single
   shared pool.
2. **`apply_kcat_constraints` light branch.** Write
   `MW_sum / (kcat * 3600)` directly as a `prot_pool` coefficient on
   each catalysed cobra reaction's stoichiometry — no per-enzyme
   `prot_<id>` mets. For reactions with multiple isozymes, pick the
   lowest-cost one (this is what MATLAB's `selectKcatValue` does for
   light models).
3. **`selectKcatValue` for light.** Verify or extend the existing
   implementation to handle the light `ec.rxns` shape (duplicate
   entries per isozyme). Probably needs an aggregation step that
   collapses per-isozyme kcats to one per cobra reaction.
4. **Tutorial port.** A `tutorials/light_ecModel/protocol.py`
   mirroring `tutorials/full_ecModel/protocol.py` but using the
   light path. The MATLAB version uses a trimmed HumanGEM; we
   could use the same, or build a light yeast-GEM as a smaller
   smoke test.
5. **`Enzyme` proxy fallback.** Relax the light-mode guards on
   read-only properties (`concentration`, `mw`, `gene`, `sequence`)
   so users can still query enzyme metadata on a light model.
   Setters stay blocked since there's no per-enzyme constraint to
   update.

## Estimated scope

Comparable to one of the larger original ports — roughly the size
of Ports 1, 5, and 8 from the recent batch combined.

- `make_ec_model` light branch: ~150 LOC + ~10 tests
- `apply_kcat_constraints` light branch: ~80 LOC + ~5 tests
- `selectKcatValue` light-aware aggregation: ~30 LOC + ~3 tests
- `Enzyme` proxy light-mode relaxation: ~20 LOC + ~3 tests
- Tutorial port + smoke test: ~200 LOC + a few hours of running

Total: ~1-2 days of focused work.

## Recommendation

Not blocking for current users — the full layout works end-to-end
including the full yeast-GEM tutorial. Worth implementing when:

- a concrete user needs a HumanGEM (or similarly-large) ecModel
  built through geckopy, or
- closing the MATLAB parity gap becomes the priority.

If you do start: **begin with the `make_ec_model` light branch**.
It forces the data-layer shape to be correct, and the downstream
functions can be filled in next as their tests demand.
