# RAVEN function inventory

GECKO MATLAB depends heavily on
[RAVEN](https://github.com/SysBioChalmers/RAVEN) for model
manipulation, FBA solves, and SBML / YAML I/O. As the geckopy port
progressed, every RAVEN call in the MATLAB source either resolved
to a direct cobra-py equivalent, was re-implemented inside
geckopy, or was identified as a gap with no Python counterpart.

This document catalogues those calls so future contributors can:

1. Quickly find the cobra-py / geckopy substitute for a RAVEN
   function when porting more MATLAB scripts.
2. Scope a future **ravenpy** package: a Python port of the RAVEN
   functions that have no cobra-py equivalent today.

The list is based on a scan of `src/geckomat/**/*.m` in GECKO
MATLAB and the existing pipeline modules in
`src/geckopy/ec_model/pipeline/`. It is not exhaustive across all
of RAVEN (RAVEN has ~100 public functions) -- only the ones GECKO
actually exercises.

## 1. Replaced by cobra-py directly

| RAVEN | cobra-py equivalent |
|---|---|
| `setParam(model, 'lb'\|'ub'\|'eq', rxn, val)` | `model.reactions.get_by_id(rxn).lower_bound = val` (or `.upper_bound`; for `obj`, `model.objective = rxn`) |
| `solveLP(model)` | `model.optimize()` |
| `importModel(filename)` | `cobra.io.read_sbml_model(filename)` |
| `exportModel(model, filename)` | `cobra.io.write_sbml_model(model, filename)` |
| `getIndexes(model, ids, 'rxns'\|'mets'\|'genes')` | `model.reactions.get_by_id(id)` (cobra-py is id-keyed) |
| `addRxns(model, rxnsToAdd)` | `model.add_reactions([cobra.Reaction(...)])` |
| `addMets(model, metsToAdd)` | `model.add_metabolites([cobra.Metabolite(...)])` |
| `removeReactions(model, rxnList)` | `model.remove_reactions(rxnList)` |
| `removeMets(model, metList)` | `model.remove_metabolites(metList)` |
| `removeGenes(model, geneList)` | `model.genes.remove(model.genes.get_by_id(g))` per id |
| `addGenes(model, geneList)` | implicitly added via reaction GPR; cobra-py does not need an explicit call |

## 2. Replaced by geckopy internals

| RAVEN | geckopy module |
|---|---|
| `expandModel` | `src/geckopy/ec_model/pipeline/expand.py::expand_model` (uses cobra-py's GPR AST instead of string manipulation) |
| `convertToIrrev` | `src/geckopy/ec_model/pipeline/preprocess.py::convert_to_irreversible` (forward + reverse pairs with `_REV` suffix) |
| pseudoreaction handling / backward-only inversion / bound validation | also `preprocess.py` (Stages 1-4 of the build pipeline) |
| `readYAMLmodel`, `writeYAMLmodel` (cobra portion) | `src/geckopy/utilities/load_ec_model.py`, `save_ec_model.py` -- but see notes in [yaml_format.md](yaml_format.md) about the canonical schema |

## 3. RAVEN-only, no cobra-py equivalent — ravenpy candidates

These are the functions worth porting to a separate **ravenpy**
package. Each provides a biochemically-informed structural or
optimization utility that is orthogonal to cobra-py's core FBA
engine.

| Function | What it does | Why cobra-py can't do it natively |
|---|---|---|
| `fillGaps` | Adds reactions from a reference model to repair an infeasible target | cobra-py has `cobra.flux_analysis.gapfill` but it differs in algorithm and signature; RAVEN's variant is what GECKO callers expect |
| `simplifyModel` | Removes dead-end metabolites and orphan reactions via iterative constraint propagation | cobra-py has no built-in equivalent; ad-hoc per-tool implementations exist |
| `contractModel` | Removes reactions via flux-based / structural analysis (model reduction) | not in cobra-py |
| `haveFlux(model, rxnList)` | Boolean per-rxn: can carry non-zero flux? (numerically tolerant) | `cobra.flux_analysis.find_blocked_reactions` exists but signature and tolerance handling differ |
| `getAllowedBounds(model)` | Tight per-reaction flux bounds via FVA under current constraints | `cobra.flux_analysis.flux_variability_analysis` exists but returns a DataFrame; not a drop-in replacement |
| `gprConsolidate(model, grRules)` | Simplifies GPRs via Boolean algebra (`(A or A) and B` -> `A and B`) | cobra-py has no GPR simplifier; `cobra.core.gene.GPR` only parses and evaluates |
| `standardizeGrRules(model)` | Normalises GPR syntax (parenthesisation, whitespace) | cobra-py's GPR parser accepts a wider syntax but does not normalise on the way out |
| `changeGeneAssociation(model, rxn, newGpr)` | Replaces a reaction's GPR while keeping the gene list consistent | partly handled by `rxn.gene_reaction_rule = ...`, but cobra-py does not garbage-collect orphan genes |
| `getGenesFromGrRules(grRules)` | Extracts gene IDs from a GPR rule string | cobra-py's `GPR.genes` works on the parsed AST, not raw strings |
| `mergeCompartments` | Consolidates two compartments into one | not in cobra-py |
| `copyToCompartments(model, mets, src, tgt)` | Replicates metabolites and reactions into another compartment | not in cobra-py |
| `removeBadRxns` | Flags reactions with structural issues for review/removal | not in cobra-py |
| `readYAMLmodel`, `writeYAMLmodel` (RAVEN dialect) | Reads/writes the RAVEN-flavoured YAML (different schema from cobra-py's) | cobra-py's YAML loader expects its own schema; reading the legacy `ecYeastGEM.yml` format needs a translator |

## 4. MATLAB-specific helpers — not needed in Python

| Function | Reason it's not ported |
|---|---|
| `solveLP(..., minFlux, hsSol)` hot-start | MATLAB-specific solver state; optlang/cobra-py handle this internally |
| `setRavenSolver(...)` | MATLAB CPLEX / Gurobi setup; cobra-py uses `model.solver = "gurobi"` |
| `importExcelModel`, `checkFileExistence`, dialog UI prompts | File I/O glue, not algorithmic |
| `convertCharArray`, `dispEM` | MATLAB cell/struct utilities |

## Suggested ravenpy scope

For a v0.1 of ravenpy, prioritise the GPR utilities
(`gprConsolidate`, `standardizeGrRules`, `getGenesFromGrRules`) +
`fillGaps` + `simplifyModel`. These are the most commonly called
in MATLAB GECKO and have the clearest "RAVEN does this, cobra-py
doesn't" signature. The compartment-manipulation functions
(`mergeCompartments`, `copyToCompartments`) and `removeBadRxns`
are next-most-valuable but used in fewer places.

The YAML reader/writer is its own concern -- the geckopy port
already chose a canonical-cobra-py-compatible schema (see
[yaml_format.md](yaml_format.md)), so the RAVEN-dialect
reader is mostly relevant for backward compatibility with legacy
files.
