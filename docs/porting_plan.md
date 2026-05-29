# GECKO to geckopy porting plan

## Target version
GECKO v3.2.5 (pinned).

## Scope
Full port of all 82 MATLAB functions in src/geckomat/ plus tutorials.

## Fidelity goal
Model structure and I/O must be comparable to the MATLAB output.
Numerical differences are acceptable where they arise from solver
or library choice.

## Data model
EcModel subclasses cobra.Model and carries an `.ec` attribute of
type EcData (dataclass) holding kcat, enzymes, rxnEnzMat, etc.

## Dependency mapping
- RAVEN: mostly replaced by cobrapy; residual functions go into
  geckopy.raven_shim, to be promoted to a separate RAVENpy package
  once its scope is known.
- COBRA: replaced by cobrapy.
- Gurobi / solvers: via cobrapy / optlang.
- BRENDA / UniProt / KEGG HTTP APIs: via `requests`.
- MATLAB .mat reference outputs: via `scipy.io.loadmat`.

## Function inventory

| Subsystem | MATLAB function | Python target | Deps (RAVEN / other) | Status |
|-----------|-----------------|---------------|----------------------|--------|
| model_adapter | ... | ... | ... | todo |