"""Shared names for the EC layer (metabolites, reactions, subsystem).

Centralises the magic strings that several modules in geckopy use to
identify enzyme metabolites, usage reactions, the shared protein pool,
and the protein-usage subsystem tag. Promoted from private module-level
constants in ``protein_pool.py`` so that other modules
(``enzyme.py``, ``pfba_enzymes``, parallel ``ec_fva``,
``relax_proteomics_greedy``, SBML I/O, ...) can share one definition.
"""

PROT_PREFIX = "prot_"
USAGE_PREFIX = "usage_prot_"
POOL_ID = "prot_pool"
POOL_EXCHANGE_ID = "prot_pool_exchange"
PROTEIN_USAGE_SUBSYSTEM = "Protein usage"
