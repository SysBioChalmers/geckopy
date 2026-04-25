"""Pipeline stages of make_ec_model, ported from GECKO MATLAB."""
from .apply_kcat import apply_kcat_constraints
from .expand import expand_model
from .populate_ec import (
    allocate_ec_for_catalyzed_reactions,
    build_rxn_enzyme_coupling,
    populate_enzyme_data,
)
from .preprocess import (
    convert_to_irreversible,
    invert_backwards_only_reactions,
    remove_pseudoreaction_gprs,
)
from .protein_pool import (
    add_protein_pool_exchange_reaction,
    add_protein_pool_pseudometabolite,
    add_protein_pseudometabolites,
    add_protein_usage_reactions,
    set_prot_pool_size,
)
from .query import get_reactions_from_enzyme
from .set_kcat import set_kcat_for_reactions

__all__ = [
    "add_protein_pool_exchange_reaction",
    "add_protein_pool_pseudometabolite",
    "add_protein_pseudometabolites",
    "add_protein_usage_reactions",
    "allocate_ec_for_catalyzed_reactions",
    "apply_kcat_constraints",
    "build_rxn_enzyme_coupling",
    "convert_to_irreversible",
    "expand_model",
    "get_reactions_from_enzyme",
    "invert_backwards_only_reactions",
    "populate_enzyme_data",
    "remove_pseudoreaction_gprs",
    "set_kcat_for_reactions",
    "set_prot_pool_size",
]
