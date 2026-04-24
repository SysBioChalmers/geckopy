"""Pipeline stages of make_ec_model, ported from GECKO MATLAB."""
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

__all__ = [
    "allocate_ec_for_catalyzed_reactions",
    "build_rxn_enzyme_coupling",
    "convert_to_irreversible",
    "expand_model",
    "invert_backwards_only_reactions",
    "populate_enzyme_data",
    "remove_pseudoreaction_gprs",
]
