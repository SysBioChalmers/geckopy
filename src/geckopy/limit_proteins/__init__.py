"""Apply enzyme-concentration constraints to ec models.

Mirrors GECKO MATLAB's `src/geckomat/limit_proteins/` directory.
"""
from .calculate_f_factor import calculate_f_factor
from .constrain_enz_concs import constrain_enz_concs
from .constrain_flux_data import (
    apply_flux_data_constraints,
    constrain_flux_data,
)
from .fill_enz_concs import fill_enz_concs
from .flexibilize_enz_concs import FlexEnzResult, flexibilize_enz_concs
from .get_conc_control_coeffs import get_conc_control_coeffs
from .relax_proteomics_greedy import (
    GreedyRelaxResult,
    RelaxationStep,
    relax_proteomics_greedy,
)

__all__ = [
    "FlexEnzResult",
    "GreedyRelaxResult",
    "RelaxationStep",
    "apply_flux_data_constraints",
    "calculate_f_factor",
    "constrain_enz_concs",
    "constrain_flux_data",  # deprecated; alias of apply_flux_data_constraints
    "fill_enz_concs",
    "flexibilize_enz_concs",
    "get_conc_control_coeffs",
    "relax_proteomics_greedy",
]
