"""geckopy: Enzyme-constrained genome-scale metabolic modeling in Python."""
from geckopy.adapter import ModelAdapter, ModelParameters
from geckopy.ec_model import EcData, EcModel, Enzyme, make_ec_model
from geckopy.ec_model.pipeline.protein_pool import set_prot_pool_size
from geckopy.kcat_sensitivity_analysis import (
    fit_sigma,
    sensitivity_tuning,
    sigma_fitter,  # deprecated; alias of fit_sigma
)
from geckopy.limit_proteins import (
    apply_flux_data_constraints,
    calculate_f_factor,
    constrain_enz_concs,
    constrain_flux_data,  # deprecated; alias of apply_flux_data_constraints
    fill_enz_concs,
    flexibilize_enz_concs,
    relax_proteomics_greedy,
)
from geckopy.utilities import (
    ec_fseof,
    ec_fva,
    enzyme_usage,
    get_enzyme_bottlenecks,
    load_conventional_gem,
    load_ec_model,
    pfba_enzymes,
    save_ec_model,
)

# SBML I/O depends on libsbml at import time; cobra already makes libsbml a
# hard dep, so this should always succeed. The try/except is belt-and-braces
# so that an environment with a broken libsbml install does not break the
# top-level `import geckopy`.
try:
    from geckopy.io import read_sbml_ec_model, write_sbml_ec_model
except ImportError:  # pragma: no cover
    read_sbml_ec_model = None  # type: ignore[assignment]
    write_sbml_ec_model = None  # type: ignore[assignment]

__all__ = [
    "EcData",
    "EcModel",
    "Enzyme",
    "ModelAdapter",
    "ModelParameters",
    "apply_flux_data_constraints",
    "calculate_f_factor",
    "constrain_enz_concs",
    "constrain_flux_data",  # deprecated; alias of apply_flux_data_constraints
    "ec_fseof",
    "ec_fva",
    "enzyme_usage",
    "fill_enz_concs",
    "fit_sigma",
    "flexibilize_enz_concs",
    "get_enzyme_bottlenecks",
    "load_conventional_gem",
    "load_ec_model",
    "make_ec_model",
    "pfba_enzymes",
    "read_sbml_ec_model",
    "relax_proteomics_greedy",
    "save_ec_model",
    "sensitivity_tuning",
    "set_prot_pool_size",
    "sigma_fitter",  # deprecated; alias of fit_sigma
    "write_sbml_ec_model",
]
