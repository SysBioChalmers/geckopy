"""SBML I/O for ecModels (GECKO 3 MW_KCAT encoding)."""
from .sbml import read_sbml_ec_model, write_sbml_ec_model

__all__ = ["read_sbml_ec_model", "write_sbml_ec_model"]
