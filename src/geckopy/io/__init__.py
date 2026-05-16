"""SBML I/O for ecModels (GECKO 3 MW_KCAT encoding)."""
from .sbml import (
    read_sbml_ec_model,
    write_sbml_ec_model,
)

# Shorter aliases. ``read_sbml`` / ``write_sbml`` are the canonical
# public names; the longer ``..._ec_model`` versions are kept for
# backward compatibility and explicit-typing situations.
read_sbml = read_sbml_ec_model
write_sbml = write_sbml_ec_model

__all__ = [
    "read_sbml",
    "read_sbml_ec_model",
    "write_sbml",
    "write_sbml_ec_model",
]
