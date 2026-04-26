"""Functions that build the enzyme-data side of an ec model.

Mirrors GECKO MATLAB's `src/geckomat/get_enzyme_data/` directory.
"""
from .ec_from_gem import get_ec_from_gem
from .ec_string import get_ec_string

__all__ = ["get_ec_from_gem", "get_ec_string"]
