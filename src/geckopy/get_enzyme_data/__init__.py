"""Functions that build the enzyme-data side of an ec model.

Mirrors GECKO MATLAB's `src/geckomat/get_enzyme_data/` directory.
"""
from .ec_from_database import get_ec_from_database
from .ec_from_gem import get_ec_from_gem
from .find_ec_in_db import find_ec_in_db

__all__ = ["find_ec_in_db", "get_ec_from_database", "get_ec_from_gem"]
