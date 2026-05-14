"""Apply enzyme-concentration constraints to ec models.

Mirrors GECKO MATLAB's `src/geckomat/limit_proteins/` directory.
"""
from .calculate_f_factor import calculate_f_factor
from .fill_enz_concs import fill_enz_concs

__all__ = ["calculate_f_factor", "fill_enz_concs"]
