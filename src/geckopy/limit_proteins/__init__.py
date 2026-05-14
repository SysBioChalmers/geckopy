"""Apply enzyme-concentration constraints to ec models.

Mirrors GECKO MATLAB's `src/geckomat/limit_proteins/` directory.
"""
from .calculate_f_factor import calculate_f_factor

__all__ = ["calculate_f_factor"]
