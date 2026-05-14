"""Sensitivity analysis on enzyme kcat values.

Mirrors GECKO MATLAB's `src/geckomat/kcat_sensitivity_analysis/`
directory.
"""
from .find_max_value import find_max_value

__all__ = ["find_max_value"]
