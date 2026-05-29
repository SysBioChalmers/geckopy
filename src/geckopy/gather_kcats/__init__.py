"""Gather kcat values for ec model reactions from various sources.

Mirrors GECKO MATLAB's `src/geckomat/gather_kcats/` directory.
"""
from .fuzzy_kcat_matching import fuzzy_kcat_matching
from .get_standard_kcat import assign_standard_kcat, get_standard_kcat
from .merge_dlkcat_and_fuzzy_kcats import merge_dlkcat_and_fuzzy_kcats
from .read_dlkcat_output import read_dlkcat_output
from .remove_standard_kcat import remove_standard_kcat
from .run_dlkcat import run_dlkcat
from .select_kcat_value import apply_kcat_list, select_kcat_value
from .write_dlkcat_input import write_dlkcat_input

__all__ = [
    "apply_kcat_list",
    "assign_standard_kcat",
    "fuzzy_kcat_matching",
    "get_standard_kcat",  # deprecated; alias of assign_standard_kcat
    "merge_dlkcat_and_fuzzy_kcats",
    "read_dlkcat_output",
    "remove_standard_kcat",
    "run_dlkcat",
    "select_kcat_value",  # deprecated; alias of apply_kcat_list
    "write_dlkcat_input",
]
