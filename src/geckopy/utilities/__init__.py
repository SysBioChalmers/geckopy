"""General utility functions for ecModels.

Mirrors GECKO MATLAB's `src/geckomat/utilities/` directory.
"""
from .enzyme_usage import EnzymeUsageResult, enzyme_usage
from .report_enzyme_usage import EnzymeUsageReport, report_enzyme_usage

__all__ = [
    "EnzymeUsageReport",
    "EnzymeUsageResult",
    "enzyme_usage",
    "report_enzyme_usage",
]
