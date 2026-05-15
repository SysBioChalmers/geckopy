"""General utility functions for ecModels.

Mirrors GECKO MATLAB's `src/geckomat/utilities/` directory.
"""
from .add_new_rxns_to_ec import (
    AddNewRxnsResult,
    NewEnzyme,
    add_new_rxns_to_ec,
)
from .ec_fseof import EcFseofResult, ec_fseof
from .ec_fva import ec_fva
from .enzyme_usage import EnzymeUsageResult, enzyme_usage
from .map_rxns_to_conv import MapRxnsResult, map_rxns_to_conv
from .report_enzyme_usage import EnzymeUsageReport, report_enzyme_usage

__all__ = [
    "AddNewRxnsResult",
    "EcFseofResult",
    "EnzymeUsageReport",
    "EnzymeUsageResult",
    "MapRxnsResult",
    "NewEnzyme",
    "add_new_rxns_to_ec",
    "ec_fseof",
    "ec_fva",
    "enzyme_usage",
    "map_rxns_to_conv",
    "report_enzyme_usage",
]
