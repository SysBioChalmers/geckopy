"""Bayesian (ABC-SMC) kcat tuning.

Modules in this package that need ``pyabc`` import it unconditionally --
it's a hard requirement here, gated one level up instead (see
``kcat_sensitivity_analysis/__init__.py``) so importing plain
``geckopy.kcat_sensitivity_analysis`` never requires the optional
``pyabc`` dependency. Install it with ``pip install geckopy[bayesian]``.
"""
from .data import BayesianData, load_bayesian_data
from .distance import (
    BIOMASS_CARBON_EQUIV,
    INFEASIBLE_PENALTY,
    bayesian_distance,
    compute_excarbon,
    dataset_rmse,
)
from .simulate import ConditionSimResult, simulate_bayesian_dataset

__all__ = [
    "BIOMASS_CARBON_EQUIV",
    "INFEASIBLE_PENALTY",
    "BayesianData",
    "ConditionSimResult",
    "bayesian_distance",
    "compute_excarbon",
    "dataset_rmse",
    "load_bayesian_data",
    "simulate_bayesian_dataset",
]
