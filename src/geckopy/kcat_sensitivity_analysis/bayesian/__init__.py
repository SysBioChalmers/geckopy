"""Bayesian kcat tuning: fitting kcats to experimental data.

:func:`~.cmaes_tuning.cmaes_kcat_tuning` (CMA-ES) is the recommended
entry point; see its module docstring for why.
:func:`~.tuning.bayesian_kcat_tuning` (ABC-SMC) is the earlier, still
supported path.

Modules in this package that need ``pyabc`` or ``cma`` import them
unconditionally -- they're hard requirements here, gated one level up
instead (see ``kcat_sensitivity_analysis/__init__.py``) so importing
plain ``geckopy.kcat_sensitivity_analysis`` never requires the optional
``bayesian`` dependencies. Install them with
``pip install geckopy[bayesian]``.
"""
from .cmaes_tuning import (
    cmaes_kcat_tuning,
    screen_kcat_leverage,
    select_tunable_mask,
)
from .data import BayesianData, load_bayesian_data
from .distance import (
    BIOMASS_CARBON_EQUIV,
    INFEASIBLE_PENALTY,
    bayesian_distance,
    compute_excarbon,
    dataset_rmse,
)
from .diagnostics import GenerationDiagnostics, GroupDiagnostics, compute_generation_diagnostics
from .priors import (
    UNLABELLED_GROUP,
    build_kcat_prior,
    build_sigma0_log,
    classify_kcat_source,
    classify_kcat_sources,
)
from .selection import (
    SelectionResult,
    next_quantile_epsilon,
    quantile_epsilon_select,
    truncation_select,
)
from .simulate import ConditionSimResult, simulate_bayesian_dataset
from .transition import GeckoTransition
from .tuning import BayesianTuningResult, bayesian_kcat_tuning

__all__ = [
    "BIOMASS_CARBON_EQUIV",
    "INFEASIBLE_PENALTY",
    "UNLABELLED_GROUP",
    "BayesianData",
    "BayesianTuningResult",
    "ConditionSimResult",
    "GeckoTransition",
    "GenerationDiagnostics",
    "GroupDiagnostics",
    "SelectionResult",
    "bayesian_distance",
    "bayesian_kcat_tuning",
    "build_kcat_prior",
    "build_sigma0_log",
    "classify_kcat_source",
    "classify_kcat_sources",
    "cmaes_kcat_tuning",
    "compute_excarbon",
    "compute_generation_diagnostics",
    "dataset_rmse",
    "load_bayesian_data",
    "next_quantile_epsilon",
    "quantile_epsilon_select",
    "screen_kcat_leverage",
    "select_tunable_mask",
    "simulate_bayesian_dataset",
    "truncation_select",
]
