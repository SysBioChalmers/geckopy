"""Bayesian kcat tuning: fitting kcats to experimental data with CMA-ES.

See ``docs/bayesian_kcat_tuning.md`` for the walkthrough:
:func:`~.tuning.screen_kcat_leverage`, then
:func:`~.tuning.select_tunable_mask`, then
:func:`~.tuning.cmaes_kcat_tuning`.

Modules in this package that need ``cma`` import it unconditionally --
it's a hard requirement here, gated one level up instead (see
``kcat_sensitivity_analysis/__init__.py``) so importing plain
``geckopy.kcat_sensitivity_analysis`` never requires the optional
``bayesian`` extra. Install it with ``pip install geckopy[bayesian]``.
"""
from .data import BayesianData, load_bayesian_data
from .distance import (
    BIOMASS_CARBON_EQUIV,
    INFEASIBLE_PENALTY,
    bayesian_distance,
    compute_excarbon,
    dataset_rmse,
)
from .priors import (
    UNLABELLED_GROUP,
    build_sigma0_log,
    classify_kcat_source,
    classify_kcat_sources,
)
from .simulate import ConditionSimResult, simulate_bayesian_dataset
from .tuning import (
    BayesianTuningResult,
    cmaes_kcat_tuning,
    screen_kcat_leverage,
    select_tunable_mask,
    tune_prior_penalty_weight,
)

__all__ = [
    "BIOMASS_CARBON_EQUIV",
    "INFEASIBLE_PENALTY",
    "UNLABELLED_GROUP",
    "BayesianData",
    "BayesianTuningResult",
    "ConditionSimResult",
    "bayesian_distance",
    "build_sigma0_log",
    "classify_kcat_source",
    "classify_kcat_sources",
    "cmaes_kcat_tuning",
    "compute_excarbon",
    "dataset_rmse",
    "load_bayesian_data",
    "screen_kcat_leverage",
    "select_tunable_mask",
    "simulate_bayesian_dataset",
    "tune_prior_penalty_weight",
]
