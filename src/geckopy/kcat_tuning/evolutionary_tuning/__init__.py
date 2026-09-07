"""Evolutionary (CMA-ES) kcat tuning: fitting kcats to experimental data.

See ``docs/evolutionary_kcat_tuning.md`` for the walkthrough:
:func:`~.tuning.screen_kcat_leverage`, then
:func:`~.tuning.select_tunable_mask`, then
:func:`~.tuning.cmaes_kcat_tuning`.

Modules in this package that need ``cma`` import it unconditionally --
it's a hard requirement here, gated one level up instead (see
``kcat_tuning/__init__.py``) so importing plain
``geckopy.kcat_tuning`` never requires the optional
``evolutionary-tuning`` extra. Install it with
``pip install geckopy[evolutionary-tuning]``.
"""
from .data import TuningData, load_tuning_data
from .distance import (
    BIOMASS_CARBON_EQUIV,
    INFEASIBLE_PENALTY,
    tuning_distance,
    compute_excarbon,
    dataset_rmse,
)
from .priors import (
    UNLABELLED_GROUP,
    build_sigma0_log,
    classify_kcat_source,
    classify_kcat_sources,
)
from .simulate import ConditionSimResult, simulate_tuning_dataset
from .tuning import (
    EvolutionaryTuningResult,
    cmaes_kcat_tuning,
    screen_kcat_leverage,
    select_tunable_mask,
    tune_prior_penalty_weight,
)

__all__ = [
    "BIOMASS_CARBON_EQUIV",
    "INFEASIBLE_PENALTY",
    "UNLABELLED_GROUP",
    "TuningData",
    "EvolutionaryTuningResult",
    "ConditionSimResult",
    "tuning_distance",
    "build_sigma0_log",
    "classify_kcat_source",
    "classify_kcat_sources",
    "cmaes_kcat_tuning",
    "compute_excarbon",
    "dataset_rmse",
    "load_tuning_data",
    "screen_kcat_leverage",
    "select_tunable_mask",
    "simulate_tuning_dataset",
    "tune_prior_penalty_weight",
]
