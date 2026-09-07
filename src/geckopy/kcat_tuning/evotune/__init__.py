"""Evotune: CMA-ES kcat tuning against experimental data.

See ``docs/evotune_kcat_tuning.md`` for the walkthrough:
:func:`~.tuning.screen_kcat_leverage`, then
:func:`~.tuning.select_tunable_mask`, then
:func:`~.tuning.cmaes_kcat_tuning`.

Modules in this package that need ``cma`` import it unconditionally --
it's a hard requirement here, gated one level up instead (see
``kcat_tuning/__init__.py``) so importing plain
``geckopy.kcat_tuning`` never requires the optional
``evotune`` extra. Install it with
``pip install geckopy[evotune]``.
"""
from .data import EvotuneData, load_evotune_data
from .distance import (
    BIOMASS_CARBON_EQUIV,
    INFEASIBLE_PENALTY,
    evotune_distance,
    compute_excarbon,
    dataset_rmse,
)
from .priors import (
    UNLABELLED_GROUP,
    build_sigma0_log,
    classify_kcat_source,
    classify_kcat_sources,
)
from .simulate import ConditionSimResult, simulate_evotune_dataset
from .tuning import (
    EvotuneResult,
    cmaes_kcat_tuning,
    screen_kcat_leverage,
    select_tunable_mask,
    tune_prior_penalty_weight,
)

__all__ = [
    "BIOMASS_CARBON_EQUIV",
    "INFEASIBLE_PENALTY",
    "UNLABELLED_GROUP",
    "EvotuneData",
    "EvotuneResult",
    "ConditionSimResult",
    "evotune_distance",
    "build_sigma0_log",
    "classify_kcat_source",
    "classify_kcat_sources",
    "cmaes_kcat_tuning",
    "compute_excarbon",
    "dataset_rmse",
    "load_evotune_data",
    "screen_kcat_leverage",
    "select_tunable_mask",
    "simulate_evotune_dataset",
    "tune_prior_penalty_weight",
]
