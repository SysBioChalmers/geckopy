"""Real-scale smoke test for Bayesian kcat tuning.

Runs all 4 selection x regularization combinations against the real
``ecYeastGEM`` tutorial model and the real GECKO
``bayesianFluxData.tsv``/``bayesianMaxGrowth.tsv``/``bayesianZeroExch.tsv``
tutorial data (ported verbatim from the GECKO MATLAB repo's
``tutorials/full_ecModel/data/``), per the plan's Test Strategy #4.

Mirrors ``tests/test_light_humangem_smoke.py``'s style: skipped
gracefully when the (gitignored, build-output) ecModel YAML isn't
present -- run ``tutorials/full_ecModel/protocol.py`` first, or point
``ECYEASTGEM_YML`` at an existing one.

This is intentionally a *small-schedule* smoke test, not a rigorous
convergence run: ``schedule_samples``/``max_generations`` are kept
small so the four combinations complete in a few minutes total rather
than the hours a real tuning run (MATLAB's default schedule reaches
1000 samples/generation) would take against ~4800 real tunable kcats x
41 real experimental conditions. The point here is "does the full
real-scale pipeline actually run and produce sane output," and a
first head-to-head timing/RMSE/diagnostic comparison -- not a
converged posterior.

geckopy has no generic ``make_anaerobic``/``change_protein_biomass``
adapter hooks yet (see ``simulate.py``'s docstring), so rows requiring
those adjustments run without them; combined with parameters this
port doesn't carry over (MATLAB's low-rank kernel, exploit/explore
mixture), RMSE is expected to be in a rough ballpark, not a match to
MATLAB's reported 0.87-0.90 benchmark (an explicit non-goal).

Run explicitly with:

    pytest tests/test_bayesian_tuning_smoke.py -m smoke -q -s
"""
from __future__ import annotations

import os
import time
from pathlib import Path

import numpy as np
import pytest

from geckopy import ModelAdapter
from geckopy.utilities.load_ec_model import load_ec_model
from geckopy.adapter.params import BayesianParams
from geckopy.kcat_sensitivity_analysis.bayesian.data import load_bayesian_data
from geckopy.kcat_sensitivity_analysis.bayesian.tuning import bayesian_kcat_tuning

_TUTORIAL_DIR = Path(__file__).resolve().parent.parent / "tutorials" / "full_ecModel"
_ECYEASTGEM_YML = Path(
    os.environ.get("ECYEASTGEM_YML", _TUTORIAL_DIR / "models" / "ecYeastGEM.yml")
)

pytestmark = [
    pytest.mark.smoke,
    pytest.mark.skipif(
        not _ECYEASTGEM_YML.is_file(),
        reason=(
            f"ecYeastGEM.yml not found at {_ECYEASTGEM_YML} (a gitignored "
            "tutorial build output). Run tutorials/full_ecModel/protocol.py "
            "first, or set ECYEASTGEM_YML to point at an existing one."
        ),
    ),
]

# Deliberately small -- see module docstring.
_SMOKE_PARAMS = BayesianParams(
    schedule_generations=[1],
    schedule_samples=[8],
    min_keep=0.3,
    max_keep=0.6,
    rmse_threshold=-1.0,  # unreachable -> always runs exactly max_generations
    max_generations=2,
)

_COMBINATIONS = [
    ("truncation", "shrinkage"),
    ("truncation", "importance_weighting"),
    ("quantile_epsilon", "shrinkage"),
    ("quantile_epsilon", "importance_weighting"),
]


def _load_model_and_data():
    adapter = ModelAdapter.from_folder(_TUTORIAL_DIR)
    model = load_ec_model("ecYeastGEM.yml", adapter=adapter)
    bay_data = load_bayesian_data(adapter)
    assert bay_data.flux_data is not None and bay_data.max_grate is not None, (
        "Expected both bayesianFluxData.tsv and bayesianMaxGrowth.tsv to be "
        "present for this smoke test."
    )
    return model, adapter, bay_data


@pytest.mark.parametrize("selection,regularization", _COMBINATIONS)
def test_bayesian_tuning_runs_at_real_scale(selection, regularization):
    model, adapter, bay_data = _load_model_and_data()
    original_kcat = model.ec.kcat.copy()
    n_tunable = int((original_kcat > 0).sum())
    assert n_tunable > 1000, (
        f"Expected thousands of tunable kcats on the real model; got {n_tunable}."
    )

    t0 = time.perf_counter()
    result = bayesian_kcat_tuning(
        model, adapter=adapter, params=_SMOKE_PARAMS, bay_data=bay_data,
        selection=selection, regularization=regularization,
        seed=0, verbose=False,
    )
    wall_time = time.perf_counter() - t0

    assert result.n_generations == 2
    assert len(result.rmse_trace) == 2
    assert all(np.isfinite(r) for r in result.rmse_trace)
    assert len(result.diagnostics_trace) == 2

    last = result.diagnostics_trace[-1]
    assert set(last.by_group) >= {"brenda", "dlkcat", "custom"}
    for group_diag in last.by_group.values():
        assert 0.0 <= group_diag.frac_active <= 1.0
        assert 0.0 <= group_diag.frac_near_prior <= 1.0

    print(
        f"\n[smoke:{selection}/{regularization}] "
        f"wall_time={wall_time:.1f}s "
        f"rmse_trace={[f'{r:.3g}' for r in result.rmse_trace]} "
        f"n_tunable={n_tunable} "
        + " | ".join(
            f"{name}: active={g.frac_active:.2f} near_prior={g.frac_near_prior:.2f} "
            f"dev={g.mean_deviation:.2f}"
            for name, g in sorted(last.by_group.items())
        )
    )
