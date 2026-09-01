"""Real-scale smoke test for Bayesian kcat tuning.

Runs both selection variants against the real
``ecYeastGEM`` tutorial model and the real GECKO
``bayesianFluxData.tsv``/``bayesianMaxGrowth.tsv``/``bayesianZeroExch.tsv``
tutorial data (ported verbatim from the GECKO MATLAB repo's
``tutorials/full_ecModel/data/``), per the plan's Test Strategy #4.

The model must be a *pre-tuning* ecModel -- the one MATLAB feeds to
``bayesianSensitivityTuning``. Note that this is NOT what
``tutorials/full_ecModel/protocol.py`` writes as ``ecYeastGEM.yml``:
that file is a Stage 3 output, i.e. already past ``sensitivity_tuning``
(which bumps the most limiting kcats, lifting max growth from
0.0013 to 0.43 /h) plus the manual ``r_0079`` curation. Bayesian
tuning is meant to *replace* that step, so starting from it would
bake in the very hand-fix under test and confound the
variant-to-variant RMSE comparison.

Default input is therefore ``ecYeastGEM_preTune_GECKOderived.yml``:
the GECKO MATLAB tutorial's own pre-Bayesian ecModel
(``GECKO/tutorials/full_ecModel/models/ecYeastGEM.yml``, ``develop4``)
copied into ``tutorials/full_ecModel/models/``. It is gitignored like
every other ``ecYeastGEM*.yml`` build output. Mirrors
``tests/test_light_humangem_smoke.py``'s style: skipped gracefully
when it isn't present -- copy it from a GECKO checkout, or point
``ECYEASTGEM_YML`` at another pre-tuning ecModel.

``_SMOKE_PARAMS`` mirrors that model's own MATLAB adapter
(``GECKO/tutorials/full_ecModel/YeastGEMAdapter.m``) rather than
``BayesianParams``' generic defaults, since the two disagree on
exactly the load-bearing settings: MATLAB configures
``kcatSources = {'OpenKineticsPredictor','brenda','custom'}`` for this
model, whereas the port's default groups are ``dlkcat``/``brenda``/
``custom`` -- which would leave this model's 1175 OKP kcats (plus
``standard``/``isozymes``) unclassified and silently on the
``*_default`` tier instead of their own sticky-prior tier. Only the
schedule is shrunk (see below).

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
from geckopy.adapter.params import BayesianParams, SourceGroupRule
from geckopy.kcat_sensitivity_analysis.bayesian.data import load_bayesian_data
from geckopy.kcat_sensitivity_analysis.bayesian.tuning import bayesian_kcat_tuning

_TUTORIAL_DIR = Path(__file__).resolve().parent.parent / "tutorials" / "full_ecModel"
_DEFAULT_MODEL = "ecYeastGEM_preTune_GECKOderived.yml"
_ECYEASTGEM_YML = Path(
    os.environ.get("ECYEASTGEM_YML", _TUTORIAL_DIR / "models" / _DEFAULT_MODEL)
)

pytestmark = [
    pytest.mark.smoke,
    pytest.mark.skipif(
        not _ECYEASTGEM_YML.is_file(),
        reason=(
            f"pre-tuning ecModel not found at {_ECYEASTGEM_YML} (a "
            "gitignored build output). Copy GECKO's own "
            "tutorials/full_ecModel/models/ecYeastGEM.yml there under that "
            "name, or set ECYEASTGEM_YML to another pre-tuning ecModel. "
            "Note protocol.py's own ecYeastGEM.yml is post-sensitivity-"
            "tuning and is NOT a valid input here -- see this module's "
            "docstring."
        ),
    ),
]

# Hyperparameters from GECKO/tutorials/full_ecModel/YeastGEMAdapter.m
# (the MATLAB adapter this very model was built with), except for the
# schedule, which is deliberately tiny -- see module docstring. MATLAB's
# dead/cosmetic knobs (targetAccept, varianceCap*, *Tight) have no
# counterpart here by design, per the plan's Q2.
_SMOKE_PARAMS = BayesianParams(
    sigma0_log_default=0.3,
    source_groups={
        # MATLAB: kcatSources = {'OpenKineticsPredictor','brenda','custom'}.
        # match_okp also catches a geckopy-built model, where OKP kcats
        # carry the raw predictor method name (OkpParams.method) instead.
        "okp": SourceGroupRule(sources=["OpenKineticsPredictor"], match_okp=True),
        "brenda": SourceGroupRule(sources=["brenda"]),
        "custom": SourceGroupRule(sources=["custom"]),
    },
    sigma0_log_source={"okp": 0.25, "brenda": 0.2, "custom": 0.1},
    shrink_thr_default=5.0,
    shrink_thr_source={"okp": 3.0, "brenda": 10.0, "custom": 12.0},
    force_prior_thr_default=3.5,
    force_prior_thr_source={"okp": 2.5, "brenda": 11.0, "custom": 13.0},
    sparsity_threshold=0.5,
    schedule_generations=[1],
    schedule_samples=[8],
    min_keep=0.3,
    max_keep=0.6,
    rmse_threshold=-1.0,  # unreachable -> always runs exactly max_generations
    max_generations=2,
)

# Scoring is parallelised across particles (one persistent EcModel per
# worker), so more workers than particles buys nothing -- cap at the
# per-generation sample count. Override with BAYESIAN_SMOKE_NPROC.
_N_PROC = int(
    os.environ.get(
        "BAYESIAN_SMOKE_NPROC",
        min(_SMOKE_PARAMS.schedule_samples[0], os.cpu_count() or 1),
    )
)

_SELECTIONS = ["truncation", "quantile_epsilon"]


def _load_model_and_data():
    adapter = ModelAdapter.from_folder(_TUTORIAL_DIR)
    model = load_ec_model(_ECYEASTGEM_YML, adapter=adapter)
    bay_data = load_bayesian_data(adapter)
    assert bay_data.flux_data is not None and bay_data.max_grate is not None, (
        "Expected both bayesianFluxData.tsv and bayesianMaxGrowth.tsv to be "
        "present for this smoke test."
    )
    return model, adapter, bay_data


@pytest.mark.parametrize("selection", _SELECTIONS)
def test_bayesian_tuning_runs_at_real_scale(selection):
    model, adapter, bay_data = _load_model_and_data()
    original_kcat = model.ec.kcat.copy()
    n_tunable = int((original_kcat > 0).sum())
    assert n_tunable > 1000, (
        f"Expected thousands of tunable kcats on the real model; got {n_tunable}."
    )

    t0 = time.perf_counter()
    result = bayesian_kcat_tuning(
        model, adapter=adapter, params=_SMOKE_PARAMS, bay_data=bay_data,
        selection=selection,
        n_proc=_N_PROC, seed=0, verbose=False,
    )
    wall_time = time.perf_counter() - t0

    assert result.n_generations == 2
    assert len(result.rmse_trace) == 2
    assert all(np.isfinite(r) for r in result.rmse_trace)
    assert len(result.diagnostics_trace) == 2

    last = result.diagnostics_trace[-1]
    # Every configured group is represented on this model (the
    # unlabelled tier -- `standard`/`isozymes` kcats -- may also appear).
    assert set(last.by_group) >= set(_SMOKE_PARAMS.source_groups)
    for group_diag in last.by_group.values():
        assert 0.0 <= group_diag.frac_active <= 1.0
        assert 0.0 <= group_diag.frac_near_prior <= 1.0

    print(
        f"\n[smoke:{selection}] "
        f"wall_time={wall_time:.1f}s n_proc={_N_PROC} "
        f"solver={model.solver.interface.__name__.rsplit('.', 1)[-1]} "
        f"rmse_trace={[f'{r:.3g}' for r in result.rmse_trace]} "
        f"n_tunable={n_tunable} "
        + " | ".join(
            f"{name}: active={g.frac_active:.2f} near_prior={g.frac_near_prior:.2f} "
            f"dev={g.mean_deviation:.2f}"
            for name, g in sorted(last.by_group.items())
        )
    )

