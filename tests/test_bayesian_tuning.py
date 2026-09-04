"""Tests for kcat_sensitivity_analysis.bayesian.tuning.

A tiny two-enzyme toy EcModel with one "brenda" (trusted, tight prior)
and one "dlkcat" (untrusted, loose prior) kcat, each gating its own
independent branch so each condition's growth depends on exactly one
of them. Both start equally "wrong" (half the true value) so any
difference in how far each one moves during tuning is attributable to
the trust-tier prior/regularization machinery, not to differing
amounts of initial error.

Run once per combination of the two axes (4 runs total), per the
plan's Test Strategy #3.
"""
from pathlib import Path

import cobra
import numpy as np
import pytest

from geckopy import EcModel, ModelAdapter
from geckopy.adapter.params import BayesianParams
from geckopy.databases.flux_data import FluxData
from geckopy.ec_model.ec_data import EcData
from geckopy.kcat_sensitivity_analysis.bayesian.data import BayesianData
from geckopy.kcat_sensitivity_analysis.bayesian.tuning import (
    BayesianTuningResult,
    bayesian_kcat_tuning,
)
from scipy import sparse

_TRUE_KCAT = 2.0
_START_KCAT = 1.0  # both branches start at half the true value
_MW = 100.0
# growth = kcat * 3600 / mw (pool size 1 mg/gDW); at the true kcat,
# growth = 72, matching both the measured growth rate and (1:1
# stoichiometry) the measured carbon uptake below.
_TRUE_GROWTH = _TRUE_KCAT * 3600.0 / _MW


def _adapter(tmp_path: Path) -> ModelAdapter:
    (tmp_path / "model_adapter.toml").write_text(
        'conv_gem = "dummy.xml"\n'
        'org_name = "test"\n'
        'bio_rxn = "biomass"\n'
    )
    return ModelAdapter.from_folder(tmp_path)


def _build_toy(adapter: ModelAdapter) -> EcModel:
    model = EcModel("toy", adapter=adapter)

    glc_e = cobra.Metabolite("glc_e", compartment="e")
    eth_e = cobra.Metabolite("eth_e", compartment="e")
    prot_pool = cobra.Metabolite("prot_pool", compartment="c")
    prot_glc = cobra.Metabolite("prot_Eglc", compartment="c")
    prot_eth = cobra.Metabolite("prot_Eeth", compartment="c")
    bio_met = cobra.Metabolite("bio_met", compartment="c")
    model.add_metabolites([glc_e, eth_e, prot_pool, prot_glc, prot_eth, bio_met])

    EX_glc = cobra.Reaction("EX_glc")
    EX_glc.add_metabolites({glc_e: -1.0})
    EX_glc.bounds = (-1000.0, 0.0)

    EX_eth = cobra.Reaction("EX_eth")
    EX_eth.add_metabolites({eth_e: -1.0})
    EX_eth.bounds = (-1000.0, 0.0)

    coeff = _MW / (_START_KCAT * 3600.0)
    R_glc = cobra.Reaction("R_glc")
    R_glc.add_metabolites({glc_e: -1.0, prot_glc: -coeff, bio_met: 1.0})
    R_glc.bounds = (0.0, 1000.0)

    R_eth = cobra.Reaction("R_eth")
    R_eth.add_metabolites({eth_e: -1.0, prot_eth: -coeff, bio_met: 1.0})
    R_eth.bounds = (0.0, 1000.0)

    BIO = cobra.Reaction("biomass")
    BIO.add_metabolites({bio_met: -1.0})
    BIO.bounds = (0.0, 1000.0)

    pool_ex = cobra.Reaction("prot_pool_exchange")
    pool_ex.add_metabolites({prot_pool: 1.0})
    pool_ex.bounds = (0.0, 1.0)

    usage_glc = cobra.Reaction("usage_prot_Eglc")
    usage_glc.add_metabolites({prot_pool: -1.0, prot_glc: 1.0})
    usage_glc.bounds = (0.0, 1000.0)

    usage_eth = cobra.Reaction("usage_prot_Eeth")
    usage_eth.add_metabolites({prot_pool: -1.0, prot_eth: 1.0})
    usage_eth.bounds = (0.0, 1000.0)

    model.add_reactions(
        [EX_glc, EX_eth, R_glc, R_eth, BIO, pool_ex, usage_glc, usage_eth]
    )
    model.objective = "biomass"

    model.ec = EcData(
        rxns=["R_glc", "R_eth"],
        kcat=np.array([_START_KCAT, _START_KCAT]),
        source=["brenda", "dlkcat"],
        notes=["", ""],
        eccodes=["", ""],
        genes=["g_glc", "g_eth"],
        enzymes=["Eglc", "Eeth"],
        mw=np.array([_MW, _MW]),
        sequence=["", ""],
        concs=np.array([np.nan, np.nan]),
        rxn_enz_mat=sparse.csr_matrix(np.eye(2)),
    )
    return model


def _bay_data() -> BayesianData:
    # max_grate (constrain=False): uptake is always opened fully
    # (-1000, ignoring the exch_fluxes values below -- they only
    # decide *which* exchange is unblocked per condition), so growth
    # is purely enzyme-limited: 36*kcat, strictly increasing with no
    # saturation. This keeps the objective unimodal with a single
    # minimum at kcat=2.0 in both directions -- unlike flux_data
    # (constrain=True), which would also compare simulated vs.
    # measured *exchange* flux for the active condition's own column;
    # since that column doubles as the uptake bound, any measured
    # value tight enough to be realistic would cap achievable growth
    # for kcat > true_kcat too, creating a flat "any kcat >= true_kcat
    # scores 0" plateau with no restoring force -- not what this test
    # wants to isolate.
    max_grate = FluxData(
        conds=["glucose", "ethanol"],
        p_tot=np.array([np.nan, np.nan]),
        gr_rate=np.array([_TRUE_GROWTH, _TRUE_GROWTH]),
        exch_fluxes=np.array(
            [
                [-1000.0, np.nan],
                [np.nan, -1000.0],
            ]
        ),
        exch_mets=["glucose", "ethanol"],
        exch_rxn_ids=["EX_glc", "EX_eth"],
    )
    return BayesianData(flux_data=None, max_grate=max_grate, zero_flux=[])


_SELECTIONS = ["truncation", "quantile_epsilon"]


_SEEDS = [0, 1, 2]


@pytest.mark.parametrize("selection", _SELECTIONS)
def test_trusted_source_moves_less_than_untrusted_for_every_combination(
    tmp_path, selection,
):
    params = BayesianParams(
        # Exaggerated trust contrast: what is under test is that
        # per-source sigma0_log transmits to how far a kcat moves, not
        # that the shipped defaults are well chosen. The shipped 0.2/0.4
        # contrast is within sampling noise on a 2-parameter toy: the
        # shrink/force-prior blend never feeds back into sampling (see
        # tuning.py's module docstring), so the prior draw and the
        # transition bandwidth are the only channel. At 0.05/1.0 both
        # selection variants separate well clear of that noise.
        sigma0_log_source={"brenda": 0.05, "dlkcat": 1.0, "custom": 0.1},
        schedule_generations=[1],
        schedule_samples=[40],
        min_keep=0.3,
        max_keep=0.6,
        rmse_threshold=-1.0,  # unreachable -> always runs exactly max_generations
        max_generations=4,
    )

    brenda_moves = []
    dlkcat_moves = []
    for seed in _SEEDS:
        adapter = _adapter(tmp_path)
        model = _build_toy(adapter)
        bay_data = _bay_data()
        apply_kcat_constraints_before = model.ec.kcat.copy()

        result = bayesian_kcat_tuning(
            model, adapter=adapter, params=params, bay_data=bay_data,
            selection=selection,
            seed=seed, verbose=False,
        )

        assert isinstance(result, BayesianTuningResult)
        assert result.rxns == ["R_glc", "R_eth"]
        assert result.groups == ["brenda", "dlkcat"]
        assert np.array_equal(result.old_kcat, apply_kcat_constraints_before)
        assert result.n_generations == 4

        # model.ec.kcat was actually mutated in place to match the result.
        assert np.array_equal(model.ec.kcat, result.new_kcat)

        # Diagnostics/rmse trace present for every generation.
        assert len(result.rmse_trace) == 4
        assert len(result.diagnostics_trace) == 4
        assert len(result.posterior_trace) == 4

        brenda_moves.append(abs(np.log(result.new_kcat[0]) - np.log(_START_KCAT)))
        dlkcat_moves.append(abs(np.log(result.new_kcat[1]) - np.log(_START_KCAT)))

    # Primary comparison, averaged over a few seeds to smooth out
    # single-run sampling noise (especially for selection=
    # "quantile_epsilon", which redraws fresh each generation with no
    # elitist carry-over): does the trusted (tight-prior) kcat move
    # less, in absolute log-space terms, than the untrusted one, given
    # an identical starting error for both?
    mean_brenda_move = float(np.mean(brenda_moves))
    mean_dlkcat_move = float(np.mean(dlkcat_moves))
    assert mean_brenda_move < mean_dlkcat_move, (
        f"[{selection}] expected brenda (trusted) to move "
        f"less than dlkcat (untrusted) in log-space, averaged over seeds "
        f"{_SEEDS}; got mean_brenda_move={mean_brenda_move:.4f}, "
        f"mean_dlkcat_move={mean_dlkcat_move:.4f} "
        f"(per-seed: brenda={brenda_moves}, dlkcat={dlkcat_moves})"
    )


@pytest.mark.parametrize("selection", _SELECTIONS)
def test_parallel_scoring_matches_serial(tmp_path, selection):
    """n_proc=2 must reproduce n_proc=1 bit-for-bit for the same seed.

    Sampling stays single-threaded in the main process (see tuning.py's
    module docstring): only the already-fixed batch of particles each
    generation gets scored in parallel, and scoring is a pure function
    of the kcat vector -- so the *set* of particles proposed, and every
    particle's RMSE, must be identical regardless of how many workers
    scored them or in what order.
    """
    params = BayesianParams(
        schedule_generations=[1],
        schedule_samples=[20],
        min_keep=0.3,
        max_keep=0.6,
        rmse_threshold=-1.0,
        max_generations=3,
    )

    def _run(n_proc):
        adapter = _adapter(tmp_path)
        model = _build_toy(adapter)
        bay_data = _bay_data()
        return bayesian_kcat_tuning(
            model, adapter=adapter, params=params, bay_data=bay_data,
            selection=selection,
            n_proc=n_proc, seed=0, verbose=False,
        )

    serial = _run(1)
    parallel = _run(2)

    assert np.array_equal(serial.new_kcat, parallel.new_kcat)
    assert serial.rmse_trace == parallel.rmse_trace
    assert serial.n_generations == parallel.n_generations
    assert serial.converged == parallel.converged


def test_no_tunable_kcats_raises(tmp_path):
    adapter = _adapter(tmp_path)
    model = _build_toy(adapter)
    model.ec.kcat[:] = 0.0
    bay_data = _bay_data()

    with pytest.raises(ValueError, match="No tunable kcats"):
        bayesian_kcat_tuning(model, adapter=adapter, bay_data=bay_data)


def test_missing_bay_data_raises(tmp_path):
    adapter = _adapter(tmp_path)
    model = _build_toy(adapter)
    empty = BayesianData(flux_data=None, max_grate=None, zero_flux=[])

    with pytest.raises(ValueError, match="nothing to tune"):
        bayesian_kcat_tuning(model, adapter=adapter, bay_data=empty)


def test_kcat_bounds_keep_the_ordinary_window():
    """An ordinary kcat proposes within 1e-2..1e4 1/s."""
    from geckopy.kcat_sensitivity_analysis.bayesian.tuning import kcat_bounds

    lo, hi = kcat_bounds(np.array([1.0, 50.0]))

    np.testing.assert_allclose(lo, [1e-2, 1e-2])
    np.testing.assert_allclose(hi, [1e4, 1e4])


def test_kcat_bounds_always_contain_the_prior():
    """A prior outside its own bounds cannot be proposed, so the search
    moves it before considering any evidence."""
    from geckopy.kcat_sensitivity_analysis.bayesian.tuning import kcat_bounds

    kcat0 = np.array([1.6e-4, 1.9e-3, 1.0, 5e4])
    lo, hi = kcat_bounds(kcat0)

    assert np.all(lo <= kcat0) and np.all(kcat0 <= hi)
    # And with room to move a hundred-fold either way.
    np.testing.assert_allclose(lo, [1.6e-6, 1.9e-5, 1e-2, 1e-2])
    np.testing.assert_allclose(hi, [1e4, 1e4, 1e4, 5e6])


def test_adapt_proposal_scale_follows_the_acceptance_rate():
    """Below target the proposal shrinks, above it grows, and the scale
    stays inside its bounds."""
    from geckopy.kcat_sensitivity_analysis.bayesian.tuning import adapt_proposal_scale

    params = BayesianParams(
        adapt_proposal_width=True, target_accept_rate=0.15,
        proposal_adaptation_rate=2.0, proposal_scale_bounds=(0.02, 2.0),
    )

    assert adapt_proposal_scale(1.0, 0.15, params) == pytest.approx(1.0)
    assert adapt_proposal_scale(1.0, 0.005, params) < 1.0
    assert adapt_proposal_scale(1.0, 0.50, params) > 1.0

    # MATLAB's own collapsing acceptance rate drives the scale down, and
    # the clamp holds it above the floor however long that runs.
    scale = 1.0
    for rate in (0.10, 0.09, 0.07, 0.05, 0.03, 0.02, 0.005, 0.005, 0.005):
        scale = adapt_proposal_scale(scale, rate, params)
    assert 0.02 <= scale < 1.0

    for _ in range(200):
        scale = adapt_proposal_scale(scale, 0.0, params)
    assert scale == pytest.approx(0.02)

    scale = 1.0
    for _ in range(200):
        scale = adapt_proposal_scale(scale, 1.0, params)
    assert scale == pytest.approx(2.0)


def test_adapt_proposal_width_is_off_by_default():
    assert BayesianParams().adapt_proposal_width is False


def test_proposal_acceptance_rate_divides_by_the_number_proposed():
    """The denominator is the proposals made, not the accepted set: the
    latter measures the accepted set's composition and sits near 1.0 in
    early generations however poor the proposals are."""
    from geckopy.kcat_sensitivity_analysis.bayesian.tuning import (
        proposal_acceptance_rate,
    )

    # 1000 proposals, 300 accepted, every one of them new.
    accepted = np.arange(300)
    assert proposal_acceptance_rate(accepted, 1000) == pytest.approx(0.30)

    # Half the accepted set carried over from previous generations.
    accepted = np.concatenate([np.arange(150), 1000 + np.arange(150)])
    assert proposal_acceptance_rate(accepted, 1000) == pytest.approx(0.15)

    # Nothing new survived.
    assert proposal_acceptance_rate(np.array([1000, 1001]), 1000) == 0.0
    assert proposal_acceptance_rate(np.array([]), 0) == 0.0


def test_prior_penalty_leaves_the_reported_rmse_plain(tmp_path):
    """Selection may run on a penalised objective, but rmse_trace must
    stay the plain RMSE so runs remain comparable across penalties and
    against MATLAB."""
    def _params(weight):
        return BayesianParams(
            schedule_generations=[1], schedule_samples=[20],
            min_keep=0.3, max_keep=0.6, rmse_threshold=-1.0,
            max_generations=3, prior_penalty_weight=weight,
        )

    def _run(weight):
        adapter = _adapter(tmp_path)
        model = _build_toy(adapter)
        return bayesian_kcat_tuning(
            model, adapter=adapter, params=_params(weight),
            bay_data=_bay_data(), selection="truncation",
            n_proc=1, seed=0, verbose=False,
        )

    unpenalised = _run(0.0)
    # With no penalty the objective *is* the RMSE.
    assert unpenalised.objective_trace == pytest.approx(unpenalised.rmse_trace)

    # The toy's RMSE is ~700, so the weight has to be on that scale to
    # bite at all; on the full ecModel (RMSE ~1) it would be ~1.
    penalised = _run(5.0e3)
    assert len(penalised.objective_trace) == len(penalised.rmse_trace)
    # The penalty is non-negative and vanishes only at the prior, so the
    # objective sits at or above the RMSE it is built from.
    assert all(
        o >= r - 1e-9
        for o, r in zip(penalised.objective_trace, penalised.rmse_trace)
    )
    # And it changed the search: a penalty that altered nothing would be
    # a silently dead knob.
    assert penalised.rmse_trace != unpenalised.rmse_trace


def test_penalised_search_moves_kcats_less(tmp_path):
    """The point of the penalty: fewer/smaller departures from prior."""
    def _run(weight):
        adapter = _adapter(tmp_path)
        model = _build_toy(adapter)
        return bayesian_kcat_tuning(
            model, adapter=adapter,
            params=BayesianParams(
                schedule_generations=[1], schedule_samples=[30],
                min_keep=0.3, max_keep=0.6, rmse_threshold=-1.0,
                max_generations=4, prior_penalty_weight=weight,
            ),
            bay_data=_bay_data(), selection="truncation",
            n_proc=1, seed=0, verbose=False,
        )

    free = _run(0.0)
    penalised = _run(5.0e4)
    free_move = np.abs(np.log(free.new_kcat / free.old_kcat)).mean()
    pen_move = np.abs(np.log(penalised.new_kcat / penalised.old_kcat)).mean()
    assert pen_move < free_move


def test_tunable_mask_holds_excluded_kcats_at_their_prior(tmp_path):
    """Parameters the data cannot speak to should not be edited at all."""
    adapter = _adapter(tmp_path)
    model = _build_toy(adapter)
    before = model.ec.kcat.copy()

    mask = np.ones(len(model.ec.rxns), dtype=bool)
    mask[0] = False                       # exclude the first tunable kcat

    result = bayesian_kcat_tuning(
        model, adapter=adapter,
        params=BayesianParams(
            schedule_generations=[1], schedule_samples=[20],
            min_keep=0.3, max_keep=0.6, rmse_threshold=-1.0, max_generations=2,
        ),
        bay_data=_bay_data(), selection="truncation",
        tunable_mask=mask, n_proc=1, seed=0, verbose=False,
    )

    # Excluded row untouched in the model, and absent from the result.
    assert model.ec.kcat[0] == pytest.approx(before[0])
    assert len(result.rxns) == int(mask.sum())
    assert model.ec.rxns[0] not in result.rxns
    assert len(result.new_kcat) == len(result.old_kcat) == int(mask.sum())


def test_tunable_mask_rejects_a_wrong_shape_and_an_empty_selection(tmp_path):
    adapter = _adapter(tmp_path)
    model = _build_toy(adapter)
    params = BayesianParams(
        schedule_generations=[1], schedule_samples=[5],
        min_keep=0.3, max_keep=0.6, rmse_threshold=-1.0, max_generations=1,
    )
    with pytest.raises(ValueError, match="tunable_mask has shape"):
        bayesian_kcat_tuning(
            model, adapter=adapter, params=params, bay_data=_bay_data(),
            tunable_mask=np.ones(3, dtype=bool), n_proc=1, seed=0, verbose=False,
        )
    with pytest.raises(ValueError, match="excludes every one"):
        bayesian_kcat_tuning(
            model, adapter=adapter, params=params, bay_data=_bay_data(),
            tunable_mask=np.zeros(len(model.ec.rxns), dtype=bool),
            n_proc=1, seed=0, verbose=False,
        )
