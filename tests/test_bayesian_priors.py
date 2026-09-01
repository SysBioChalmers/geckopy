"""Tests for kcat_sensitivity_analysis.bayesian.priors."""
import numpy as np
import pytest

from geckopy.adapter.params import BayesianParams, SourceGroupRule
from geckopy.kcat_sensitivity_analysis.bayesian.priors import (
    UNLABELLED_GROUP,
    SpikeSlabRV,
    build_kcat_prior,
    build_kcat_sparsity_prior,
    build_sigma0_log,
    classify_kcat_source,
    classify_kcat_sources,
)


def _params() -> BayesianParams:
    return BayesianParams(
        sigma0_log_default=0.5,
        source_groups={
            "dlkcat": SourceGroupRule(sources=["dlkcat"]),
            "brenda": SourceGroupRule(sources=["brenda"]),
            "okp": SourceGroupRule(match_okp=True),
        },
        sigma0_log_source={"dlkcat": 0.4, "brenda": 0.2, "okp": 0.3},
        shrink_thr_source={"dlkcat": 1.5, "brenda": 3.5, "okp": 2.0},
        force_prior_thr_source={"dlkcat": -1.0, "brenda": 4.0, "okp": 1.0},
    )


# --------------------------------------------------------------------------- #
# classify_kcat_source(s)
# --------------------------------------------------------------------------- #

def test_classify_matches_explicit_source_case_insensitively():
    params = _params()
    assert classify_kcat_source("DLKcat", params) == "dlkcat"
    assert classify_kcat_source("brenda", params) == "brenda"
    assert classify_kcat_source("BRENDA", params) == "brenda"


def test_classify_unmatched_source_is_unlabelled():
    params = _params()
    assert classify_kcat_source("some_other_source", params) == UNLABELLED_GROUP


def test_classify_matches_okp_method():
    params = _params()
    assert classify_kcat_source("CataPro", params, okp_method="CataPro") == "okp"
    # Case-insensitive here too.
    assert classify_kcat_source("catapro", params, okp_method="CataPro") == "okp"
    # No okp_method configured -> falls through to unlabelled.
    assert classify_kcat_source("CataPro", params, okp_method=None) == UNLABELLED_GROUP


def test_classify_kcat_sources_vectorised():
    params = _params()
    groups = classify_kcat_sources(
        ["dlkcat", "brenda", "CataPro", "unknown"], params, okp_method="CataPro",
    )
    assert list(groups) == ["dlkcat", "brenda", "okp", UNLABELLED_GROUP]


# --------------------------------------------------------------------------- #
# build_sigma0_log
# --------------------------------------------------------------------------- #

def test_build_sigma0_log_uses_group_values_and_default_fallback():
    params = _params()
    groups = np.array(["dlkcat", "brenda", UNLABELLED_GROUP, "okp"], dtype=object)

    sigma0_log = build_sigma0_log(groups, params)

    assert sigma0_log.tolist() == pytest.approx([0.4, 0.2, 0.5, 0.3])


# --------------------------------------------------------------------------- #
# build_kcat_prior
# --------------------------------------------------------------------------- #

def test_build_kcat_prior_mean_matches_kcat0():
    kcat0 = np.array([1.0, 10.0])
    sigma0_log = np.array([0.4, 0.2])

    prior = build_kcat_prior(kcat0, sigma0_log)

    # The prior's mean (not median) should sit at kcat0 -- this is
    # MATLAB's `muLog = log(kcats) - 0.5*sigma0log^2` bias correction,
    # verified via the underlying scipy distribution's .mean().
    for i, k0 in enumerate(kcat0):
        rv = prior[f"k{i}"]
        assert rv.distribution.mean() == pytest.approx(k0, rel=1e-6)
        # And the median sits strictly below the mean (a lognormal is
        # right-skewed), confirming the correction isn't a no-op.
        assert rv.distribution.median() < k0


def test_build_kcat_prior_sampling_is_lognormal_and_positive():
    rng_state = np.random.get_state()
    try:
        np.random.seed(0)
        kcat0 = np.array([5.0])
        sigma0_log = np.array([0.3])
        prior = build_kcat_prior(kcat0, sigma0_log)
        samples = np.array([prior.rvs()["k0"] for _ in range(2000)])
    finally:
        np.random.set_state(rng_state)

    assert np.all(samples > 0)
    # log-samples should be approximately normal around log(kcat0) - 0.5*sigma^2.
    expected_mu = np.log(5.0) - 0.5 * 0.3**2
    assert np.mean(np.log(samples)) == pytest.approx(expected_mu, abs=0.05)
    assert np.std(np.log(samples)) == pytest.approx(0.3, abs=0.05)


def test_build_kcat_prior_rejects_nonpositive_kcat0():
    with pytest.raises(ValueError, match="strictly positive"):
        build_kcat_prior(np.array([1.0, 0.0]), np.array([0.4, 0.4]))


def test_build_kcat_prior_rejects_shape_mismatch():
    with pytest.raises(ValueError, match="same shape"):
        build_kcat_prior(np.array([1.0, 2.0]), np.array([0.4]))


# --------------------------------------------------------------------------- #
# SpikeSlabRV / build_kcat_sparsity_prior
# --------------------------------------------------------------------------- #

def test_spike_slab_rv_pdf_is_weighted_mixture():
    rv = SpikeSlabRV(2.0, sigma_slab_log=0.4, spike_weight=0.7, spike_sigma_log=0.02)

    x = 2.0
    expected = 0.7 * rv._spike.pdf(x) + 0.3 * rv._slab.pdf(x)
    assert rv.pdf(x) == pytest.approx(expected)


def test_spike_slab_rv_concentrates_near_kcat0_when_spike_dominant():
    rng_state = np.random.get_state()
    try:
        np.random.seed(0)
        rv = SpikeSlabRV(3.0, sigma_slab_log=0.5, spike_weight=0.95, spike_sigma_log=0.01)
        samples = np.array([rv.rvs() for _ in range(1000)])
    finally:
        np.random.set_state(rng_state)

    # With 95% spike weight and a tight spike, most samples should sit
    # very close to kcat0=3.0.
    frac_near = np.mean(np.abs(samples - 3.0) < 0.2)
    assert frac_near > 0.85


def test_spike_slab_rv_rejects_invalid_weight():
    with pytest.raises(ValueError, match="spike_weight"):
        SpikeSlabRV(1.0, sigma_slab_log=0.4, spike_weight=1.5)


def test_build_kcat_sparsity_prior_shape_and_positivity():
    kcat0 = np.array([1.0, 2.0, 3.0])
    sigma0_log = np.array([0.4, 0.3, 0.2])

    prior = build_kcat_sparsity_prior(kcat0, sigma0_log, spike_weight=0.5)

    assert set(prior.keys()) == {"k0", "k1", "k2"}
    for key in prior.keys():
        assert isinstance(prior[key], SpikeSlabRV)
    sample = prior.rvs()
    assert all(sample[k] > 0 for k in ("k0", "k1", "k2"))
