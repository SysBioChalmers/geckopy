"""Tests for kcat_sensitivity_analysis.bayesian.priors."""
import numpy as np
import pytest

from geckopy.adapter.params import BayesianParams, SourceGroupRule
from geckopy.kcat_sensitivity_analysis.bayesian.priors import (
    UNLABELLED_GROUP,
    build_kcat_prior,
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


def test_proposal_sigma_log_defaults_to_the_prior_width():
    from geckopy.adapter.params import BayesianParams, SourceGroupRule
    from geckopy.kcat_sensitivity_analysis.bayesian.priors import (
        build_proposal_sigma_log, build_sigma0_log)

    params = BayesianParams(
        sigma0_log_default=0.3,
        source_groups={"brenda": SourceGroupRule(sources=["brenda"]),
                       "custom": SourceGroupRule(sources=["custom"])},
        sigma0_log_source={"brenda": 0.2, "custom": 0.1},
        shrink_thr_source={"brenda": 3.5, "custom": 5.5},
        force_prior_thr_source={"brenda": 4.0, "custom": 8.0},
    )
    groups = np.array(["brenda", "custom", "unlabelled"])
    assert np.array_equal(build_proposal_sigma_log(groups, params),
                          build_sigma0_log(groups, params))


def test_proposal_sigma_log_overrides_leave_the_prior_alone():
    """Widening the prior widens every proposal unless the two are split."""
    from geckopy.adapter.params import BayesianParams, SourceGroupRule
    from geckopy.kcat_sensitivity_analysis.bayesian.priors import (
        build_proposal_sigma_log, build_sigma0_log)

    params = BayesianParams(
        sigma0_log_default=2.5,
        source_groups={"brenda": SourceGroupRule(sources=["brenda"]),
                       "custom": SourceGroupRule(sources=["custom"])},
        sigma0_log_source={"brenda": 1.5, "custom": 0.4},
        shrink_thr_source={"brenda": 3.5, "custom": 5.5},
        force_prior_thr_source={"brenda": 4.0, "custom": 8.0},
        proposal_sigma_log_default=0.3,
        proposal_sigma_log_source={"brenda": 0.2},
    )
    groups = np.array(["brenda", "custom", "unlabelled"])

    # The prior keeps the measured widths.
    assert list(build_sigma0_log(groups, params)) == [1.5, 0.4, 2.5]
    # Proposals use the override, and custom -- absent from it -- keeps
    # its sigma0_log_source value rather than the proposal default.
    assert list(build_proposal_sigma_log(groups, params)) == [0.2, 0.4, 0.3]


def test_proposal_sigma_log_source_rejects_unknown_groups():
    import pytest
    from geckopy.adapter.params import BayesianParams, SourceGroupRule

    with pytest.raises(ValueError, match="not source_groups keys"):
        BayesianParams(
            source_groups={"brenda": SourceGroupRule(sources=["brenda"])},
            sigma0_log_source={"brenda": 0.2},
            shrink_thr_source={"brenda": 3.5},
            force_prior_thr_source={"brenda": 4.0},
            proposal_sigma_log_source={"typo": 0.2},
        )
