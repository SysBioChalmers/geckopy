"""Tests for kcat_sensitivity_analysis.bayesian.posterior (Axis 2, variant A).

Uses 3 synthetic parameters -- one per trust tier -- so this doubles
as the first real trusted-vs-untrusted comparison point (per the
plan's Test Strategy #1): does a trusted source's point estimate move
less than an untrusted source's, given an *identical* raw pull from
the accepted samples?
"""
import numpy as np
import pytest

from geckopy.kcat_sensitivity_analysis.bayesian.posterior import (
    PosteriorUpdate,
    update_posterior_shrinkage,
)

KCAT0 = np.array([1.0, 1.0, 1.0])
SIGMA0_LOG = np.array([0.1, 0.5, 0.3])  # trusted, untrusted, unlabelled
GROUPS = np.array(["trusted", "untrusted", "unlabelled"], dtype=object)


def _kcat_top_for_deviation(dev_in_sigmas: float, n_accepted: int = 4) -> np.ndarray:
    """Every accepted particle identical, at exactly `dev_in_sigmas`
    prior-sigmas above kcat0 in log-space -- gives a deterministic
    mean and zero std, so the expected output is hand-computable
    without needing to reproduce update_posterior_shrinkage's own math."""
    mu_log = np.log(KCAT0) + dev_in_sigmas * SIGMA0_LOG
    column = np.exp(mu_log)
    return np.tile(column.reshape(-1, 1), (1, n_accepted))


def test_trusted_source_moves_less_than_untrusted_for_equal_raw_pull():
    kcat_top = _kcat_top_for_deviation(dev_in_sigmas=2.0)

    result = update_posterior_shrinkage(
        kcat_top, KCAT0, SIGMA0_LOG, GROUPS,
        shrink_thr_default=1.0,
        shrink_thr_source={"trusted": 5.0, "untrusted": 0.5},
        force_prior_thr_default=-1.0,
        force_prior_thr_source={"trusted": -1.0, "untrusted": -1.0},
        sparsity_threshold=0.1,
    )

    assert isinstance(result, PosteriorUpdate)
    rel_change = np.abs(result.kcats - KCAT0) / KCAT0
    # trusted (sticky, high shrink threshold) moves noticeably less
    # than untrusted (loose, low shrink threshold) for the same 2-sigma
    # raw pull.
    assert rel_change[0] < rel_change[1]

    # Hand-computed expected values (see module docstring's worked
    # example): shrink_weight = min(dev/thr, 1).
    expected_shrink = np.array([min(2.0 / 5.0, 1.0), min(2.0 / 0.5, 1.0), min(2.0 / 1.0, 1.0)])
    assert result.shrink_weight == pytest.approx(expected_shrink)

    mu_log = np.log(KCAT0) + 2.0 * SIGMA0_LOG
    expected_kcats = np.exp(
        expected_shrink * mu_log + (1 - expected_shrink) * np.log(KCAT0)
    )
    assert result.kcats == pytest.approx(expected_kcats)
    assert not result.snapped_to_prior.any()
    assert not result.forced_to_prior.any()


def test_force_prior_threshold_overrides_shrink_weight_to_zero():
    # A small deviation (0.5 sigma) that's still below the untrusted
    # group's force-prior threshold (1.0) should force shrink_weight
    # to exactly 0, discarding whatever the raw shrink-weight formula
    # would have given.
    kcat_top = _kcat_top_for_deviation(dev_in_sigmas=0.5)

    result = update_posterior_shrinkage(
        kcat_top, KCAT0, SIGMA0_LOG, GROUPS,
        shrink_thr_default=1.0,
        shrink_thr_source={"trusted": 5.0, "untrusted": 0.5},
        force_prior_thr_default=-1.0,
        force_prior_thr_source={"trusted": -1.0, "untrusted": 1.0},
        sparsity_threshold=0.0,  # disabled, isolate the force-prior effect
    )

    assert result.forced_to_prior[1]  # untrusted row forced
    assert result.shrink_weight[1] == pytest.approx(0.0)
    assert result.kcats[1] == pytest.approx(KCAT0[1])
    assert result.sigma_log[1] == pytest.approx(SIGMA0_LOG[1])
    # trusted row's force-prior threshold is disabled (-1.0) -> not forced.
    assert not result.forced_to_prior[0]


def test_sparsity_threshold_snaps_small_changes_to_exact_prior():
    # A tiny deviation, well within sparsity_threshold * sigma0_log,
    # should snap the blended estimate back to *exactly* kcat0/sigma0_log
    # even though shrink_weight itself is nonzero.
    kcat_top = _kcat_top_for_deviation(dev_in_sigmas=0.01)

    result = update_posterior_shrinkage(
        kcat_top, KCAT0, SIGMA0_LOG, GROUPS,
        shrink_thr_default=1.0,
        shrink_thr_source={"trusted": 5.0, "untrusted": 0.5},
        force_prior_thr_default=-1.0,
        force_prior_thr_source={"trusted": -1.0, "untrusted": -1.0},
        sparsity_threshold=0.5,  # generous -> should catch this tiny change
    )

    assert result.snapped_to_prior.all()
    assert result.kcats == pytest.approx(KCAT0)
    assert result.sigma_log == pytest.approx(SIGMA0_LOG)


def test_unlabelled_group_uses_default_thresholds():
    kcat_top = _kcat_top_for_deviation(dev_in_sigmas=2.0)

    result = update_posterior_shrinkage(
        kcat_top, KCAT0, SIGMA0_LOG, GROUPS,
        shrink_thr_default=10.0,  # very sticky default
        shrink_thr_source={"trusted": 5.0, "untrusted": 0.5},
        force_prior_thr_default=-1.0,
        force_prior_thr_source={"trusted": -1.0, "untrusted": -1.0},
        sparsity_threshold=0.0,
    )

    # unlabelled (index 2) should use shrink_thr_default=10.0 -> small
    # shrink weight -> stays close to kcat0.
    assert result.shrink_weight[2] == pytest.approx(min(2.0 / 10.0, 1.0))


def test_rejects_shape_mismatch():
    kcat_top = _kcat_top_for_deviation(dev_in_sigmas=1.0)
    with pytest.raises(ValueError, match="kcat0"):
        update_posterior_shrinkage(
            kcat_top, np.array([1.0, 1.0]), SIGMA0_LOG, GROUPS,
            shrink_thr_default=1.0, shrink_thr_source={},
            force_prior_thr_default=-1.0, force_prior_thr_source={},
            sparsity_threshold=0.1,
        )


def test_rejects_no_accepted_particles():
    empty = np.empty((3, 0))
    with pytest.raises(ValueError, match="no accepted particles"):
        update_posterior_shrinkage(
            empty, KCAT0, SIGMA0_LOG, GROUPS,
            shrink_thr_default=1.0, shrink_thr_source={},
            force_prior_thr_default=-1.0, force_prior_thr_source={},
            sparsity_threshold=0.1,
        )
