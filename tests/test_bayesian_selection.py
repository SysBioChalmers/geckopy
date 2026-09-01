"""Tests for kcat_sensitivity_analysis.bayesian.selection."""
import numpy as np
import pytest

from geckopy.kcat_sensitivity_analysis.bayesian.selection import (
    SelectionResult,
    next_quantile_epsilon,
    quantile_epsilon_select,
    truncation_select,
)


# --------------------------------------------------------------------------- #
# truncation_select
# --------------------------------------------------------------------------- #

def test_truncation_select_keeps_exact_top_fraction():
    distances = np.array([5.0, 1.0, 4.0, 2.0, 3.0, 9.0, 8.0, 7.0, 6.0, 0.5])  # n=10

    result = truncation_select(distances, min_keep=0.3)

    assert isinstance(result, SelectionResult)
    assert len(result.accepted_idx) == 3  # floor(0.3 * 10)
    accepted_values = sorted(distances[result.accepted_idx])
    assert accepted_values == [0.5, 1.0, 2.0]
    assert result.epsilon == pytest.approx(2.0)


def test_truncation_select_always_keeps_at_least_one():
    distances = np.array([3.0, 1.0, 2.0])

    result = truncation_select(distances, min_keep=0.01)

    assert len(result.accepted_idx) == 1
    assert distances[result.accepted_idx[0]] == pytest.approx(1.0)


def test_truncation_select_min_keep_one_keeps_everything():
    distances = np.array([3.0, 1.0, 2.0])

    result = truncation_select(distances, min_keep=1.0)

    assert sorted(result.accepted_idx.tolist()) == [0, 1, 2]
    assert result.epsilon == pytest.approx(3.0)


def test_truncation_select_rejects_invalid_min_keep():
    with pytest.raises(ValueError, match="min_keep"):
        truncation_select(np.array([1.0, 2.0]), min_keep=0.0)
    with pytest.raises(ValueError, match="min_keep"):
        truncation_select(np.array([1.0, 2.0]), min_keep=1.5)


def test_truncation_select_rejects_empty_distances():
    with pytest.raises(ValueError, match="empty"):
        truncation_select(np.array([]), min_keep=0.3)


def test_truncation_select_matches_matlab_worked_example():
    """The MATLAB acceptance step, with targetAccept's dead percentile
    cut removed: rmse = [newRmse; rmseTop], keep the
    floor(minKeep*numel(rmse)) smallest, epsilon = the worst kept."""
    new_rmse = np.array([0.9, 0.5, 1.2, 0.3, 0.8])
    rmse_top = np.array([0.4, 0.6])
    combined = np.concatenate([new_rmse, rmse_top])  # n=7
    min_keep = 0.3  # floor(0.3*7) = 2

    result = truncation_select(combined, min_keep=min_keep)

    assert len(result.accepted_idx) == 2
    assert sorted(combined[result.accepted_idx]) == [0.3, 0.4]
    assert result.epsilon == pytest.approx(0.4)


# --------------------------------------------------------------------------- #
# quantile_epsilon_select
# --------------------------------------------------------------------------- #

def test_quantile_epsilon_select_accepts_below_fixed_threshold():
    distances = np.array([0.1, 0.5, 0.9, 1.5, 2.0])

    result = quantile_epsilon_select(distances, epsilon=1.0)

    assert sorted(result.accepted_idx.tolist()) == [0, 1, 2]
    assert result.epsilon == pytest.approx(1.0)


def test_quantile_epsilon_select_accepted_count_is_not_fixed():
    """Unlike truncation_select, the accepted count here is whatever
    naturally passes -- it can be empty or everything, since epsilon
    was fixed *before* seeing this batch (the "prospective" property)."""
    distances = np.array([5.0, 6.0, 7.0])
    assert len(quantile_epsilon_select(distances, epsilon=0.1).accepted_idx) == 0
    assert len(quantile_epsilon_select(distances, epsilon=100.0).accepted_idx) == 3


# --------------------------------------------------------------------------- #
# next_quantile_epsilon
# --------------------------------------------------------------------------- #

def test_next_quantile_epsilon_matches_unweighted_median():
    accepted = np.array([1.0, 2.0, 3.0, 4.0, 5.0])

    eps = next_quantile_epsilon(accepted, alpha=0.5)

    assert eps == pytest.approx(np.median(accepted))


def test_next_quantile_epsilon_weighted_shifts_toward_heavier_weight():
    accepted = np.array([1.0, 10.0])
    unweighted = next_quantile_epsilon(accepted, alpha=0.5)
    weighted = next_quantile_epsilon(
        accepted, alpha=0.5, weights=np.array([0.9, 0.1]),
    )

    # Heavier weight on the smaller value should pull the weighted
    # median down relative to the unweighted one.
    assert weighted < unweighted


def test_next_quantile_epsilon_rejects_empty():
    with pytest.raises(ValueError, match="empty"):
        next_quantile_epsilon(np.array([]))


# --------------------------------------------------------------------------- #
# Prospective vs retrospective: the essential Axis-1 semantic difference
# --------------------------------------------------------------------------- #

def test_prospective_epsilon_is_independent_of_current_batch():
    """quantile_epsilon_select's threshold comes from the *previous*
    generation and does not adapt to whatever this batch looks like --
    the defining difference from truncation_select's retrospective cut,
    which always looks at the current batch's own distribution."""
    prev_gen_accepted = np.array([0.2, 0.3, 0.4])
    epsilon = next_quantile_epsilon(prev_gen_accepted, alpha=0.5)  # 0.3

    # A batch that's entirely worse than the fixed epsilon accepts nothing,
    # even though truncation_select would still keep its best fraction.
    bad_batch = np.array([10.0, 20.0, 30.0])
    assert len(quantile_epsilon_select(bad_batch, epsilon=epsilon).accepted_idx) == 0
    assert len(truncation_select(bad_batch, min_keep=0.3).accepted_idx) == 1
