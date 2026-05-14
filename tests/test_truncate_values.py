"""Tests for truncate_values."""
import numpy as np
import pytest

from geckopy.kcat_sensitivity_analysis import truncate_values
from geckopy.kcat_sensitivity_analysis.truncate_values import _truncate_scalar


# --------------------------------------------------------------------------- #
# Scalar helper
# --------------------------------------------------------------------------- #

def test_scalar_zero_passes_through():
    assert _truncate_scalar(0.0) == 0.0


def test_scalar_value_below_one():
    """For |v| < 10, order_magn = 0; round to 6 decimals."""
    assert _truncate_scalar(0.1234567890) == pytest.approx(0.123457)


def test_scalar_value_around_one():
    """For 1 <= |v| < 10, log10 in [0, 1) so ceil = 1; round to 5 decimals."""
    assert _truncate_scalar(1.234567890) == pytest.approx(1.23457)


def test_scalar_two_digits():
    """For 10 <= |v| < 100, order_magn = 2; round to 4 decimals."""
    assert _truncate_scalar(12.34567890) == pytest.approx(12.3457)


def test_scalar_exact_six_digit_integer():
    """For 100000 <= v < 1e6, order_magn = 6; round to 0 decimals."""
    assert _truncate_scalar(123456.789) == pytest.approx(123457.0)


def test_scalar_seven_digit_drops_to_negative_decimals():
    """For 1234567, order_magn = 7; round to -1 decimals -> nearest 10."""
    assert _truncate_scalar(1234567.0) == pytest.approx(1234570.0)


def test_scalar_negative_handled():
    assert _truncate_scalar(-1.234567890) == pytest.approx(-1.23457)
    assert _truncate_scalar(-1234567.0) == pytest.approx(-1234570.0)


def test_scalar_nan_passes_through():
    assert np.isnan(_truncate_scalar(float("nan")))


def test_scalar_inf_passes_through():
    assert np.isinf(_truncate_scalar(float("inf")))
    assert np.isinf(_truncate_scalar(float("-inf")))


def test_scalar_very_small_value():
    """log10 of small values is very negative; order_magn capped at 0;
    rounded to 6 decimals."""
    assert _truncate_scalar(0.0001234567) == pytest.approx(0.000123)


# --------------------------------------------------------------------------- #
# 1-D array
# --------------------------------------------------------------------------- #

def test_1d_array_each_value_processed():
    arr = np.array([0.123456789, 1.234567890, 12.345678])
    out = truncate_values(arr)
    np.testing.assert_allclose(
        out, [0.123457, 1.23457, 12.3457],
    )


def test_1d_columns_ignored():
    """`columns` is ignored for 1-D input."""
    arr = np.array([0.1234567, 1.234567, 12.34567])
    out = truncate_values(arr, columns=[0])  # ignored
    np.testing.assert_allclose(out, [0.123457, 1.23457, 12.3457])


def test_1d_empty_array():
    out = truncate_values(np.array([], dtype=float))
    assert out.size == 0


# --------------------------------------------------------------------------- #
# 2-D array
# --------------------------------------------------------------------------- #

def test_2d_default_processes_all_columns():
    arr = np.array([
        [0.1234567, 1.234567],
        [12.34567, 123.4567],
    ])
    out = truncate_values(arr)
    np.testing.assert_allclose(
        out,
        [[0.123457, 1.23457],
         [12.3457, 123.457]],
    )


def test_2d_column_subset_only_specified_processed():
    """Columns not in the list are passed through unchanged."""
    arr = np.array([
        [0.1234567890, 999999.999999],  # second col not processed
        [1.234567890, 999999.999999],
    ])
    out = truncate_values(arr, columns=[0])
    np.testing.assert_allclose(
        out[:, 0], [0.123457, 1.23457],
    )
    # Second column unchanged.
    np.testing.assert_array_equal(out[:, 1], arr[:, 1])


def test_2d_empty_array():
    out = truncate_values(np.empty((0, 3), dtype=float))
    assert out.shape == (0, 3)


# --------------------------------------------------------------------------- #
# Non-mutation
# --------------------------------------------------------------------------- #

def test_input_array_not_mutated():
    arr = np.array([1.234567890, 12.345678])
    snapshot = arr.copy()
    truncate_values(arr)
    np.testing.assert_array_equal(arr, snapshot)


# --------------------------------------------------------------------------- #
# Edge: 3-D rejected
# --------------------------------------------------------------------------- #

def test_3d_array_raises():
    arr = np.zeros((2, 2, 2))
    with pytest.raises(ValueError, match="1-D or 2-D"):
        truncate_values(arr)
