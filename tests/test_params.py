"""Tests for adapter parameter schemas (validation)."""
import pytest
from pydantic import ValidationError

from geckopy.adapter.params import BayesianParams


def test_defaults_are_consistent():
    # The shipped defaults must satisfy the parallel-list validator.
    bp = BayesianParams()
    assert len(bp.sigma0_log_source) == len(bp.kcat_sources)


def test_mismatched_source_list_raises():
    with pytest.raises(ValidationError, match="sigma0_log_source"):
        BayesianParams(
            kcat_sources=["dlkcat", "brenda", "custom"],
            sigma0_log_source=[0.4, 0.2],  # one short
        )


def test_mismatched_schedule_lengths_raises():
    with pytest.raises(ValidationError, match="schedule_generations"):
        BayesianParams(
            schedule_generations=[1, 2, 9],
            schedule_samples=[1000, 800],
        )


def test_consistent_custom_lists_ok():
    bp = BayesianParams(
        kcat_sources=["a", "b"],
        sigma0_log_source=[0.1, 0.2],
        shrink_thr_source=[1.0, 2.0],
        variance_cap_source=[5.0, 6.0],
        force_prior_thr_source=[0.0, 1.0],
    )
    assert bp.kcat_sources == ["a", "b"]
