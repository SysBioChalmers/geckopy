"""Tests for adapter parameter schemas (validation)."""
import pytest
from pydantic import ValidationError

from geckopy.adapter.params import BayesianParams, SourceGroupRule


def test_defaults_are_consistent():
    # The shipped defaults must satisfy the group-key validator.
    bp = BayesianParams()
    assert set(bp.sigma0_log_source) == set(bp.source_groups)
    assert set(bp.shrink_thr_source) == set(bp.source_groups)
    assert set(bp.force_prior_thr_source) == set(bp.source_groups)


def test_mismatched_source_dict_raises():
    with pytest.raises(ValidationError, match="sigma0_log_source"):
        BayesianParams(
            source_groups={
                "dlkcat": SourceGroupRule(sources=["dlkcat"]),
                "brenda": SourceGroupRule(sources=["brenda"]),
            },
            sigma0_log_source={"dlkcat": 0.4},  # missing "brenda"
        )


def test_mismatched_schedule_lengths_raises():
    with pytest.raises(ValidationError, match="schedule_generations"):
        BayesianParams(
            schedule_generations=[1, 2, 9],
            schedule_samples=[1000, 800],
        )


def test_consistent_custom_groups_ok():
    bp = BayesianParams(
        source_groups={
            "a": SourceGroupRule(sources=["src_a"]),
            "b": SourceGroupRule(sources=["src_b"], match_okp=True),
        },
        sigma0_log_source={"a": 0.1, "b": 0.2},
        shrink_thr_source={"a": 1.0, "b": 2.0},
        force_prior_thr_source={"a": 0.0, "b": 1.0},
    )
    assert set(bp.source_groups) == {"a", "b"}
    assert bp.source_groups["b"].match_okp is True


def test_dropped_fields_are_rejected():
    """target_accept and variance_cap_* were confirmed dead/cosmetic
    in MATLAB and are not part of the schema; extra="forbid" should
    reject them outright rather than silently ignoring them."""
    with pytest.raises(ValidationError):
        BayesianParams(target_accept=10.0)
    with pytest.raises(ValidationError):
        BayesianParams(variance_cap_default=10.0)


def test_sparsity_threshold_default_is_not_the_confirmed_bad_range():
    bp = BayesianParams()
    assert bp.sparsity_threshold == pytest.approx(0.5)
