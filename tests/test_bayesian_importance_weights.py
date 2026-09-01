"""Tests for kcat_sensitivity_analysis.bayesian.importance_weights
(Axis 2, variant B).

Uses plain multivariate-normal prior/transition callables (not tied to
the kcat-specific log-space meaning) since this module's formula is
generic; the kcat-specific prior/transition pdfs are `priors.py`'s and
`transition.py`'s concern, wired together in `tuning.py`.
"""
import numpy as np
import pytest

from geckopy.kcat_sensitivity_analysis.bayesian.importance_weights import (
    compute_importance_weights,
)


def _normal_logpdf(x: np.ndarray, mean: np.ndarray, std: float) -> float:
    d = x.size
    return float(
        -0.5 * d * np.log(2 * np.pi * std**2) - 0.5 * np.sum(((x - mean) / std) ** 2)
    )


def _prior_logpdf(x: np.ndarray) -> float:
    return _normal_logpdf(x, mean=np.zeros_like(x), std=1.0)


def _make_transition_logpdf(sigma_k: float):
    def transition_logpdf(x: np.ndarray, parent: np.ndarray) -> float:
        return _normal_logpdf(x, mean=parent, std=sigma_k)
    return transition_logpdf


def test_generation_one_returns_uniform_weights():
    particles = np.array([[0.1, 0.2, 0.3], [1.0, 1.1, 1.2]])  # (n_params=2, n_particles=3)

    weights = compute_importance_weights(particles, _prior_logpdf)

    assert weights == pytest.approx([1 / 3, 1 / 3, 1 / 3])


def test_single_parent_matches_prior_over_transition_ratio():
    parent = np.array([0.5, 0.5])
    parents = parent.reshape(-1, 1)
    parent_weights = np.array([1.0])
    sigma_k = 0.3
    transition_logpdf = _make_transition_logpdf(sigma_k)

    particles = np.array([[0.0, 0.5, 1.0], [0.0, 0.5, 1.0]])  # 3 particles, 2 params

    weights = compute_importance_weights(
        particles, _prior_logpdf,
        parents=parents, parent_weights=parent_weights,
        transition_logpdf=transition_logpdf,
    )

    # With a single parent (weight 1), the denominator is exactly that
    # parent's transition density -- no mixture ambiguity, so the
    # un-normalised weight is exactly prior/transition per particle.
    raw = np.array([
        np.exp(_prior_logpdf(particles[:, i]) - transition_logpdf(particles[:, i], parent))
        for i in range(3)
    ])
    expected = raw / raw.sum()
    assert weights == pytest.approx(expected)
    assert weights.sum() == pytest.approx(1.0)


def test_multiple_parents_denominator_is_weighted_mixture():
    parents = np.array([[0.0, 2.0], [0.0, 2.0]])  # 2 parents, 2 params
    parent_weights = np.array([0.3, 0.7])
    sigma_k = 0.5
    transition_logpdf = _make_transition_logpdf(sigma_k)
    particles = np.array([[0.1, 1.9], [0.1, 1.9]])  # 2 particles

    weights = compute_importance_weights(
        particles, _prior_logpdf,
        parents=parents, parent_weights=parent_weights,
        transition_logpdf=transition_logpdf,
    )

    expected_raw = []
    for i in range(2):
        theta = particles[:, i]
        denom = sum(
            parent_weights[j] * np.exp(transition_logpdf(theta, parents[:, j]))
            for j in range(2)
        )
        expected_raw.append(np.exp(_prior_logpdf(theta)) / denom)
    expected_raw = np.array(expected_raw)
    expected = expected_raw / expected_raw.sum()

    assert weights == pytest.approx(expected, rel=1e-6)


def test_zero_weight_parent_is_ignored():
    parents_with_dead = np.array([[0.0, 5.0], [0.0, 5.0]])
    weights_with_dead = np.array([1.0, 0.0])
    parents_without = np.array([[0.0], [0.0]])
    weights_without = np.array([1.0])
    sigma_k = 0.4
    transition_logpdf = _make_transition_logpdf(sigma_k)
    particles = np.array([[0.2, -0.3], [0.2, -0.3]])

    w1 = compute_importance_weights(
        particles, _prior_logpdf,
        parents=parents_with_dead, parent_weights=weights_with_dead,
        transition_logpdf=transition_logpdf,
    )
    w2 = compute_importance_weights(
        particles, _prior_logpdf,
        parents=parents_without, parent_weights=weights_without,
        transition_logpdf=transition_logpdf,
    )

    assert w1 == pytest.approx(w2)


def test_rejects_missing_parent_weights_or_transition():
    particles = np.array([[0.0, 1.0]])
    parents = np.array([[0.0, 1.0]])
    with pytest.raises(ValueError, match="parent_weights"):
        compute_importance_weights(particles, _prior_logpdf, parents=parents)


def test_rejects_shape_mismatch():
    particles = np.array([[0.0, 1.0], [0.0, 1.0]])  # 2 params
    parents = np.array([[0.0, 1.0]])  # 1 param -- mismatch
    with pytest.raises(ValueError, match="parents"):
        compute_importance_weights(
            particles, _prior_logpdf,
            parents=parents, parent_weights=np.array([0.5, 0.5]),
            transition_logpdf=_make_transition_logpdf(0.3),
        )
