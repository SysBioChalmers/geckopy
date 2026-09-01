"""Tests for kcat_sensitivity_analysis.bayesian.transition.

Synthetic (X, w) arrays throughout -- no FBA, no EcModel.
"""
import numpy as np
import pandas as pd
import pytest
import scipy.stats

from geckopy.kcat_sensitivity_analysis.bayesian.transition import GeckoTransition


def _fitted(
    n_params: int = 2, n_particles: int = 200, seed: int = 0,
    sigma0_log: float = 0.4,
) -> tuple[GeckoTransition, pd.DataFrame, np.ndarray]:
    rng = np.random.default_rng(seed)
    log_x = rng.normal(loc=0.0, scale=sigma0_log, size=(n_particles, n_params))
    X = pd.DataFrame(
        np.exp(log_x), columns=[f"k{i}" for i in range(n_params)],
    )
    w = np.full(n_particles, 1.0 / n_particles)
    transition = GeckoTransition(np.full(n_params, sigma0_log))
    transition.fit(X, w)
    return transition, X, w


def test_fit_requires_matching_column_count():
    transition = GeckoTransition(np.array([0.4, 0.4, 0.4]))
    X = pd.DataFrame({"k0": [1.0, 2.0], "k1": [1.0, 2.0]})
    with pytest.raises(ValueError, match="columns"):
        transition.fit(X, np.array([0.5, 0.5]))


def test_rvs_single_before_fit_raises():
    transition = GeckoTransition(np.array([0.4]))
    with pytest.raises(RuntimeError, match="fit"):
        transition.rvs_single()


def test_bandwidth_blends_observed_std_with_prior():
    # All mass on a single point (zero observed std): with
    # adapt_frac_early=0.5, blended = 0.5*0 + 0.5*sigma0_log = 0.25,
    # which already exceeds the floor (0.15*sigma0_log = 0.075) -> the
    # floor shouldn't bind here, only the blend matters.
    X = pd.DataFrame({"k0": [2.0, 2.0, 2.0]})
    w = np.array([1.0, 1.0, 1.0])
    sigma0_log = np.array([0.5])
    transition = GeckoTransition(sigma0_log, adapt_frac_early=0.5, sigma_floor_frac=0.15)

    transition.fit(X, w)

    assert transition._log_bandwidth[0] == pytest.approx(0.25)


def test_bandwidth_floor_binds_when_blended_value_is_smaller():
    X = pd.DataFrame({"k0": [2.0, 2.0, 2.0]})  # zero observed std
    w = np.array([1.0, 1.0, 1.0])
    sigma0_log = np.array([0.5])
    # adapt_frac_early=1.0 -> blended = 1.0*std_obs + 0.0*sigma0_log = 0
    transition = GeckoTransition(sigma0_log, adapt_frac_early=1.0, sigma_floor_frac=0.15)

    transition.fit(X, w)

    assert transition._log_bandwidth[0] == pytest.approx(0.15 * 0.5)


def test_rvs_single_samples_are_positive_and_near_parents():
    transition, X, w = _fitted(n_params=2, n_particles=50)

    samples = [transition.rvs_single() for _ in range(200)]
    values = np.array([[s["k0"], s["k1"]] for s in samples])

    assert np.all(values > 0)
    # Samples should stay within a fairly generous band of the
    # fitted particles' own range (a perturbation kernel shouldn't
    # wander arbitrarily far in a single step).
    assert values.max() < X.to_numpy().max() * 20
    assert values.min() > X.to_numpy().min() / 20


def test_pdf_matches_hand_computed_mixture_for_series_input():
    X = pd.DataFrame({"k0": [1.0, 3.0]})
    w = np.array([0.25, 0.75])
    sigma0_log = np.array([0.3])
    transition = GeckoTransition(sigma0_log, adapt_frac_early=0.0, sigma_floor_frac=1.0)
    # adapt_frac_early=0, sigma_floor_frac=1.0 -> bandwidth == sigma0_log exactly.
    transition.fit(X, w)

    from pyabc.parameters import Parameter
    x = Parameter({"k0": 2.0})

    expected = (
        0.25 * scipy.stats.lognorm(s=0.3, scale=1.0).pdf(2.0)
        + 0.75 * scipy.stats.lognorm(s=0.3, scale=3.0).pdf(2.0)
    )
    assert transition.pdf(x) == pytest.approx(expected, rel=1e-6)


def test_pdf_dataframe_input_returns_array():
    transition, X, w = _fitted(n_params=1, n_particles=20)
    query = pd.DataFrame({"k0": [0.5, 1.0, 2.0]})

    densities = transition.pdf(query)

    assert isinstance(densities, np.ndarray)
    assert densities.shape == (3,)
    assert np.all(densities >= 0)


def test_component_logpdf_matches_lognorm_directly():
    transition, X, w = _fitted(n_params=2, n_particles=10)
    x = np.array([1.5, 0.8])
    parent = np.array([1.0, 1.0])

    result = transition.component_logpdf(x, parent)

    expected = sum(
        scipy.stats.lognorm.logpdf(xi, s=h, scale=pi)
        for xi, h, pi in zip(x, transition._log_bandwidth, parent)
    )
    assert result == pytest.approx(expected)


def test_component_logpdf_before_fit_raises():
    transition = GeckoTransition(np.array([0.4]))
    with pytest.raises(RuntimeError, match="fit"):
        transition.component_logpdf(np.array([1.0]), np.array([1.0]))


def test_rvs_single_draws_the_intended_lognormal_around_its_parent():
    """Each coordinate is its parent's value scaled by exp(bandwidth * Z).

    Fitting a single particle removes the parent-choice mixture, so the
    log-ratio to that parent is exactly N(0, log_bandwidth) per
    coordinate.
    """
    sigma0_log = np.array([0.25, 0.6])
    X = pd.DataFrame({"k0": [2.0], "k1": [5.0]})
    transition = GeckoTransition(
        sigma0_log, adapt_frac_early=0.0, sigma_floor_frac=1.0,
    )
    transition.fit(X, np.array([1.0]))

    np.random.seed(0)
    values = np.array(
        [[s["k0"], s["k1"]] for s in (transition.rvs_single() for _ in range(4000))]
    )
    log_ratio = np.log(values / X.to_numpy(dtype=float))

    assert np.all(values > 0)
    np.testing.assert_allclose(log_ratio.mean(axis=0), [0.0, 0.0], atol=0.03)
    np.testing.assert_allclose(log_ratio.std(axis=0), sigma0_log, rtol=0.05)
