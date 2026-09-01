"""Axis 2, variant B: proper SMC-ABC importance weighting.

Not a MATLAB port -- MATLAB's tuning loop never computes real particle
weights (its selection step just keeps a truncated set with implicit
uniform weight, see ``selection.truncation_select``). This module
implements the standard sequential-importance-sampling particle weight
update instead (Toni et al. 2009; Beaumont et al. 2009), the
theoretically proper alternative to MATLAB's hand-tuned shrink-weight
blend (``posterior.py``): with a tight per-group prior
(``sigma0_log_source``), shrinkage-to-prior falls out automatically
from the weighting, no hand-tuned thresholds needed.

For generation t>1, particle i (drawn by perturbing a parent particle
from generation t-1 through a transition kernel) gets weight

    w_i^(t) ∝ prior_pdf(theta_i^(t))
              / sum_j( parent_weight_j^(t-1) * transition_pdf(theta_i^(t) | parent_j^(t-1)) )

Generation 1 has no parents/transition kernel -- particles are drawn
directly from the prior, so the standard SMC-ABC update degenerates to
uniform weights.

Pure numpy plus caller-supplied pdf callables -- no direct ``pyabc``
dependency (even though those callables typically wrap a
``pyabc.Distribution``/``GeckoTransition`` in practice), matching
``posterior.py``'s testability.
"""
from __future__ import annotations

from typing import Callable, Optional

import numpy as np


def compute_importance_weights(
    particles: np.ndarray,
    prior_logpdf: Callable[[np.ndarray], float],
    *,
    parents: Optional[np.ndarray] = None,
    parent_weights: Optional[np.ndarray] = None,
    transition_logpdf: Optional[Callable[[np.ndarray, np.ndarray], float]] = None,
) -> np.ndarray:
    """Normalised SMC-ABC importance weights for one generation.

    Parameters
    ----------
    particles
        This generation's accepted particles, shape ``(n_params,
        n_particles)``.
    prior_logpdf
        Callable taking one particle (shape ``(n_params,)``) and
        returning its log prior density -- e.g. a
        ``pyabc.Distribution`` wrapped to take a plain array and
        return ``log(pdf(...))``.
    parents
        Previous generation's particles, shape ``(n_params,
        n_parents)``. ``None`` for generation 1.
    parent_weights
        Previous generation's *normalised* weights, shape
        ``(n_parents,)``. Required iff ``parents`` is given.
    transition_logpdf
        Callable taking ``(particle, parent)`` and returning the log
        transition density of proposing ``particle`` from ``parent``
        -- e.g. a fitted ``transition.GeckoTransition.pdf`` wrapped
        the same way. Required iff ``parents`` is given.

    Returns
    -------
    numpy.ndarray, shape ``(n_particles,)``
        Normalised weights (sum to 1).

    Raises
    ------
    ValueError
        On shape mismatches, missing required arguments, or if every
        particle ends up with zero weight (numerically degenerate).
    """
    n_params, n_particles = particles.shape

    if parents is None:
        return np.full(n_particles, 1.0 / n_particles)

    if parent_weights is None or transition_logpdf is None:
        raise ValueError(
            "parent_weights and transition_logpdf are required when "
            "parents is given."
        )
    n_parents = parents.shape[1]
    if parents.shape[0] != n_params:
        raise ValueError(
            f"parents has {parents.shape[0]} parameter rows; expected "
            f"{n_params} to match particles."
        )
    if parent_weights.shape != (n_parents,):
        raise ValueError(
            f"parent_weights shape {parent_weights.shape} must be "
            f"({n_parents},) to match parents' second dimension."
        )

    log_weights = np.empty(n_particles, dtype=float)
    log_parent_weights = np.full(n_parents, -np.inf)
    positive = parent_weights > 0
    log_parent_weights[positive] = np.log(parent_weights[positive])

    for i in range(n_particles):
        theta = particles[:, i]
        log_prior = prior_logpdf(theta)

        log_terms = np.array([
            log_parent_weights[j] + transition_logpdf(theta, parents[:, j])
            for j in range(n_parents)
            if positive[j]
        ])
        if log_terms.size == 0 or np.all(np.isneginf(log_terms)):
            log_weights[i] = -np.inf
            continue
        m = np.max(log_terms)
        log_denom = m + np.log(np.sum(np.exp(log_terms - m)))
        log_weights[i] = log_prior - log_denom

    finite = np.isfinite(log_weights)
    if not finite.any():
        raise ValueError(
            "Every particle got zero importance weight (all -inf) -- "
            "check that the transition kernel can actually reach these "
            "particles from the given parents."
        )
    m = np.max(log_weights[finite])
    weights = np.where(finite, np.exp(log_weights - m), 0.0)
    return weights / weights.sum()
