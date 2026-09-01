"""Axis 1: which particles survive each ABC-SMC generation.

Two variants, per the plan's comparison design -- both take a batch of
already-simulated distances and return which indices survive plus the
generation's effective epsilon, so ``tuning.py`` can swap between them
behind one interface:

- :func:`truncation_select` -- MATLAB-faithful. Ported from
  ``bayesianSensitivityTuning.m``'s acceptance step (fixed batch,
  retrospective top-``min_keep``-fraction cut by distance, uniform
  weights afterward), verified against ``develop4``'s current source
  -- but *without* ``targetAccept``'s percentile pre-cut: confirmed
  dead for the shipped defaults (``minKeep``'s floor consistently
  bound first) and dropped from ``BayesianParams`` entirely, so this
  implements the mechanism that actually runs, not the vestigial one.
- :func:`quantile_epsilon_select` -- pyABC-native. Accepts against a
  *fixed, prospective* epsilon (computed from the *previous*
  generation via :func:`next_quantile_epsilon`, which reuses pyABC's
  own ``weighted_quantile`` -- the same computation
  ``pyabc.epsilon.QuantileEpsilon``/``MedianEpsilon`` perform
  internally, without needing to drive their ``ABCSMC``-coupled
  lifecycle). Unlike truncation's fixed accepted-count, however many
  candidates pass is however many are returned -- pyABC's actual
  "streaming, unbounded simulation budget" semantics is a *calling
  pattern* (keep simulating until enough pass), which is
  ``tuning.py``'s responsibility since it owns the simulator calls;
  this function is the accept/reject test that pattern uses.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np
from pyabc.weighted_statistics import weighted_quantile


@dataclass
class SelectionResult:
    """One generation's selection outcome.

    Attributes
    ----------
    accepted_idx
        Indices into the candidate ``distances`` array that survived
        (not guaranteed sorted).
    epsilon
        The effective distance threshold for this generation: the
        worst accepted distance, for :func:`truncation_select`; the
        fixed threshold that was applied, for
        :func:`quantile_epsilon_select`.
    """

    accepted_idx: np.ndarray
    epsilon: float


def truncation_select(distances: np.ndarray, *, min_keep: float) -> SelectionResult:
    """MATLAB-faithful top-``min_keep``-fraction truncation.

    Parameters
    ----------
    distances
        This generation's candidate distances (lower is better),
        combining new proposals with the previous generation's
        accepted set, per MATLAB's SMC step (``rmse = [newRmse;
        rmseTop]``).
    min_keep
        Fraction of ``distances`` to keep, in ``(0, 1]``. At least one
        candidate is always kept.

    Returns
    -------
    SelectionResult
    """
    if not 0.0 < min_keep <= 1.0:
        raise ValueError(f"min_keep must be in (0, 1]; got {min_keep}.")
    n = len(distances)
    if n == 0:
        raise ValueError("distances is empty.")
    keep_n = max(1, int(np.floor(min_keep * n)))
    order = np.argsort(distances, kind="stable")
    accepted_idx = order[:keep_n]
    epsilon = float(distances[accepted_idx[-1]])
    return SelectionResult(accepted_idx=accepted_idx, epsilon=epsilon)


def quantile_epsilon_select(
    distances: np.ndarray, *, epsilon: float,
) -> SelectionResult:
    """Accept every candidate within a pre-computed, fixed epsilon.

    Parameters
    ----------
    distances
        This generation's candidate distances.
    epsilon
        The threshold to apply, from :func:`next_quantile_epsilon`
        applied to the *previous* generation's accepted distances (or
        an initial value for generation 1). Unlike
        :func:`truncation_select`, this is fixed *before* seeing
        ``distances`` -- the "prospective" half of pyABC's scheme.

    Returns
    -------
    SelectionResult
        ``accepted_idx`` may be empty (if nothing in this batch beat
        ``epsilon``) or include every candidate -- the accepted count
        isn't forced to any particular size, unlike
        :func:`truncation_select`.
    """
    accepted_idx = np.flatnonzero(np.asarray(distances) <= epsilon)
    return SelectionResult(accepted_idx=accepted_idx, epsilon=float(epsilon))


def next_quantile_epsilon(
    accepted_distances: np.ndarray,
    *,
    alpha: float = 0.5,
    weights: Optional[np.ndarray] = None,
) -> float:
    """Adaptive epsilon for the next generation (pyABC's scheme).

    Reuses ``pyabc.weighted_statistics.weighted_quantile`` directly --
    the same computation ``pyabc.epsilon.QuantileEpsilon``/
    ``MedianEpsilon`` perform internally as part of a full
    ``ABCSMC.run()``, without needing to drive that class's
    ``ABCSMC``-coupled ``initialize``/``update`` lifecycle.

    Parameters
    ----------
    accepted_distances
        The current generation's accepted distances.
    alpha
        The quantile to use; ``0.5`` (the default) reproduces
        ``MedianEpsilon``.
    weights
        Optional per-particle weights (e.g. from
        ``importance_weights.py`` under Axis 2 variant B). ``None``
        weights every accepted particle equally.

    Returns
    -------
    float
    """
    if len(accepted_distances) == 0:
        raise ValueError("accepted_distances is empty.")
    if weights is not None:
        weights = np.asarray(weights, dtype=float)
        weights = weights / weights.sum()  # weighted_quantile asserts normalisation
    return float(
        weighted_quantile(
            np.asarray(accepted_distances, dtype=float), weights, alpha=alpha,
        )
    )
