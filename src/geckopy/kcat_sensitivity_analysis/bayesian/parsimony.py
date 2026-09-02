"""Trading fit against how many kcats were moved.

A tuned kcat vector that fits well by rewriting every parameter is not
a useful answer: changes should be few, and concentrated in the
parameters we are least confident about. Movement is measured here in
units of each kcat's own prior sigma -- MATLAB's ``devFromPrior`` --
which weights by confidence automatically, since ``sigma0_log`` is
smaller for trusted sources.

The tuner returns the best-RMSE particle, which on the full ecYeastGEM
model moves ~99% of kcats by ~8 prior sigmas. Reverting the least-moved
of those to their priors costs nothing up to about 3 sigma and beats
the shrunk "blend" estimate at every sparsity level, so a frontier is a
more honest report than either single vector. See
``docs/internal/matlab_replication_results.md``.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Sequence

import numpy as np


@dataclass(frozen=True)
class FrontierPoint:
    """One (parsimony, fit) trade-off, from :func:`parsimony_frontier`."""

    threshold: float
    kcat: np.ndarray
    n_changed: int
    mean_dev: float
    rmse: float


def movement_in_sigma(
    kcat: np.ndarray, kcat0: np.ndarray, sigma0_log: np.ndarray,
) -> np.ndarray:
    """Per-kcat displacement from the prior, in prior standard deviations.

    ``|log(kcat / kcat0)| / sigma0_log``, matching MATLAB's
    ``devFromPrior = abs(muLog - log(kcat0)) ./ sigma0log``.
    """
    return np.abs(np.log(np.asarray(kcat, dtype=float) / kcat0)) / sigma0_log


def revert_below(
    kcat: np.ndarray,
    kcat0: np.ndarray,
    sigma0_log: np.ndarray,
    threshold: float,
) -> np.ndarray:
    """Put every kcat that moved less than ``threshold`` sigma back on
    its prior, leaving the rest at their tuned values."""
    keep = movement_in_sigma(kcat, kcat0, sigma0_log) >= threshold
    return np.where(keep, np.asarray(kcat, dtype=float), kcat0)


def n_changed(
    kcat: np.ndarray, kcat0: np.ndarray, rel_tol: float = 0.02,
) -> int:
    """How many kcats differ from their prior by more than ``rel_tol``.

    The 2% default matches the threshold MATLAB reports its per-source
    "unchanged" percentages at.
    """
    return int(
        np.count_nonzero(
            np.abs(np.log(np.asarray(kcat, dtype=float) / kcat0))
            > np.log1p(rel_tol)
        )
    )


def parsimony_frontier(
    score: Callable[[np.ndarray], float],
    kcat: np.ndarray,
    kcat0: np.ndarray,
    sigma0_log: np.ndarray,
    thresholds: Sequence[float] = (0.0, 1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0, 15.0),
) -> list[FrontierPoint]:
    """Score ``kcat`` with progressively more of it reverted to prior.

    Parameters
    ----------
    score
        Scores one kcat vector; one call per threshold.
    kcat, kcat0, sigma0_log
        Tuned vector, prior vector, and per-kcat prior sigma.
    thresholds
        Sigma cutoffs to sweep, ascending. ``0.0`` reverts nothing and
        so reproduces ``kcat``'s own score.

    Returns
    -------
    list of FrontierPoint, in the order given.
    """
    kcat0 = np.asarray(kcat0, dtype=float)
    sigma0_log = np.asarray(sigma0_log, dtype=float)
    points = []
    for t in thresholds:
        vec = revert_below(kcat, kcat0, sigma0_log, t)
        points.append(FrontierPoint(
            threshold=float(t),
            kcat=vec,
            n_changed=n_changed(vec, kcat0),
            mean_dev=float(movement_in_sigma(vec, kcat0, sigma0_log).mean()),
            rmse=float(score(vec)),
        ))
    return points


def best_parsimonious(
    points: Sequence[FrontierPoint], tolerance: float = 0.02,
) -> FrontierPoint:
    """The fewest-changes point whose RMSE is within ``tolerance`` of the
    best on the frontier.

    Defaults to 2%, which is the scale of the alternate-optimum noise
    measured on this model -- below that, two vectors are not
    meaningfully different in fit, so the sparser one wins.
    """
    if not points:
        raise ValueError("points is empty.")
    best_rmse = min(p.rmse for p in points)
    within = [p for p in points if p.rmse <= best_rmse * (1.0 + tolerance)]
    return min(within, key=lambda p: p.n_changed)
