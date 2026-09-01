"""Per-generation, per-source-group diagnostics.

Ported design intent from GECKO MATLAB's per-source diagnostic block
in ``bayesianSensitivityTuning.m`` (near-prior counts, mean deviation,
variance ratio, "active" counts) -- but computed directly from a
generation's accepted particle population (``kcat_top`` + ``weights``)
rather than from ``posterior.PosteriorUpdate``'s
``shrink_weight``/``sigma_log`` outputs. Nothing here reads
``PosteriorUpdate``, so these metrics stay comparable across selection
variants and remain meaningful even though the shrinkage blend never
feeds back into sampling.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass
class GroupDiagnostics:
    """One source-group's diagnostics for one generation.

    Attributes
    ----------
    group
        Source-group name (or :data:`~.priors.UNLABELLED_GROUP`).
    n
        Number of tunable parameters in this group.
    mean_deviation
        Mean, over this group's parameters, of ``|weighted mean
        log-kcat - log(kcat0)| / sigma0_log`` -- how many prior-sigmas
        the population's centre of mass has moved.
    variance_ratio
        Mean, over this group's parameters, of ``weighted std(log-kcat)
        / sigma0_log`` -- how the population's spread compares to the
        prior's (1.0 = unchanged, <1 = contracted, >1 = expanded).
    frac_active
        Weighted fraction of (parameter, particle) pairs whose
        per-particle deviation from ``kcat0`` exceeds
        ``active_threshold`` prior-sigmas.
    frac_near_prior
        Fraction of this group's parameters whose spread stayed within
        ``near_prior_tol`` (absolute, log-space) of the prior's.
    """

    group: str
    n: int
    mean_deviation: float
    variance_ratio: float
    frac_active: float
    frac_near_prior: float


@dataclass
class GenerationDiagnostics:
    """One generation's full diagnostic snapshot."""

    generation: int
    n_total: int
    n_accepted: int
    acceptance_rate: float
    best_rmse: float
    mean_rmse: float
    median_rmse: float
    by_group: dict[str, GroupDiagnostics] = field(default_factory=dict)


def compute_generation_diagnostics(
    generation: int,
    kcat_top: np.ndarray,
    weights: np.ndarray,
    rmse_top: np.ndarray,
    n_total_this_gen: int,
    kcat0: np.ndarray,
    sigma0_log: np.ndarray,
    groups: np.ndarray,
    *,
    active_threshold: float = 0.3,
    near_prior_tol: float = 0.1,
) -> GenerationDiagnostics:
    """Diagnostics for one generation's accepted population.

    Parameters
    ----------
    generation
        Generation number (1-indexed).
    kcat_top
        This generation's accepted particles, shape ``(n_params,
        n_accepted)``.
    weights
        This generation's particle weights (need not be pre-normalised),
        shape ``(n_accepted,)``.
    rmse_top
        This generation's accepted particles' distances, shape
        ``(n_accepted,)``.
    n_total_this_gen
        Total candidates evaluated this generation (before selection)
        -- for the acceptance rate.
    kcat0, sigma0_log, groups
        As in ``posterior.update_posterior_shrinkage``.
    active_threshold
        Per-particle log-space deviation (in prior-sigma units) above
        which a (parameter, particle) pair counts as "active".
    near_prior_tol
        Absolute log-space tolerance for "near prior" spread.

    Returns
    -------
    GenerationDiagnostics
    """
    n_params = kcat_top.shape[0]
    if kcat0.shape != (n_params,) or sigma0_log.shape != (n_params,) or groups.shape != (n_params,):
        raise ValueError(
            "kcat0, sigma0_log, and groups must each have shape "
            f"({n_params},) to match kcat_top's first dimension."
        )

    w_norm = np.asarray(weights, dtype=float)
    w_norm = w_norm / w_norm.sum()

    log_kcat_top = np.log(kcat_top)
    log_kcat0 = np.log(kcat0)

    mean_log = log_kcat_top @ w_norm
    var_log = ((log_kcat_top - mean_log[:, None]) ** 2) @ w_norm
    std_log = np.sqrt(var_log)

    deviation = np.abs(mean_log - log_kcat0) / sigma0_log
    variance_ratio = std_log / sigma0_log

    per_particle_dev = np.abs(log_kcat_top - log_kcat0[:, None]) / sigma0_log[:, None]
    active_frac_per_param = (per_particle_dev > active_threshold).astype(float) @ w_norm
    near_prior_per_param = np.abs(std_log - sigma0_log) < near_prior_tol

    by_group: dict[str, GroupDiagnostics] = {}
    for name in np.unique(groups):
        mask = groups == name
        by_group[str(name)] = GroupDiagnostics(
            group=str(name),
            n=int(mask.sum()),
            mean_deviation=float(np.mean(deviation[mask])),
            variance_ratio=float(np.mean(variance_ratio[mask])),
            frac_active=float(np.mean(active_frac_per_param[mask])),
            frac_near_prior=float(np.mean(near_prior_per_param[mask])),
        )

    return GenerationDiagnostics(
        generation=generation,
        n_total=n_total_this_gen,
        n_accepted=len(rmse_top),
        acceptance_rate=len(rmse_top) / n_total_this_gen,
        best_rmse=float(np.min(rmse_top)),
        mean_rmse=float(np.mean(rmse_top)),
        median_rmse=float(np.median(rmse_top)),
        by_group=by_group,
    )
