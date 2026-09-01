"""Axis 2, variant A: MATLAB-faithful shrink-weight/force-to-prior/
sparsity-snap posterior update.

Ported from GECKO MATLAB:
src/geckomat/kcat_sensitivity_analysis/Bayesian/bayesianSensitivityTuning.m
(the per-generation "Update posterior kcat and sigmaLog" block).
Verified against ``develop4``'s current source.

Pure numpy -- zero dependency on ``EcModel``/``pyabc``, independently
unit-testable (see ``tests/test_bayesian_posterior.py``).

MATLAB's per-source variance cap (``varianceCapDefault``/
``varianceCapSource``) is deliberately *not* ported here: tracing its
actual data flow shows the capped ``kcatSigmaLog`` feeds only
diagnostics/reporting -- ``buildLowRankLogProposal`` (which sets the
*actual* next-generation proposal width) recomputes its own std
directly from the raw accepted samples, never reading the capped
value. Confirmed cosmetic, matching the plan's resolved decision to
drop ``variance_cap_*`` from ``BayesianParams`` outright.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

#: Must match `priors.UNLABELLED_GROUP`. Duplicated as a plain literal
#: (rather than importing from `.priors`) so this module stays free of
#: even a transitive `pyabc` import -- `priors.py` requires pyabc at
#: module level, but this module's own logic never touches it.
_UNLABELLED_GROUP = "unlabelled"


@dataclass
class PosteriorUpdate:
    """One generation's Axis-2-variant-A posterior update.

    Attributes
    ----------
    kcats
        Updated point-estimate kcat values, one per tunable row.
    sigma_log
        Updated per-row log-space std dev. For reporting/diagnostics
        only -- MATLAB never feeds this back into the next
        generation's proposal width (see module docstring), so
        neither does this port.
    shrink_weight
        Per-row blend weight actually applied (0 = fully prior, 1 =
        fully the accepted-sample posterior mean). Exposed for
        diagnostics.
    forced_to_prior
        Per-row boolean: True where a force-prior threshold snapped
        the shrink weight to 0 before blending.
    snapped_to_prior
        Per-row boolean: True where the post-hoc sparsity threshold
        snapped the blended estimate back to the *exact* prior value.
    """

    kcats: np.ndarray
    sigma_log: np.ndarray
    shrink_weight: np.ndarray
    forced_to_prior: np.ndarray
    snapped_to_prior: np.ndarray


def update_posterior_shrinkage(
    kcat_top: np.ndarray,
    kcat0: np.ndarray,
    sigma0_log: np.ndarray,
    groups: np.ndarray,
    *,
    shrink_thr_default: float,
    shrink_thr_source: dict[str, float],
    force_prior_thr_default: float,
    force_prior_thr_source: dict[str, float],
    sparsity_threshold: float,
) -> PosteriorUpdate:
    """Blend one generation's accepted samples toward the prior.

    Parameters
    ----------
    kcat_top
        Accepted particles this generation, shape ``(n_params,
        n_accepted)`` -- MATLAB's ``kcatTop``.
    kcat0
        Current (pre-generation) point estimate per parameter, shape
        ``(n_params,)``.
    sigma0_log
        Per-parameter prior std dev in log-space (from
        ``priors.build_sigma0_log``), shape ``(n_params,)``.
    groups
        Per-parameter source-group name (from
        ``priors.classify_kcat_sources``), shape ``(n_params,)``.
        Entries equal to ``"unlabelled"`` (``priors.UNLABELLED_GROUP``)
        use the ``*_default`` thresholds.
    shrink_thr_default, shrink_thr_source, force_prior_thr_default,
    force_prior_thr_source, sparsity_threshold
        From ``BayesianParams``.

    Returns
    -------
    PosteriorUpdate
    """
    n_params, n_accepted = kcat_top.shape
    if kcat0.shape != (n_params,):
        raise ValueError(
            f"kcat0 shape {kcat0.shape} must be ({n_params},) to match "
            f"kcat_top's first dimension."
        )
    if sigma0_log.shape != (n_params,):
        raise ValueError(
            f"sigma0_log shape {sigma0_log.shape} must be ({n_params},)."
        )
    if groups.shape != (n_params,):
        raise ValueError(f"groups shape {groups.shape} must be ({n_params},).")
    if n_accepted == 0:
        raise ValueError("kcat_top has no accepted particles (dimension 1 is 0).")

    log_kcat_top = np.log(kcat_top)
    mu_log = log_kcat_top.mean(axis=1)
    # MATLAB's std(X, 1, dim): normalise by N, not N-1.
    kcat_sigma_log = log_kcat_top.std(axis=1, ddof=0)

    log_kcat0 = np.log(kcat0)
    dev_from_prior = np.abs(mu_log - log_kcat0) / sigma0_log

    shrink_weight = np.minimum(dev_from_prior / shrink_thr_default, 1.0)
    for name, thr in shrink_thr_source.items():
        mask = groups == name
        shrink_weight[mask] = np.minimum(dev_from_prior[mask] / thr, 1.0)

    forced_to_prior = np.zeros(n_params, dtype=bool)
    for name, thr in force_prior_thr_source.items():
        if thr <= 0:  # "only apply if threshold is positive" (MATLAB)
            continue
        mask = groups == name
        forced_to_prior[mask] = dev_from_prior[mask] < thr
    unlabelled_mask = groups == _UNLABELLED_GROUP
    if force_prior_thr_default > 0:
        forced_to_prior[unlabelled_mask] = (
            dev_from_prior[unlabelled_mask] < force_prior_thr_default
        )

    shrink_weight = np.where(forced_to_prior, 0.0, shrink_weight)

    sigma_log_raw = shrink_weight * kcat_sigma_log + (1 - shrink_weight) * sigma0_log
    kcats_raw = np.exp(shrink_weight * mu_log + (1 - shrink_weight) * log_kcat0)

    snapped_to_prior = (
        np.abs(np.log(kcats_raw) - log_kcat0) < sparsity_threshold * sigma0_log
    )
    kcats = np.where(snapped_to_prior, kcat0, kcats_raw)
    sigma_log = np.where(snapped_to_prior, sigma0_log, sigma_log_raw)

    return PosteriorUpdate(
        kcats=kcats,
        sigma_log=sigma_log,
        shrink_weight=shrink_weight,
        forced_to_prior=forced_to_prior,
        snapped_to_prior=snapped_to_prior,
    )
