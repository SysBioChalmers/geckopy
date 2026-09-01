"""Per-kcat lognormal priors for Bayesian kcat tuning.

Ported from GECKO MATLAB:
src/geckomat/kcat_sensitivity_analysis/Bayesian/bayesianSensitivityTuning.m
(the ``kcatSourceIdx``/``sigma0log`` construction and generation 1's
inline ``lognrnd`` sampling). Verified against ``develop4``'s current
source.

Two prior builders are provided, matching the plan's Axis 2 (how
trust-tier priors constrain movement):

- :func:`build_kcat_prior` -- the plain per-source lognormal prior,
  shared by both regularization variants (variant A's shrink-weight
  blend consumes it as a fixed reference point; variant B's importance
  weighting consumes it as the actual prior density).
- :func:`build_kcat_sparsity_prior` -- variant B's natural complement:
  a genuine sparsity-inducing spike-and-slab mixture per parameter,
  so an insignificant parameter's *posterior* concentrates at the
  prior on its own, rather than needing MATLAB's post-hoc "snap
  negligible changes to prior" step.
"""
from __future__ import annotations

from typing import Optional

import numpy as np
import pyabc
import scipy.stats
from pyabc.random_variables.random_variables import RVBase

from ...adapter.params import BayesianParams

#: Sentinel group name for an ``ec.source`` value matched by no
#: ``source_groups`` entry -- MATLAB's ``noKcatSource``.
UNLABELLED_GROUP = "unlabelled"


def classify_kcat_source(
    source: str,
    params: BayesianParams,
    *,
    okp_method: Optional[str] = None,
) -> str:
    """Map one ``ec.source`` string to a trust-tier group name.

    Matching is case-insensitive, mirroring MATLAB's
    ``strcmpi(ecModel.ec.source, kcatSources{i})``.

    Parameters
    ----------
    source
        One ``ec.source`` entry.
    params
        Supplies ``source_groups``.
    okp_method
        The project's configured OpenKineticsPredictor method
        (``OkpParams.method``), if any. A group with ``match_okp=True``
        also matches ``source`` when it equals this (case-insensitive).

    Returns
    -------
    str
        A key of ``params.source_groups``, or :data:`UNLABELLED_GROUP`
        if nothing matches -- callers then fall back to
        ``params.sigma0_log_default`` / ``shrink_thr_default`` /
        ``force_prior_thr_default``.
    """
    source_lower = source.lower()
    for name, rule in params.source_groups.items():
        if any(source_lower == s.lower() for s in rule.sources):
            return name
        if (
            rule.match_okp
            and okp_method is not None
            and source_lower == okp_method.lower()
        ):
            return name
    return UNLABELLED_GROUP


def classify_kcat_sources(
    ec_sources: list[str],
    params: BayesianParams,
    *,
    okp_method: Optional[str] = None,
) -> np.ndarray:
    """Vectorised :func:`classify_kcat_source` over ``ec.source``.

    Returns
    -------
    numpy.ndarray of str, shape ``(len(ec_sources),)``
    """
    return np.array(
        [
            classify_kcat_source(s, params, okp_method=okp_method)
            for s in ec_sources
        ],
        dtype=object,
    )


def build_sigma0_log(groups: np.ndarray, params: BayesianParams) -> np.ndarray:
    """Per-row prior std dev in log-space, from source-group membership.

    Mirrors MATLAB's ``sigma0log = sigma0logDefault * ones(...);
    sigma0log(uniqKcatParams) = sigma0logSource(kcatSourceIdx(...))``:
    every row starts at ``params.sigma0_log_default``, then rows whose
    group is a real (non-:data:`UNLABELLED_GROUP`) entry get
    overridden by that group's ``sigma0_log_source`` value.

    Parameters
    ----------
    groups
        Per-row group names from :func:`classify_kcat_sources`.
    params
        Supplies ``sigma0_log_default`` and ``sigma0_log_source``.

    Returns
    -------
    numpy.ndarray of float, same shape as ``groups``.
    """
    sigma0_log = np.full(len(groups), params.sigma0_log_default, dtype=float)
    for name, value in params.sigma0_log_source.items():
        sigma0_log[groups == name] = value
    return sigma0_log


def _check_shapes(kcat0: np.ndarray, sigma0_log: np.ndarray) -> None:
    if kcat0.shape != sigma0_log.shape:
        raise ValueError(
            f"kcat0 and sigma0_log must have the same shape; got "
            f"{kcat0.shape} vs {sigma0_log.shape}."
        )
    if np.any(kcat0 <= 0):
        raise ValueError(
            "kcat0 must be strictly positive for every entry -- filter "
            "out unassigned (kcat <= 0) rows before building a prior; "
            "they aren't tunable parameters."
        )


def _lognorm_scale_at(k0: float, sigma_log: float) -> float:
    """The ``scale`` (``exp(mu)``) that puts a lognormal's *mean* --
    not its median -- at ``k0``.

    Matches MATLAB's ``muLog = log(kcats) - 0.5 .* sigma0log.^2`` bias
    correction before ``lognrnd(muLog, sigma0log)``: for a lognormal
    with underlying-normal mean ``mu`` and std ``sigma_log``, the
    *median* is ``exp(mu)`` but the *mean* is ``exp(mu + sigma_log**2
    / 2)``. Subtracting ``0.5 * sigma_log**2`` from ``mu = log(k0)``
    up front cancels that offset, so the constructed lognormal's
    *mean* -- not its median -- sits at ``k0``: without the
    correction, the mean would drift above ``k0`` by a
    sigma-dependent factor.
    """
    mu = np.log(k0) - 0.5 * sigma_log**2
    return float(np.exp(mu))


def _lognorm_mean_at(k0: float, sigma_log: float) -> "scipy.stats.rv_frozen":
    return scipy.stats.lognorm(s=sigma_log, scale=_lognorm_scale_at(k0, sigma_log))


def build_kcat_prior(
    kcat0: np.ndarray, sigma0_log: np.ndarray,
) -> pyabc.Distribution:
    """Independent per-kcat lognormal prior, centred on ``kcat0``.

    Mirrors generation 1's inline sampling in
    ``bayesianSensitivityTuning.m``.

    Parameters
    ----------
    kcat0
        Current (pre-tuning) kcat values for every *tunable* row
        (strictly positive; filter out unassigned kcats first).
    sigma0_log
        Per-row prior std dev in log-space, from :func:`build_sigma0_log`
        (already filtered to the same rows as ``kcat0``).

    Returns
    -------
    pyabc.Distribution
        One lognormal RV per row, keyed ``"k0"``, ``"k1"``, ... in
        row order.
    """
    _check_shapes(kcat0, sigma0_log)
    rvs = {
        f"k{i}": pyabc.RV(
            "lognorm",
            s=float(sigma0_log[i]),
            scale=_lognorm_scale_at(kcat0[i], sigma0_log[i]),
        )
        for i in range(len(kcat0))
    }
    return pyabc.Distribution(**rvs)


def kcat_prior_logpdf(
    theta: np.ndarray, kcat0: np.ndarray, sigma0_log: np.ndarray,
) -> float:
    """Log prior density of a raw kcat vector.

    The same independent per-parameter lognormal prior
    :func:`build_kcat_prior` constructs, evaluated directly via scipy
    rather than through ``pyabc.Distribution``'s dict-based interface
    -- for use as the ``prior_logpdf`` callable in
    ``importance_weights.compute_importance_weights`` (Axis 2 variant
    B), which needs a plain array in, scalar out function.

    Parameters
    ----------
    theta
        One kcat vector, shape ``(n_params,)``.
    kcat0, sigma0_log
        As in :func:`build_kcat_prior`.

    Returns
    -------
    float
    """
    _check_shapes(kcat0, sigma0_log)
    if theta.shape != kcat0.shape:
        raise ValueError(
            f"theta shape {theta.shape} must match kcat0 shape {kcat0.shape}."
        )
    return float(np.sum([
        scipy.stats.lognorm.logpdf(
            theta[i], s=sigma0_log[i], scale=_lognorm_scale_at(kcat0[i], sigma0_log[i]),
        )
        for i in range(len(kcat0))
    ]))


class SpikeSlabRV(RVBase):
    """Spike-and-slab prior for one kcat (Axis 2, variant B).

    With probability ``spike_weight`` the parameter is drawn from a
    tight "spike" concentrated at ``kcat0`` (log-space std
    ``spike_sigma_log``, a small fraction of the slab's); otherwise
    from the ordinary "slab" -- the same lognormal prior
    :func:`build_kcat_prior` would use. A tight spike lets the
    posterior concentrate mass at "no real change" for parameters the
    data doesn't inform, the probabilistic complement to Axis 1
    variant B's importance weighting: instead of MATLAB's post-hoc
    "snap negligible changes to prior" step, an insignificant
    parameter's *posterior* already sits there.

    Not built via ``pyabc.RV`` (which only wraps ``scipy.stats``
    distributions looked up by name): a genuine two-component mixture
    has no single named ``scipy.stats`` distribution to delegate to,
    so this implements pyABC's ``RVBase`` interface directly --
    exactly the extension point pyABC itself provides for this case.
    """

    def __init__(
        self,
        kcat0: float,
        sigma_slab_log: float,
        *,
        spike_weight: float = 0.5,
        spike_sigma_log: float = 0.025,
    ):
        if not 0.0 <= spike_weight <= 1.0:
            raise ValueError(f"spike_weight must be in [0, 1]; got {spike_weight}.")
        if kcat0 <= 0:
            raise ValueError(f"kcat0 must be strictly positive; got {kcat0}.")
        self.kcat0 = float(kcat0)
        self.sigma_slab_log = float(sigma_slab_log)
        self.spike_weight = float(spike_weight)
        self.spike_sigma_log = float(spike_sigma_log)

    @property
    def _spike(self):
        return _lognorm_mean_at(self.kcat0, self.spike_sigma_log)

    @property
    def _slab(self):
        return _lognorm_mean_at(self.kcat0, self.sigma_slab_log)

    def copy(self) -> "SpikeSlabRV":
        return SpikeSlabRV(
            self.kcat0, self.sigma_slab_log,
            spike_weight=self.spike_weight, spike_sigma_log=self.spike_sigma_log,
        )

    def rvs(self, *args, **kwargs) -> float:  # noqa: ARG002
        component = self._spike if np.random.rand() < self.spike_weight else self._slab
        return float(component.rvs())

    def pdf(self, x: float, *args, **kwargs) -> float:  # noqa: ARG002
        return (
            self.spike_weight * self._spike.pdf(x)
            + (1 - self.spike_weight) * self._slab.pdf(x)
        )

    def pmf(self, x: float, *args, **kwargs) -> float:  # noqa: ARG002
        return self.pdf(x)

    def cdf(self, x: float, *args, **kwargs) -> float:  # noqa: ARG002
        return (
            self.spike_weight * self._spike.cdf(x)
            + (1 - self.spike_weight) * self._slab.cdf(x)
        )

    def __repr__(self) -> str:
        return (
            f"<SpikeSlabRV kcat0={self.kcat0}, sigma_slab_log={self.sigma_slab_log}, "
            f"spike_weight={self.spike_weight}, spike_sigma_log={self.spike_sigma_log}>"
        )


def build_kcat_sparsity_prior(
    kcat0: np.ndarray,
    sigma0_log: np.ndarray,
    *,
    spike_weight: float = 0.5,
    spike_sigma_frac: float = 0.05,
) -> pyabc.Distribution:
    """Independent per-kcat spike-and-slab prior (Axis 2, variant B).

    Parameters
    ----------
    kcat0, sigma0_log
        As in :func:`build_kcat_prior`.
    spike_weight
        Prior probability of drawing from the "no real change" spike
        rather than the ordinary slab, per parameter.
    spike_sigma_frac
        The spike's log-space std dev, as a fraction of the slab's
        (``sigma0_log[i]``) -- how tight "no real change" is.

    Returns
    -------
    pyabc.Distribution
        One :class:`SpikeSlabRV` per row, keyed ``"k0"``, ``"k1"``,
        ... in row order.
    """
    _check_shapes(kcat0, sigma0_log)
    rvs = {
        f"k{i}": SpikeSlabRV(
            float(kcat0[i]), float(sigma0_log[i]),
            spike_weight=spike_weight,
            spike_sigma_log=float(sigma0_log[i]) * spike_sigma_frac,
        )
        for i in range(len(kcat0))
    }
    return pyabc.Distribution(**rvs)
