"""GeckoTransition: diagonal log-space perturbation kernel.

Wires GECKO's proposal-kernel design into pyABC's ``Transition``
contract (``fit``/``rvs_single``/``pdf``). Ported design intent from
GECKO MATLAB's ``buildLowRankLogProposal`` (the diagonal, non-PCA path
that's actually reached, per the resolved decision -- see
docs/internal/bayesian_tuning_plan.md's "Not building parallel
variants for: Proposal kernel"): the fitted per-parameter bandwidth
blends the accepted particles' own observed spread with the prior's
sigma0_log, floored at a fraction of sigma0_log, exactly matching
``adaptFracEarly``/``sigmaFloorFrac`` -- MATLAB's own comments mark
these "FIXED ALGORITHM PARAMETERS (rarely changed)", not
project-configurable, hence constructor defaults here rather than new
``BayesianParams`` fields.

MATLAB's abandoned low-rank PCA kernel and its explicit
exploit/explore mixture are deliberately not ported: a diagonal kernel
is also the right shape at genome scale (``ec.kcat`` has thousands of
entries vs. a few hundred particles per generation, so a full/low-rank
covariance would be rank-deficient).

``fit()`` is variant-agnostic: it just consumes whatever ``(X, w)`` it
is given. Axis 2 variant A (``posterior.py``) calls it with uniform
weights over its blended point estimate; variant B
(``importance_weights.py``) calls it with the particles' actual
importance weights.
"""
from __future__ import annotations

from typing import Union

import numpy as np
import pandas as pd
import scipy.stats
from pyabc.parameters import Parameter
from pyabc.transition import Transition


class GeckoTransition(Transition):
    """Diagonal per-parameter lognormal perturbation kernel.

    Working in log-space (kcats are strictly positive and
    lognormally distributed) means "diagonal Normal kernel in
    log-space" is exactly a diagonal *lognormal* kernel in kcat-space;
    :func:`scipy.stats.lognorm` already encodes the Jacobian for that
    change of variables, so no separate correction term is needed.
    """

    def __init__(
        self,
        sigma0_log: np.ndarray,
        *,
        adapt_frac_early: float = 0.5,
        sigma_floor_frac: float = 0.15,
        bandwidth_scale: float = 1.0,
    ):
        """
        Parameters
        ----------
        sigma0_log
            Per-parameter prior std dev in log-space (from
            ``priors.build_sigma0_log``), shape ``(n_params,)``, in
            the same order as ``fit()``'s ``X`` columns will be.
        adapt_frac_early
            Blend weight between the accepted particles' own observed
            log-space std and ``sigma0_log`` (MATLAB's
            ``adaptFracEarly``, a fixed constant there).
        sigma_floor_frac
            Floor on the fitted bandwidth, as a fraction of
            ``sigma0_log`` (MATLAB's ``sigmaFloorFrac``).
        bandwidth_scale
            Extra multiplier on the fitted bandwidth -- widen (>1) or
            narrow (<1) relative to the blended value, e.g. for an
            explore/exploit-style adjustment layered on top.
        """
        self.sigma0_log = np.asarray(sigma0_log, dtype=float)
        self.adapt_frac_early = adapt_frac_early
        self.sigma_floor_frac = sigma_floor_frac
        self.bandwidth_scale = bandwidth_scale
        self._columns: list[str] | None = None
        self._log_bandwidth: np.ndarray | None = None

    def fit(self, X: pd.DataFrame, w: np.ndarray) -> None:
        if X.shape[1] != len(self.sigma0_log):
            raise ValueError(
                f"X has {X.shape[1]} columns; expected {len(self.sigma0_log)} "
                f"to match sigma0_log."
            )
        self.X = X
        self.w = np.asarray(w, dtype=float)
        self._columns = list(X.columns)

        log_x = np.log(X.to_numpy(dtype=float))
        w_norm = self.w / self.w.sum()
        mean = np.average(log_x, axis=0, weights=w_norm)
        var = np.average((log_x - mean) ** 2, axis=0, weights=w_norm)
        std_obs = np.sqrt(var)

        blended = (
            self.adapt_frac_early * std_obs
            + (1 - self.adapt_frac_early) * self.sigma0_log
        )
        floored = np.maximum(blended, self.sigma_floor_frac * self.sigma0_log)
        self._log_bandwidth = floored * self.bandwidth_scale

    def rvs_single(self) -> Parameter:
        if self._log_bandwidth is None:
            raise RuntimeError("GeckoTransition.rvs_single() called before fit().")
        w_norm = self.w / self.w.sum()
        parent_idx = np.random.choice(len(self.X), p=w_norm)
        parent = self.X.iloc[parent_idx].to_numpy(dtype=float)
        sample = np.array([
            scipy.stats.lognorm(s=h, scale=p).rvs()
            for h, p in zip(self._log_bandwidth, parent)
        ])
        return Parameter(dict(zip(self._columns, sample)))

    def component_logpdf(self, x: np.ndarray, parent: np.ndarray) -> float:
        """Log-density of the single component that perturbs
        ``parent`` to reach ``x`` -- one mixture term of :meth:`pdf`,
        exposed on its own for Axis 2 variant B's importance-weight
        formula (``importance_weights.compute_importance_weights``'
        ``transition_logpdf`` callable takes one parent at a time,
        summing the weighted mixture itself)."""
        if self._log_bandwidth is None:
            raise RuntimeError("GeckoTransition.component_logpdf() called before fit().")
        return float(np.sum([
            scipy.stats.lognorm.logpdf(xi, s=h, scale=pi)
            for xi, h, pi in zip(x, self._log_bandwidth, parent)
        ]))

    def pdf(
        self, x: Union[Parameter, "pd.Series", pd.DataFrame],
    ) -> Union[float, np.ndarray]:
        if self._log_bandwidth is None:
            raise RuntimeError("GeckoTransition.pdf() called before fit().")
        if isinstance(x, pd.DataFrame):
            return np.array([
                self._pdf_single(row.to_numpy(dtype=float))
                for _, row in x.iterrows()
            ])
        x_arr = np.array([x[c] for c in self._columns], dtype=float)
        return self._pdf_single(x_arr)

    def _pdf_single(self, x_arr: np.ndarray) -> float:
        w_norm = self.w / self.w.sum()
        log_components = np.array([
            self.component_logpdf(x_arr, self.X.iloc[j].to_numpy(dtype=float))
            for j in range(len(self.X))
        ])
        m = np.max(log_components)
        mixture = np.sum(w_norm * np.exp(log_components - m))
        return float(mixture * np.exp(m))
