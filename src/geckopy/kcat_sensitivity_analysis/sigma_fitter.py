"""Fit the average enzyme saturation factor (``sigma``).

The protein pool's upper bound is set to ``P_tot * f * sigma``,
where ``sigma`` is the average enzyme saturation factor — a
fudge factor between 0 and 1 capturing the fact that enzymes
don't usually run at their full Vmax in vivo. The default value
is 0.5; this function fits a better one.

The algorithm: try a range of sigma values, set the protein-pool
bound at each, solve the model, compare the predicted growth
rate to the experimental one, pick the sigma that minimises the
difference.

Ported from GECKO MATLAB:
src/geckomat/kcat_sensitivity_analysis/sigmaFitter.m.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Optional

import numpy as np

from ..ec_model.pipeline.protein_pool import set_prot_pool_size

if TYPE_CHECKING:
    from ..ec_model.ec_model import EcModel


@dataclass
class SigmaFitterResult:
    """Outcome of a sigma scan.

    Attributes
    ----------
    sigma
        The sigma value that minimised the absolute relative error
        between predicted and experimental growth.
    sigma_grid
        Sigma values that were tried.
    growth_grid
        LP-optimal growth at each sigma.
    error_grid
        Absolute relative error at each sigma, in percent.
    """

    sigma: float = 0.0
    sigma_grid: np.ndarray = field(
        default_factory=lambda: np.empty(0, dtype=float)
    )
    growth_grid: np.ndarray = field(
        default_factory=lambda: np.empty(0, dtype=float)
    )
    error_grid: np.ndarray = field(
        default_factory=lambda: np.empty(0, dtype=float)
    )


def fit_sigma(
    model: "EcModel",
    *,
    growth_rate: Optional[float] = None,
    p_tot: Optional[float] = None,
    f: Optional[float] = None,
    n_sigma_steps: int = 100,
) -> SigmaFitterResult:
    """Sweep sigma and pick the one matching ``growth_rate`` best.

    For each sigma in ``[1/n, 2/n, ..., 1.0]``, the protein pool size
    is set to ``p_tot * f * sigma`` and the LP is solved. The sigma
    minimising ``|relative_error|`` between predicted and target
    growth is recorded and re-applied to ``model`` before return.

    Ported from GECKO MATLAB:
    src/geckomat/kcat_sensitivity_analysis/sigmaFitter.m.

    MATLAB-COMPAT: GECKO MATLAB leaves the model at the LAST trial
    (sigma = 1.0) even though its docstring claims to return the
    model adapted to the optimal sigma. geckopy re-applies the best
    sigma at the end. Tracked in ``docs/future_improvements.md``.

    MATLAB-COMPAT: GECKO MATLAB takes a ``modelAdapter`` arg and a
    ``makePlot`` flag. geckopy reads the adapter from
    ``model.adapter`` and returns the diagnostic grids in a
    dataclass; callers plot via matplotlib if wanted.

    Parameters
    ----------
    model
        EcModel with the protein pool machinery installed and the
        objective set (typically the biomass reaction). Mutated in
        place: the optimal sigma is applied via
        ``set_prot_pool_size`` before return.
    growth_rate
        Experimental growth rate to match. Defaults to
        ``model.adapter.params.gr_exp``.
    p_tot
        Total protein content (g/gDCW). Defaults to
        ``model.adapter.params.p_tot``.
    f
        Mass fraction of model enzymes. Defaults to
        ``model.adapter.params.f``.
    n_sigma_steps
        Number of sigma values to try in ``(0, 1]`` (``i / n`` for
        ``i = 1..n``). Default 100.

    Returns
    -------
    SigmaFitterResult

    Raises
    ------
    ValueError
        If ``model.adapter`` is None and any default is needed; if
        ``n_sigma_steps`` is < 1.
    """
    if n_sigma_steps < 1:
        raise ValueError(f"n_sigma_steps must be >= 1, got {n_sigma_steps}")

    if (
        (growth_rate is None or p_tot is None or f is None)
        and model.adapter is None
    ):
        raise ValueError(
            "model.adapter is None and one of growth_rate / p_tot / f "
            "is not provided."
        )
    params = model.adapter.params if model.adapter is not None else None

    if growth_rate is None:
        growth_rate = float(params.gr_exp)
    if p_tot is None:
        p_tot = float(params.p_tot)
    if f is None:
        f = float(params.f)

    sigma_grid = np.array(
        [(i + 1) / n_sigma_steps for i in range(n_sigma_steps)], dtype=float,
    )
    growth_grid = np.zeros(n_sigma_steps, dtype=float)
    error_grid = np.zeros(n_sigma_steps, dtype=float)

    for k, sigma in enumerate(sigma_grid):
        with model:
            set_prot_pool_size(
                model, p_tot=p_tot, f=f, sigma=float(sigma),
            )
            sol = model.optimize()
            growth = float(sol.objective_value or 0.0)
        growth_grid[k] = growth
        if growth_rate != 0:
            error_grid[k] = abs(
                (growth_rate - growth) / growth_rate
            ) * 100.0
        else:
            error_grid[k] = abs(growth) * 100.0

    best_k = int(np.argmin(error_grid))
    best_sigma = float(sigma_grid[best_k])

    # Apply the best sigma to the model permanently (geckopy divergence).
    set_prot_pool_size(model, p_tot=p_tot, f=f, sigma=best_sigma)

    return SigmaFitterResult(
        sigma=best_sigma,
        sigma_grid=sigma_grid,
        growth_grid=growth_grid,
        error_grid=error_grid,
    )


def sigma_fitter(
    model: "EcModel",
    *,
    growth_rate: Optional[float] = None,
    p_tot: Optional[float] = None,
    f: Optional[float] = None,
    n_sigma_steps: int = 100,
) -> SigmaFitterResult:
    """Deprecated alias for :func:`fit_sigma`.

    Kept for backward compatibility. Will be removed in a future
    release; switch to ``fit_sigma``.
    """
    import warnings

    warnings.warn(
        "sigma_fitter is deprecated; use fit_sigma instead. "
        "The old name will be removed in a future release.",
        DeprecationWarning,
        stacklevel=2,
    )
    return fit_sigma(
        model,
        growth_rate=growth_rate, p_tot=p_tot, f=f,
        n_sigma_steps=n_sigma_steps,
    )
