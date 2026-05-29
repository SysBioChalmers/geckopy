"""ecFVA cumulative-distribution plot.

Python port of tutorials/full_ecModel/code/plotEcFVA.m (an inline
helper inside the larger MATLAB script). Plots empirical CDFs of
the flux variability range (max - min) for one or more models on
a log-x axis. Reactions whose |min| and |max| are both below
``zero_threshold`` are excluded (effectively unused at the
solved point).

Functional equivalence with the MATLAB plot.
"""
from __future__ import annotations

from pathlib import Path
from typing import Optional, Sequence

import matplotlib.pyplot as plt
import numpy as np


def plot_ec_fva(
    min_flux: np.ndarray,
    max_flux: np.ndarray,
    *,
    labels: Optional[Sequence[str]] = None,
    save_path: Optional[Path] = None,
    zero_threshold: float = 1e-10,
) -> None:
    """Plot CDFs of flux variability ranges.

    Parameters
    ----------
    min_flux, max_flux
        2-D arrays shaped ``(n_rxns, n_models)`` of FVA results,
        one column per model. 1-D arrays are accepted and treated
        as a single model.
    labels
        Per-model label strings. Defaults to ``Model #1``, ...
    save_path
        If supplied, the figure is saved to this path.
    zero_threshold
        Reactions where ``|min| < zero_threshold`` AND
        ``|max| < zero_threshold`` are treated as unused and
        excluded from the CDF.
    """
    min_flux = np.atleast_2d(min_flux)
    max_flux = np.atleast_2d(max_flux)
    if min_flux.shape[0] == 1 and min_flux.shape[1] > 1:
        min_flux = min_flux.T
    if max_flux.shape[0] == 1 and max_flux.shape[1] > 1:
        max_flux = max_flux.T

    n_models = min_flux.shape[1]
    if labels is None:
        labels = [f"Model #{i + 1}" for i in range(n_models)]

    fig, ax = plt.subplots(figsize=(7, 5))
    for i in range(n_models):
        zero = (
            (np.abs(min_flux[:, i]) < zero_threshold)
            & (np.abs(max_flux[:, i]) < zero_threshold)
        )
        flux_range = max_flux[:, i] - min_flux[:, i]
        flux_range = flux_range[~zero]
        flux_range = flux_range[~np.isnan(flux_range)]
        if flux_range.size == 0:
            continue
        sorted_range = np.sort(flux_range)
        y = np.arange(1, len(sorted_range) + 1) / len(sorted_range)
        median = float(np.median(sorted_range))
        ax.plot(
            sorted_range, y, linewidth=2,
            label=f"{labels[i]} (median: {median:.2g})",
        )

    ax.set_xscale("log")
    ax.set_xlim(1e-7, 1e4)
    ax.set_xlabel("Variability range (mmol/gDCWh)")
    ax.set_ylabel("Cumulative distribution")
    ax.set_title("Flux variability (cumulative distribution)")
    ax.legend(loc="upper left")
    fig.tight_layout()
    if save_path is not None:
        fig.savefig(save_path)
