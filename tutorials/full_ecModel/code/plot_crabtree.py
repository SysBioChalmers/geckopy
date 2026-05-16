"""Crabtree-effect simulation and plot.

Python port of tutorials/full_ecModel/code/plotCrabtree.m from
GECKO MATLAB. Sweeps growth rate, finds the FBA solution that
maximises substrate uptake (and then minimises protein-pool
usage subject to a 1% slack on the solved uptake), and plots:

  Left:  exchange fluxes (glucose, O2, CO2, ethanol) vs growth
         rate, with experimental data points (van Hoek 1998).
  Right: fraction of protein pool used vs growth rate.

Functional equivalence with the MATLAB plot, not pixel parity.
"""
from __future__ import annotations

from pathlib import Path
from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from geckopy import EcModel


_GROWTH_RXN = "r_2111"           # biomass-equation-coupled growth rxn
_GLC_EXCHANGE = "r_1714"         # glucose uptake
_O2_EXCHANGE = "r_1992"          # oxygen uptake
_CO2_EXCHANGE = "r_1672"         # CO2 secretion
_ETHANOL_EXCHANGE = "r_1761"     # ethanol secretion
_POOL_EXCHANGE = "prot_pool_exchange"


def plot_crabtree(
    ec_model: EcModel,
    *,
    data_path: Optional[Path] = None,
    save_path: Optional[Path] = None,
    growth_rates: Optional[np.ndarray] = None,
) -> tuple[pd.DataFrame, np.ndarray]:
    """Run the Crabtree simulation and plot the result.

    Parameters
    ----------
    ec_model
        An ecModel with a populated protein pool.
    data_path
        Path to ``vanHoek1998.tsv`` (the experimental reference);
        if None, no experimental scatter is added.
    save_path
        If supplied, the figure is saved to this path (PDF/PNG).
    growth_rates
        Vector of growth rates to sweep. Default: ``np.arange(0,
        0.401, 0.025)`` (matching MATLAB).

    Returns
    -------
    (fluxes_df, growth_rates)
        ``fluxes_df``: DataFrame with one row per reaction and one
        column per growth rate (NaN where infeasible). ``growth_rates``:
        the swept growth rates.
    """
    if growth_rates is None:
        growth_rates = np.arange(0.0, 0.401, 0.025)

    pool_lb_abs = abs(ec_model.reactions.get_by_id(_POOL_EXCHANGE).lower_bound)
    pool_ub_abs = abs(ec_model.reactions.get_by_id(_POOL_EXCHANGE).upper_bound)
    total_protein = pool_ub_abs if pool_ub_abs > 0 else pool_lb_abs

    rxn_ids = [r.id for r in ec_model.reactions]
    out = np.full((len(rxn_ids), len(growth_rates)), np.nan, dtype=float)

    for i, gr in enumerate(growth_rates):
        with ec_model:
            ec_model.reactions.get_by_id(_GROWTH_RXN).lower_bound = float(gr)
            # Step 1: maximise glucose uptake (= least-negative flux).
            ec_model.objective = _GLC_EXCHANGE
            sol = ec_model.optimize()
            if sol.status != "optimal":
                continue
            glc = sol.fluxes[_GLC_EXCHANGE]
            # Step 2: lock glucose at +1% slack, minimise pool usage.
            ec_model.reactions.get_by_id(_GLC_EXCHANGE).lower_bound = (
                glc * 1.01
            )
            ec_model.objective = {
                ec_model.reactions.get_by_id(_POOL_EXCHANGE): -1.0,
            }
            sol = ec_model.optimize()
            if sol.status != "optimal":
                continue
            out[:, i] = sol.fluxes.values

    fluxes = pd.DataFrame(
        out, index=rxn_ids, columns=[f"{gr:.3f}" for gr in growth_rates],
    )

    fig, (ax_left, ax_right) = plt.subplots(1, 2, figsize=(10, 4))

    plot_ids = [_O2_EXCHANGE, _CO2_EXCHANGE, _GLC_EXCHANGE, _ETHANOL_EXCHANGE]
    plot_names = ["O2 (r_1992)", "CO2 (r_1672)",
                  "glucose (r_1714)", "ethanol (r_1761)"]
    for rid, name, color in zip(
        plot_ids, plot_names,
        ["#1f77b4", "#d62728", "#e6a300", "#9467bd"],
    ):
        ax_left.plot(
            growth_rates, np.abs(fluxes.loc[rid].values),
            label=name, color=color, linewidth=2,
        )

    if data_path is not None and Path(data_path).is_file():
        exp = pd.read_csv(data_path, sep=";", skiprows=2, header=None)
        exp.columns = [
            "growth", "qO2", "qCO2", "qglucose", "qethanol",
            "qacetate", "qpyruvate", "qglycerol",
        ]
        for col, color in zip(
            ["qO2", "qCO2", "qglucose", "qethanol"],
            ["#1f77b4", "#d62728", "#e6a300", "#9467bd"],
        ):
            ax_left.scatter(exp["growth"], exp[col], color=color, s=30)

    ax_left.set_xlabel("Growth rate (/hour)")
    ax_left.set_ylabel("Absolute flux (mmol/gDCWh)")
    ax_left.set_ylim(0, 20)
    ax_left.legend(loc="upper left", fontsize=9)

    pool_used = np.abs(fluxes.loc[_POOL_EXCHANGE].values) / total_protein
    ax_right.plot(growth_rates, pool_used, linewidth=2)
    ax_right.set_xlabel("Growth rate (/hour)")
    ax_right.set_ylabel("Fraction of protein pool used")

    fig.tight_layout()
    if save_path is not None:
        fig.savefig(save_path)
    return fluxes, growth_rates
