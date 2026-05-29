"""FluxData container for proteomics-paired flux measurements.

Used by ``constrain_flux_data``. A future ``load_flux_data`` loader
(parallel to ``load_pax_db``) would parse a TSV into this shape.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass
class FluxData:
    """Per-condition flux measurements with per-reaction exchange data.

    Attributes
    ----------
    conds
        Names of the experimental conditions / samples.
    p_tot
        Total protein content per condition, g/gDCW. Shape
        ``(n_conds,)``.
    gr_rate
        Measured growth rate per condition, 1/h. Shape ``(n_conds,)``.
    exch_fluxes
        Exchange flux per (condition, reaction), mmol/gDCW/h. Shape
        ``(n_conds, n_rxns)``. NaN means "no measurement"; ``+/-1000``
        is a sentinel for "unconstrained" (used by
        ``constrain_flux_data``).
    exch_mets
        Metabolite names matching the columns of ``exch_fluxes``.
    exch_rxn_ids
        Exchange reaction IDs in the model, matching the columns of
        ``exch_fluxes`` and ``exch_mets``.
    """

    conds: list[str] = field(default_factory=list)
    p_tot: np.ndarray = field(
        default_factory=lambda: np.empty(0, dtype=float)
    )
    gr_rate: np.ndarray = field(
        default_factory=lambda: np.empty(0, dtype=float)
    )
    exch_fluxes: np.ndarray = field(
        default_factory=lambda: np.empty((0, 0), dtype=float)
    )
    exch_mets: list[str] = field(default_factory=list)
    exch_rxn_ids: list[str] = field(default_factory=list)
