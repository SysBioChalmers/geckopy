"""FluxData container for proteomics-paired flux measurements,
plus a ``load_flux_data`` parser for the GECKO `fluxData.tsv` format.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import numpy as np


# `<met_name> (<rxn_id>)` -> capture met_name and rxn_id.
_EXCH_HEADER_RE = re.compile(r"^(.*)\s*\(([^()]+)\)\s*$")


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
    bayesian_rmse_weight
        Optional per-condition weights for the Bayesian-kcat-tuning
        RMSE term. ``None`` if the source file did not have the
        ``bayesianRMSEweight`` column.
    source
        Optional per-condition free-text describing where the data
        came from. ``None`` if the source file did not have a
        ``source`` column.
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
    bayesian_rmse_weight: Optional[np.ndarray] = None
    source: Optional[list[str]] = None


def load_flux_data(path: str | Path) -> FluxData:
    """Parse a GECKO ``fluxData.tsv`` into a ``FluxData``.

    File format (tab-delimited, one header row):

        condition  Ptot  grRate  <met1> (<rxn1>)  <met2> (<rxn2>)  ...
        c1         0.5   0.4     -10.0            5.0              ...
        c2         0.6   0.5     -8.5             4.2              ...

    Optional extra columns (in any position):

    * ``bayesianRMSEweight`` (float): weights for Bayesian kcat tuning.
    * ``source`` (str): free-text describing the data origin.

    Both columns are dropped from the matrix and surfaced on the
    returned ``FluxData`` as ``bayesian_rmse_weight`` and ``source``.

    Ported from GECKO MATLAB:
    src/geckomat/utilities/loadFluxData.m.

    MATLAB-COMPAT: GECKO MATLAB takes a ``modelAdapter`` and defaults
    the path to ``adapter.params.path/data/fluxData.tsv``. geckopy
    requires the path explicitly; the caller resolves it.

    MATLAB-COMPAT: GECKO MATLAB stores the optional columns under
    camelCase keys (``bayesianRMSEweight``, ``source``); geckopy
    uses snake_case (``bayesian_rmse_weight``, ``source``).

    Parameters
    ----------
    path
        Path to the TSV file.

    Returns
    -------
    FluxData

    Raises
    ------
    FileNotFoundError
        If ``path`` does not exist.
    ValueError
        If the file has no header row or no exchange columns.
    """
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"fluxData file not found: {path}")

    rows: list[list[str]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        rows.append(line.rstrip("\n").split("\t"))

    if not rows:
        raise ValueError(f"{path} is empty.")
    if len(rows[0]) < 3:
        raise ValueError(
            f"{path}: header has only {len(rows[0])} column(s); "
            f"expected at least 'condition Ptot grRate'."
        )

    header = rows[0]
    data = rows[1:]

    # Pull out optional columns first; collect indices to drop.
    optional_indices: list[int] = []
    bayesian_idx = _index_of(header, "bayesianRMSEweight")
    source_idx = _index_of(header, "source")

    bayesian_weights: Optional[np.ndarray] = None
    source_values: Optional[list[str]] = None

    if bayesian_idx is not None:
        bayesian_weights = np.array(
            [_parse_float(r[bayesian_idx]) if bayesian_idx < len(r) else float("nan")
             for r in data],
            dtype=float,
        )
        optional_indices.append(bayesian_idx)
    if source_idx is not None:
        source_values = [
            r[source_idx] if source_idx < len(r) else "" for r in data
        ]
        optional_indices.append(source_idx)

    # Build the kept-column index list (preserving order).
    kept = [j for j in range(len(header)) if j not in optional_indices]
    kept_header = [header[j] for j in kept]

    # First three columns are condition / Ptot / grRate. Remaining are exchanges.
    if len(kept_header) < 4:
        raise ValueError(
            f"{path}: no exchange-flux columns found after the first three."
        )

    exch_headers = kept_header[3:]
    exch_mets, exch_rxn_ids = _split_exch_headers(exch_headers)

    conds = [r[kept[0]] if kept[0] < len(r) else "" for r in data]
    p_tot = np.array(
        [_parse_float(r[kept[1]]) if kept[1] < len(r) else float("nan")
         for r in data],
        dtype=float,
    )
    gr_rate = np.array(
        [_parse_float(r[kept[2]]) if kept[2] < len(r) else float("nan")
         for r in data],
        dtype=float,
    )
    exch_fluxes = np.full((len(data), len(exch_headers)), np.nan, dtype=float)
    for i, r in enumerate(data):
        for ci, j in enumerate(kept[3:]):
            if j < len(r):
                exch_fluxes[i, ci] = _parse_float(r[j])

    return FluxData(
        conds=conds,
        p_tot=p_tot,
        gr_rate=gr_rate,
        exch_fluxes=exch_fluxes,
        exch_mets=exch_mets,
        exch_rxn_ids=exch_rxn_ids,
        bayesian_rmse_weight=bayesian_weights,
        source=source_values,
    )


def _index_of(header: list[str], name: str) -> Optional[int]:
    """Return the index of `name` in `header`, or None if absent."""
    for j, h in enumerate(header):
        if h.strip() == name:
            return j
    return None


def _parse_float(s: str) -> float:
    """Parse a numeric cell; empty / non-numeric -> NaN."""
    s = s.strip()
    if not s:
        return float("nan")
    try:
        return float(s)
    except ValueError:
        return float("nan")


def _split_exch_headers(headers: list[str]) -> tuple[list[str], list[str]]:
    """Split each header `'<met> (<rxn>)'` into (met, rxn) pairs.

    Headers without parentheses pass through with the rxn_id equal
    to the full header (matching MATLAB's regex behaviour: the
    capture-group falls through unchanged).
    """
    mets: list[str] = []
    rxns: list[str] = []
    for h in headers:
        m = _EXCH_HEADER_RE.match(h.strip())
        if m:
            mets.append(m.group(1).strip())
            rxns.append(m.group(2).strip())
        else:
            mets.append(h.strip())
            rxns.append(h.strip())
    return mets, rxns
