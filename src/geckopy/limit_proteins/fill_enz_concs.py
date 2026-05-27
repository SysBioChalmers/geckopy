"""Fill ``model.ec.concs`` from a proteomics measurement.

Takes a ``ProtData`` (a parsed proteomics table from
``load_prot_data`` / ``load_pax_db``) and copies the measured
concentrations into ``model.ec.concs``, matching by uniprot id.
Enzymes that aren't in the proteomics data keep their default
``NaN`` (no measurement available).

This only writes the numbers into the data array. To actually
make the LP respect them, call ``constrain_enz_concs`` next.

Ported from GECKO MATLAB:
src/geckomat/limit_proteins/fillEnzConcs.m.
"""
from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from ..databases.pax_db_loader import ProtData
    from ..ec_model.ec_model import EcModel

logger = logging.getLogger(__name__)


def fill_enz_concs(
    model: "EcModel",
    prot_data: "ProtData",
    *,
    data_col: int = 0,
) -> None:
    """Fill ``model.ec.concs`` from ``prot_data.abundances`` (mg/gDCW).

    Resets ``model.ec.concs`` to all-NaN of length ``n_enzymes``,
    then writes each (uniprot_id, abundance) pair into the slot
    matching ``model.ec.enzymes``. UniProt IDs absent from
    ``model.ec.enzymes`` are silently ignored. Enzymes absent from
    ``prot_data`` keep their NaN.

    For 2-D ``prot_data.abundances`` (multiple conditions /
    samples), ``data_col`` selects which column to use.

    Ported from GECKO MATLAB:
    src/geckomat/limit_proteins/fillEnzConcs.m.

    MATLAB-COMPAT: MATLAB's ``dataCol`` is 1-indexed and defaults to
    1. geckopy's ``data_col`` is 0-indexed (Python convention) and
    defaults to 0.

    Parameters
    ----------
    model
        EcModel; ``model.ec.concs`` is replaced in place.
    prot_data
        Pre-loaded proteomics data (e.g. from ``load_pax_db`` or a
        manual TSV reader). Abundances are expected in mg/gDCW.
    data_col
        0-indexed column of ``prot_data.abundances`` to read. Must be
        ``0`` for 1-D abundances.

    Raises
    ------
    IndexError
        If ``data_col`` is out of range for the supplied abundances.
    """
    n_enz = model.ec.n_enzymes
    model.ec.concs = np.full(n_enz, np.nan, dtype=float)

    if prot_data.abundances.ndim == 2:
        if not 0 <= data_col < prot_data.abundances.shape[1]:
            raise IndexError(
                f"data_col {data_col} out of range for abundances "
                f"with {prot_data.abundances.shape[1]} column(s)."
            )
        column = prot_data.abundances[:, data_col]
    else:
        if data_col != 0:
            raise IndexError(
                f"prot_data.abundances is 1-D; data_col must be 0, "
                f"got {data_col}."
            )
        column = prot_data.abundances

    if n_enz == 0 or len(prot_data.uniprot_ids) == 0:
        return

    enzyme_to_idx = {enz: i for i, enz in enumerate(model.ec.enzymes)}
    for uid, conc in zip(prot_data.uniprot_ids, column):
        idx = enzyme_to_idx.get(uid)
        if idx is None:
            continue
        value = float(conc)
        if value < 0:
            # A measured concentration can't be negative; skip it so it
            # doesn't become an infeasible bound in constrain_enz_concs.
            logger.warning(
                "fill_enz_concs: negative concentration %g for %r; skipped.",
                value, uid,
            )
            continue
        model.ec.concs[idx] = value
