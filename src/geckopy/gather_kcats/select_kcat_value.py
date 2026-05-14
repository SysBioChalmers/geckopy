"""Pick one kcat per reaction from a kcat list and write to model.ec.kcat.

Ported from GECKO MATLAB:
src/geckomat/gather_kcats/selectKcatValue.m.
"""
from __future__ import annotations

from typing import TYPE_CHECKING, Literal, Union

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    from ..ec_model.ec_model import EcModel


_CRITERIA = ("max", "min", "median", "mean")
_OverwriteMode = Union[bool, Literal["if_higher"]]


def select_kcat_value(
    model: "EcModel",
    kcat_list: pd.DataFrame,
    *,
    criteria: Literal["max", "min", "median", "mean"] = "max",
    overwrite: _OverwriteMode = True,
) -> list[str]:
    """Choose one kcat value per reaction and write it to ``model.ec.kcat``.

    Ported from GECKO MATLAB:
    src/geckomat/gather_kcats/selectKcatValue.m.

    For each reaction with one or more entries in ``kcat_list``,
    aggregate the candidate kcats per ``criteria`` (``"max"`` is the
    default and most common) and write the result into
    ``model.ec.kcat`` at the matching ``model.ec.rxns`` position.
    The contributing entry's ``source`` is recorded in
    ``model.ec.source``.

    Rows with ``kcat == 0`` are treated as no-match sentinels and
    dropped before aggregation, matching MATLAB.

    MATLAB-COMPAT: GECKO MATLAB takes the source via
    ``[v, j] = median(...)`` / ``mean(...)``, which doesn't actually
    return a meaningful index in standard MATLAB. geckopy attributes
    the source to the **first row** of each group when ``criteria``
    is ``"median"`` or ``"mean"``. For ``"max"`` and ``"min"`` the
    source is unambiguous (the row whose kcat won).

    MATLAB-COMPAT: GECKO MATLAB's ``kcatList`` has both a scalar
    ``source`` field and an optional per-row ``kcatSource`` cell;
    geckopy collapses both into a single per-row ``source`` column,
    matching the schema produced by ``fuzzy_kcat_matching`` and the
    upcoming DLKcat/manual loaders.

    MATLAB-COMPAT: GECKO MATLAB returns ``(model, rxnIdx)`` where
    ``rxnIdx`` is a list of integer positions. geckopy returns a
    list of reaction IDs (strings) which are more directly useful
    for downstream callers.

    Parameters
    ----------
    model
        EcModel with ``model.ec.rxns``, ``model.ec.kcat``, and
        ``model.ec.source`` already allocated. Mutated in place.
    kcat_list
        DataFrame with at least the columns ``rxn_id``, ``kcat``
        (1/s, ``0`` treated as no-match), and ``source`` (per-row).
        Multiple rows for the same reaction are aggregated per
        ``criteria``.
    criteria
        How to aggregate multiple candidates per reaction:
        ``"max"``, ``"min"``, ``"median"``, or ``"mean"``.
    overwrite
        ``True`` (default) overwrites unconditionally. ``False``
        only fills entries currently 0 or NaN. ``"if_higher"`` only
        updates when the new kcat exceeds the existing one (with
        0/NaN counted as "infinitely small" so they always update).

    Returns
    -------
    list of str
        Reaction IDs whose kcat was updated.

    Raises
    ------
    ValueError
        If ``criteria`` or ``overwrite`` is invalid, or if any
        ``rxn_id`` in ``kcat_list`` is absent from ``model.ec.rxns``.
    KeyError
        If ``kcat_list`` is missing the ``rxn_id``, ``kcat``, or
        ``source`` column.
    """
    if criteria not in _CRITERIA:
        raise ValueError(
            f"criteria must be one of {_CRITERIA}; got {criteria!r}"
        )
    if overwrite not in (True, False, "if_higher"):
        raise ValueError(
            f"overwrite must be True, False, or 'if_higher'; got {overwrite!r}"
        )

    for required in ("rxn_id", "kcat", "source"):
        if required not in kcat_list.columns:
            raise KeyError(
                f"kcat_list missing required column {required!r}"
            )

    nonzero = kcat_list[kcat_list["kcat"] != 0]
    if nonzero.empty:
        return []

    ec_rxn_ids = set(model.ec.rxns)
    unknown = sorted(set(nonzero["rxn_id"].astype(str)) - ec_rxn_ids)
    if unknown:
        preview = unknown[:5]
        raise ValueError(
            f"{len(unknown)} reaction ID(s) in kcat_list are not present "
            f"in model.ec.rxns (examples: {preview})"
        )

    rxn_to_idx = {rid: i for i, rid in enumerate(model.ec.rxns)}

    selected: dict[str, tuple[float, str]] = {}
    for rxn_id, group in nonzero.groupby("rxn_id", sort=False):
        kcats = group["kcat"].values.astype(float)
        sources = group["source"].astype(str).values
        if criteria == "max":
            j = int(np.argmax(kcats))
            kcat_val, source_val = float(kcats[j]), str(sources[j])
        elif criteria == "min":
            j = int(np.argmin(kcats))
            kcat_val, source_val = float(kcats[j]), str(sources[j])
        elif criteria == "median":
            kcat_val = float(np.median(kcats))
            source_val = str(sources[0])
        else:  # "mean"
            kcat_val = float(np.mean(kcats))
            source_val = str(sources[0])
        selected[str(rxn_id)] = (kcat_val, source_val)

    updated: list[str] = []
    for rxn_id, (kcat_val, source_val) in selected.items():
        i = rxn_to_idx[rxn_id]
        if _should_overwrite(model.ec.kcat[i], kcat_val, overwrite):
            model.ec.kcat[i] = kcat_val
            model.ec.source[i] = source_val
            updated.append(rxn_id)

    return updated


def _should_overwrite(
    current: float, new: float, mode: _OverwriteMode,
) -> bool:
    """Decide whether to write `new` over `current` per the overwrite mode."""
    is_unset = np.isnan(current) or current == 0
    if mode is True:
        return True
    if mode is False:
        return is_unset
    # "if_higher"
    return is_unset or new > current
