"""Merge DLKcat and fuzzy-BRENDA kcat lists by priority.

Ported from GECKO MATLAB:
src/geckomat/gather_kcats/mergeDLKcatAndFuzzyKcats.m.
"""
from __future__ import annotations

import pandas as pd


_REQUIRED_COLUMNS = [
    "rxn_id", "source", "eccode", "substrates", "genes",
    "kcat", "wildcard_level", "origin",
]


def merge_dlkcat_and_fuzzy_kcats(
    dlkcat_kcats: pd.DataFrame,
    fuzzy_kcats: pd.DataFrame,
    *,
    top_origin_limit: int = 6,
    bottom_origin_limit: int = 6,
    wildcard_limit: int = 3,
) -> pd.DataFrame:
    """Merge DLKcat predictions with fuzzy BRENDA matches.

    Ported from GECKO MATLAB:
    src/geckomat/gather_kcats/mergeDLKcatAndFuzzyKcats.m.

    Priority (highest first):

    1. **Fuzzy (prio1):** ``wildcard_level == 0`` AND
       ``origin <= top_origin_limit``. The "best" BRENDA matches.
    2. **DLKcat (prio2):** any DLKcat row whose reaction is NOT
       already covered by prio1.
    3. **Fuzzy (prio3):** ``(wc == 0 AND top < origin <= bottom)``
       OR ``(0 < wc <= wildcard_limit AND origin <= bottom)``, AND
       the reaction is NOT covered by prio2. The wildcard-vs-origin
       prioritization is already done by ``fuzzy_kcat_matching``,
       so the two branches are joined here.

    Output rows are arranged as ``[prio1 fuzzy; prio3 fuzzy; prio2
    dlkcat]`` (fuzzy followed by DLKcat). Reaction order is not
    preserved (matching MATLAB).

    MATLAB-COMPAT: MATLAB's scalar ``source = "Merged DLKcat and
    fuzzy"`` field is dropped; the per-row ``source`` column already
    distinguishes ``"brenda"`` (fuzzy) from ``"DLKcat"``.

    MATLAB-COMPAT: MATLAB uses ``1000`` as a sentinel for NaN
    wildcard/origin so unmatched rows fail every priority
    comparison. geckopy uses ``fillna(1000)`` to the same effect.

    Parameters
    ----------
    dlkcat_kcats
        DataFrame produced by ``read_dlkcat_output``.
    fuzzy_kcats
        DataFrame produced by ``fuzzy_kcat_matching``.
    top_origin_limit
        Origin upper bound for prio1 (BRENDA preferred over DLKcat).
        Must be in ``[1, 6]``.
    bottom_origin_limit
        Origin upper bound for prio3 (BRENDA as fallback). Must be
        in ``[1, 6]``.
    wildcard_limit
        Maximum allowed wildcard count for prio3's wildcard branch.
        Must be in ``[0, 3]``.

    Returns
    -------
    pandas.DataFrame
        Concatenated kcat list in the canonical 8-column schema.

    Raises
    ------
    ValueError
        If any limit is outside its allowed range, or if either
        input is missing required columns.
    """
    if not (1 <= top_origin_limit <= 6):
        raise ValueError(
            f"top_origin_limit must be in [1, 6], got {top_origin_limit}"
        )
    if not (1 <= bottom_origin_limit <= 6):
        raise ValueError(
            f"bottom_origin_limit must be in [1, 6], got {bottom_origin_limit}"
        )
    if not (0 <= wildcard_limit <= 3):
        raise ValueError(
            f"wildcard_limit must be in [0, 3], got {wildcard_limit}"
        )

    for label, df in (("dlkcat_kcats", dlkcat_kcats),
                      ("fuzzy_kcats", fuzzy_kcats)):
        missing = [c for c in _REQUIRED_COLUMNS if c not in df.columns]
        if missing:
            raise ValueError(
                f"{label} missing required column(s): {missing}"
            )

    fuzzy_wc = fuzzy_kcats["wildcard_level"].fillna(1000)
    fuzzy_origin = fuzzy_kcats["origin"].fillna(1000)

    prio1_mask = (fuzzy_wc == 0) & (fuzzy_origin <= top_origin_limit)
    prio3_candidate = (
        ((fuzzy_wc == 0)
         & (fuzzy_origin > top_origin_limit)
         & (fuzzy_origin <= bottom_origin_limit))
        | ((fuzzy_wc > 0)
           & (fuzzy_wc <= wildcard_limit)
           & (fuzzy_origin <= bottom_origin_limit))
    )

    prio1_rxns = set(fuzzy_kcats.loc[prio1_mask, "rxn_id"])
    prio2_mask = ~dlkcat_kcats["rxn_id"].isin(prio1_rxns)

    prio2_rxns = set(dlkcat_kcats.loc[prio2_mask, "rxn_id"])
    prio3_mask = prio3_candidate & ~fuzzy_kcats["rxn_id"].isin(prio2_rxns)

    fuzzy_kept = fuzzy_kcats[prio1_mask | prio3_mask][_REQUIRED_COLUMNS]
    dlkcat_kept = dlkcat_kcats[prio2_mask][_REQUIRED_COLUMNS]

    # Filter out empty parts before concat to avoid the pandas
    # "concatenation with empty or all-NA entries" FutureWarning.
    if fuzzy_kept.empty and dlkcat_kept.empty:
        return fuzzy_kept.copy()
    if fuzzy_kept.empty:
        return dlkcat_kept.reset_index(drop=True)
    if dlkcat_kept.empty:
        return fuzzy_kept.reset_index(drop=True)
    return pd.concat([fuzzy_kept, dlkcat_kept], ignore_index=True)
