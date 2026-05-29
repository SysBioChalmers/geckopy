"""Tests for merge_dlkcat_and_fuzzy_kcats."""
import pandas as pd
import pytest

from geckopy.gather_kcats import merge_dlkcat_and_fuzzy_kcats


# --------------------------------------------------------------------------- #
# Fixture builders
# --------------------------------------------------------------------------- #

_CANONICAL_COLUMNS = [
    "rxn_id", "source", "eccode", "substrates", "genes",
    "kcat", "wildcard_level", "origin",
]


def _fuzzy_row(
    rxn_id: str,
    kcat: float,
    *,
    wildcard_level: object,
    origin: object,
    eccode: str = "1.1.1.1",
    substrates: list[str] | None = None,
) -> dict:
    return {
        "rxn_id": rxn_id,
        "source": "brenda",
        "eccode": eccode,
        "substrates": substrates if substrates is not None else ["alpha"],
        "genes": [],
        "kcat": kcat,
        "wildcard_level": wildcard_level,
        "origin": origin,
    }


def _dlkcat_row(rxn_id: str, kcat: float, gene: str = "g1") -> dict:
    return {
        "rxn_id": rxn_id,
        "source": "DLKcat",
        "eccode": "",
        "substrates": ["alpha"],
        "genes": [gene],
        "kcat": kcat,
        "wildcard_level": pd.NA,
        "origin": pd.NA,
    }


def _fuzzy_df(rows: list[dict]) -> pd.DataFrame:
    df = pd.DataFrame(rows, columns=_CANONICAL_COLUMNS) if rows else pd.DataFrame(
        columns=_CANONICAL_COLUMNS,
    )
    if "wildcard_level" in df:
        df["wildcard_level"] = df["wildcard_level"].astype("Int64")
    if "origin" in df:
        df["origin"] = df["origin"].astype("Int64")
    return df


def _dlkcat_df(rows: list[dict]) -> pd.DataFrame:
    df = pd.DataFrame(rows, columns=_CANONICAL_COLUMNS) if rows else pd.DataFrame(
        columns=_CANONICAL_COLUMNS,
    )
    if "wildcard_level" in df:
        df["wildcard_level"] = df["wildcard_level"].astype("Int64")
    if "origin" in df:
        df["origin"] = df["origin"].astype("Int64")
    return df


# --------------------------------------------------------------------------- #
# Validation
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("top", [0, 7, -1])
def test_top_origin_out_of_range_raises(top):
    with pytest.raises(ValueError, match="top_origin_limit"):
        merge_dlkcat_and_fuzzy_kcats(
            _dlkcat_df([]), _fuzzy_df([]), top_origin_limit=top,
        )


@pytest.mark.parametrize("bot", [0, 7])
def test_bottom_origin_out_of_range_raises(bot):
    with pytest.raises(ValueError, match="bottom_origin_limit"):
        merge_dlkcat_and_fuzzy_kcats(
            _dlkcat_df([]), _fuzzy_df([]), bottom_origin_limit=bot,
        )


@pytest.mark.parametrize("wc", [-1, 4])
def test_wildcard_limit_out_of_range_raises(wc):
    with pytest.raises(ValueError, match="wildcard_limit"):
        merge_dlkcat_and_fuzzy_kcats(
            _dlkcat_df([]), _fuzzy_df([]), wildcard_limit=wc,
        )


def test_missing_required_column_raises():
    bad_df = pd.DataFrame({"rxn_id": ["r1"]})  # missing most columns
    with pytest.raises(ValueError, match="missing required"):
        merge_dlkcat_and_fuzzy_kcats(bad_df, _fuzzy_df([]))


# --------------------------------------------------------------------------- #
# Trivial cases
# --------------------------------------------------------------------------- #

def test_both_empty_returns_empty():
    out = merge_dlkcat_and_fuzzy_kcats(_dlkcat_df([]), _fuzzy_df([]))
    assert out.empty
    assert list(out.columns) == _CANONICAL_COLUMNS


def test_only_dlkcat_input():
    out = merge_dlkcat_and_fuzzy_kcats(
        _dlkcat_df([_dlkcat_row("r1", 5.0), _dlkcat_row("r2", 7.0)]),
        _fuzzy_df([]),
    )
    assert sorted(out["rxn_id"]) == ["r1", "r2"]
    assert all(out["source"] == "DLKcat")


def test_only_fuzzy_input():
    out = merge_dlkcat_and_fuzzy_kcats(
        _dlkcat_df([]),
        _fuzzy_df([_fuzzy_row("r1", 5.0, wildcard_level=0, origin=1)]),
    )
    assert list(out["rxn_id"]) == ["r1"]
    assert list(out["source"]) == ["brenda"]


# --------------------------------------------------------------------------- #
# Prio 1: fuzzy wc=0 + origin<=top wins over dlkcat
# --------------------------------------------------------------------------- #

def test_fuzzy_prio1_wins_when_reaction_also_in_dlkcat():
    out = merge_dlkcat_and_fuzzy_kcats(
        _dlkcat_df([_dlkcat_row("r1", 100.0)]),
        _fuzzy_df([_fuzzy_row("r1", 5.0, wildcard_level=0, origin=1)]),
        top_origin_limit=3,
    )
    # Fuzzy prio1 wins: only one row for r1, from fuzzy.
    assert len(out) == 1
    row = out.iloc[0]
    assert row["source"] == "brenda"
    assert row["kcat"] == 5.0


def test_fuzzy_prio1_uses_default_top_origin_limit_6():
    """With default top_origin_limit=6, origin=6 still qualifies for prio1."""
    out = merge_dlkcat_and_fuzzy_kcats(
        _dlkcat_df([_dlkcat_row("r1", 100.0)]),
        _fuzzy_df([_fuzzy_row("r1", 5.0, wildcard_level=0, origin=6)]),
    )
    assert len(out) == 1
    assert out.iloc[0]["source"] == "brenda"


def test_fuzzy_with_wildcard_does_not_qualify_for_prio1():
    """Wildcard fuzzy is NOT prio1; DLKcat wins."""
    out = merge_dlkcat_and_fuzzy_kcats(
        _dlkcat_df([_dlkcat_row("r1", 100.0)]),
        _fuzzy_df([_fuzzy_row("r1", 5.0, wildcard_level=1, origin=1)]),
        top_origin_limit=3,
        wildcard_limit=0,  # prevent prio3 by setting wildcard_limit too low
    )
    # No prio1 (wc != 0), no prio3 (wc > wildcard_limit) -> only DLKcat.
    assert len(out) == 1
    assert out.iloc[0]["source"] == "DLKcat"


# --------------------------------------------------------------------------- #
# Prio 2: DLKcat fills in reactions not in fuzzy prio1
# --------------------------------------------------------------------------- #

def test_dlkcat_prio2_when_fuzzy_origin_too_high():
    """Fuzzy origin > top_origin_limit -> not prio1.
    With bottom_origin_limit < fuzzy origin -> not prio3 either.
    So DLKcat wins (prio2)."""
    out = merge_dlkcat_and_fuzzy_kcats(
        _dlkcat_df([_dlkcat_row("r1", 100.0)]),
        _fuzzy_df([_fuzzy_row("r1", 5.0, wildcard_level=0, origin=6)]),
        top_origin_limit=3,
        bottom_origin_limit=3,
    )
    assert len(out) == 1
    assert out.iloc[0]["source"] == "DLKcat"


def test_dlkcat_prio2_when_fuzzy_missing():
    out = merge_dlkcat_and_fuzzy_kcats(
        _dlkcat_df([_dlkcat_row("r1", 100.0)]),
        _fuzzy_df([_fuzzy_row("r2", 5.0, wildcard_level=0, origin=1)]),
    )
    # r1 only in DLKcat, r2 only in fuzzy prio1 -> both kept.
    assert sorted(out["rxn_id"]) == ["r1", "r2"]


# --------------------------------------------------------------------------- #
# Prio 3: fuzzy fallback when DLKcat doesn't cover it
# --------------------------------------------------------------------------- #

def test_fuzzy_prio3_via_origin_branch():
    """top < origin <= bottom, wc = 0 -> prio3.
    DLKcat does NOT cover r1, so fuzzy prio3 wins."""
    out = merge_dlkcat_and_fuzzy_kcats(
        _dlkcat_df([]),
        _fuzzy_df([_fuzzy_row("r1", 5.0, wildcard_level=0, origin=5)]),
        top_origin_limit=3,
        bottom_origin_limit=5,
    )
    assert len(out) == 1
    assert out.iloc[0]["source"] == "brenda"
    assert out.iloc[0]["kcat"] == 5.0


def test_fuzzy_prio3_via_wildcard_branch():
    """0 < wc <= wildcard_limit AND origin <= bottom -> prio3."""
    out = merge_dlkcat_and_fuzzy_kcats(
        _dlkcat_df([]),
        _fuzzy_df([_fuzzy_row("r1", 5.0, wildcard_level=2, origin=3)]),
        top_origin_limit=3,
        bottom_origin_limit=6,
        wildcard_limit=3,
    )
    assert len(out) == 1
    assert out.iloc[0]["source"] == "brenda"


def test_fuzzy_prio3_excluded_when_wildcard_too_high():
    out = merge_dlkcat_and_fuzzy_kcats(
        _dlkcat_df([]),
        _fuzzy_df([_fuzzy_row("r1", 5.0, wildcard_level=3, origin=1)]),
        wildcard_limit=2,
    )
    assert out.empty


def test_fuzzy_prio3_excluded_when_origin_above_bottom():
    out = merge_dlkcat_and_fuzzy_kcats(
        _dlkcat_df([]),
        _fuzzy_df([_fuzzy_row("r1", 5.0, wildcard_level=1, origin=6)]),
        bottom_origin_limit=3,
    )
    assert out.empty


def test_dlkcat_prio2_wins_over_fuzzy_prio3():
    """When the same reaction has both a prio3 fuzzy row AND a DLKcat row,
    DLKcat (prio2) wins."""
    out = merge_dlkcat_and_fuzzy_kcats(
        _dlkcat_df([_dlkcat_row("r1", 100.0)]),
        _fuzzy_df([_fuzzy_row("r1", 5.0, wildcard_level=1, origin=3)]),
        top_origin_limit=3,
        bottom_origin_limit=6,
        wildcard_limit=3,
    )
    assert len(out) == 1
    assert out.iloc[0]["source"] == "DLKcat"


# --------------------------------------------------------------------------- #
# NaN handling
# --------------------------------------------------------------------------- #

def test_fuzzy_row_with_na_origin_excluded_from_all_priorities():
    """A fuzzy row that didn't match anything (NaN origin / NaN wildcard)
    must be excluded from prio1 and prio3 alike."""
    out = merge_dlkcat_and_fuzzy_kcats(
        _dlkcat_df([_dlkcat_row("r1", 100.0)]),
        _fuzzy_df([_fuzzy_row("r1", 0.0, wildcard_level=pd.NA, origin=pd.NA)]),
    )
    assert len(out) == 1
    assert out.iloc[0]["source"] == "DLKcat"


# --------------------------------------------------------------------------- #
# Ordering and schema
# --------------------------------------------------------------------------- #

def test_output_schema_matches_inputs():
    out = merge_dlkcat_and_fuzzy_kcats(
        _dlkcat_df([_dlkcat_row("r1", 100.0)]),
        _fuzzy_df([_fuzzy_row("r2", 5.0, wildcard_level=0, origin=1)]),
    )
    assert list(out.columns) == _CANONICAL_COLUMNS


def test_fuzzy_rows_come_before_dlkcat_rows():
    out = merge_dlkcat_and_fuzzy_kcats(
        _dlkcat_df([_dlkcat_row("r1", 100.0)]),
        _fuzzy_df([_fuzzy_row("r2", 5.0, wildcard_level=0, origin=1)]),
    )
    sources = list(out["source"])
    # First all "brenda" entries, then all "DLKcat" entries.
    first_dlkcat = next((i for i, s in enumerate(sources) if s == "DLKcat"), len(sources))
    assert all(s == "brenda" for s in sources[:first_dlkcat])
    assert all(s == "DLKcat" for s in sources[first_dlkcat:])


# --------------------------------------------------------------------------- #
# Multiple rows per reaction
# --------------------------------------------------------------------------- #

def test_multiple_dlkcat_rows_per_reaction_all_kept_when_prio2():
    out = merge_dlkcat_and_fuzzy_kcats(
        _dlkcat_df([
            _dlkcat_row("r1", 100.0),
            _dlkcat_row("r1", 200.0, gene="g2"),
        ]),
        _fuzzy_df([]),
    )
    assert len(out) == 2
    assert all(out["rxn_id"] == "r1")


def test_multiple_fuzzy_rows_per_reaction_all_kept_when_prio1():
    out = merge_dlkcat_and_fuzzy_kcats(
        _dlkcat_df([]),
        _fuzzy_df([
            _fuzzy_row("r1", 5.0, wildcard_level=0, origin=1, eccode="1.1.1.1"),
            _fuzzy_row("r1", 8.0, wildcard_level=0, origin=1, eccode="2.7.7.7"),
        ]),
    )
    assert len(out) == 2
    assert all(out["rxn_id"] == "r1")
    assert sorted(out["kcat"]) == [5.0, 8.0]


def test_fuzzy_prio1_row_blocks_dlkcat_rows_for_same_reaction():
    """One fuzzy prio1 row for r1 makes ALL DLKcat rows for r1
    drop out, even though there are multiple."""
    out = merge_dlkcat_and_fuzzy_kcats(
        _dlkcat_df([
            _dlkcat_row("r1", 100.0),
            _dlkcat_row("r1", 200.0, gene="g2"),
        ]),
        _fuzzy_df([_fuzzy_row("r1", 5.0, wildcard_level=0, origin=1)]),
    )
    assert len(out) == 1
    assert out.iloc[0]["source"] == "brenda"


# --------------------------------------------------------------------------- #
# Realistic end-to-end with apply_kcat_list
# --------------------------------------------------------------------------- #

def test_output_can_feed_apply_kcat_list(tmp_path):
    """End-to-end: a merged DataFrame should pass through
    apply_kcat_list without column adjustments."""
    import numpy as np
    from scipy import sparse

    from geckopy.ec_model import EcModel
    from geckopy.ec_model.ec_data import EcData
    from geckopy.gather_kcats import apply_kcat_list

    merged = merge_dlkcat_and_fuzzy_kcats(
        _dlkcat_df([_dlkcat_row("r2", 100.0)]),
        _fuzzy_df([_fuzzy_row("r1", 5.0, wildcard_level=0, origin=1)]),
    )

    model = EcModel("test")
    model.ec = EcData(
        rxns=["r1", "r2"],
        kcat=np.array([np.nan, np.nan]),
        source=["", ""],
        notes=["", ""],
        eccodes=["", ""],
        rxn_enz_mat=sparse.csr_matrix((2, 0), dtype=float),
    )
    updated = apply_kcat_list(model, merged)
    assert sorted(updated) == ["r1", "r2"]
    assert model.ec.source[0] == "brenda"
    assert model.ec.source[1] == "DLKcat"
