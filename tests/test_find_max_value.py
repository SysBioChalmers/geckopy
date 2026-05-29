"""Tests for find_max_value."""
import pandas as pd

from geckopy.databases import BrendaData
from geckopy.kcat_sensitivity_analysis import find_max_value


# --------------------------------------------------------------------------- #
# Fixture builders
# --------------------------------------------------------------------------- #

def _brenda(
    kcat_rows: list[tuple[str, str, str, float]] | None = None,
    sa_rows: list[tuple[str, str, float, float]] | None = None,
) -> BrendaData:
    """kcat_rows: (ec_code, substrate, organism, kcat).
    sa_rows: (ec_code, organism, kcat (= SA*MW already), mw)."""
    kcat_rows = kcat_rows or []
    sa_rows = sa_rows or []
    # The same rows fill both the max and median views for these tests;
    # find_max_value only consults the max view (via the .kcat back-compat
    # alias) so the median rows are inert.
    kcat_df = pd.DataFrame(
        kcat_rows, columns=["ec_code", "substrate", "organism", "kcat"]
    )
    sa_df = pd.DataFrame(
        sa_rows, columns=["ec_code", "organism", "kcat", "mw"]
    )
    return BrendaData(
        kcat_max=kcat_df, kcat_median=kcat_df,
        sa_max=sa_df, sa_median=sa_df,
    )


# --------------------------------------------------------------------------- #
# Trivial cases
# --------------------------------------------------------------------------- #

def test_empty_string_returns_empty_result():
    val, org, param = find_max_value("", _brenda())
    assert val == 0.0
    assert org == ""
    assert param == ""


def test_no_match_in_either_table_returns_empty():
    val, org, param = find_max_value(
        "EC9.9.9.9",
        _brenda(kcat_rows=[("1.1.1.1", "*", "yeast", 5.0)]),
    )
    assert val == 0.0
    assert org == ""
    assert param == ""


def test_only_kcat_match():
    val, org, param = find_max_value(
        "EC1.1.1.1",
        _brenda(kcat_rows=[("1.1.1.1", "*", "yeast", 5.0)]),
    )
    assert val == 5.0
    assert org == "yeast"
    assert param == "K_cat"


def test_only_sa_match():
    val, org, param = find_max_value(
        "EC1.1.1.1",
        _brenda(sa_rows=[("1.1.1.1", "yeast", 8.0, 50.0)]),
    )
    assert val == 8.0
    assert org == "yeast"
    assert param == "SA*Mw"


def test_kcat_wins_over_lower_sa():
    val, _, param = find_max_value(
        "EC1.1.1.1",
        _brenda(
            kcat_rows=[("1.1.1.1", "*", "yeast", 10.0)],
            sa_rows=[("1.1.1.1", "yeast", 5.0, 50.0)],
        ),
    )
    assert val == 10.0
    assert param == "K_cat"


def test_sa_wins_over_lower_kcat():
    val, _, param = find_max_value(
        "EC1.1.1.1",
        _brenda(
            kcat_rows=[("1.1.1.1", "*", "yeast", 5.0)],
            sa_rows=[("1.1.1.1", "yeast", 10.0, 50.0)],
        ),
    )
    assert val == 10.0
    assert param == "SA*Mw"


# --------------------------------------------------------------------------- #
# Multiple ECs
# --------------------------------------------------------------------------- #

def test_multiple_ecs_picks_overall_max():
    val, org, _ = find_max_value(
        "EC1.1.1.1 EC2.7.7.7",
        _brenda(kcat_rows=[
            ("1.1.1.1", "*", "yeast", 5.0),
            ("2.7.7.7", "*", "ecoli", 100.0),
        ]),
    )
    assert val == 100.0
    assert org == "ecoli"


def test_multiple_ecs_some_missing():
    val, org, _ = find_max_value(
        "EC9.9.9.9 EC1.1.1.1",
        _brenda(kcat_rows=[("1.1.1.1", "*", "yeast", 5.0)]),
    )
    assert val == 5.0
    assert org == "yeast"


# --------------------------------------------------------------------------- #
# EC prefix tolerance
# --------------------------------------------------------------------------- #

def test_token_with_ec_prefix_matches():
    val, _, _ = find_max_value(
        "EC1.1.1.1",
        _brenda(kcat_rows=[("1.1.1.1", "*", "yeast", 5.0)]),
    )
    assert val == 5.0


def test_token_without_ec_prefix_also_matches():
    """Tolerate the prefix-stripped form."""
    val, _, _ = find_max_value(
        "1.1.1.1",
        _brenda(kcat_rows=[("1.1.1.1", "*", "yeast", 5.0)]),
    )
    assert val == 5.0


def test_lowercase_ec_prefix_tolerated():
    val, _, _ = find_max_value(
        "ec1.1.1.1",
        _brenda(kcat_rows=[("1.1.1.1", "*", "yeast", 5.0)]),
    )
    assert val == 5.0


# --------------------------------------------------------------------------- #
# Wildcard prefix matching (geckopy implementation of MATLAB intent)
# --------------------------------------------------------------------------- #

def test_wildcard_prefix_matches_all_subclasses():
    val, _, _ = find_max_value(
        "EC1.1.1.-",
        _brenda(kcat_rows=[
            ("1.1.1.1", "*", "yeast", 5.0),
            ("1.1.1.2", "*", "ecoli", 99.0),
            ("1.1.2.1", "*", "fugu", 1000.0),  # NOT under "1.1.1."
        ]),
    )
    # Only 1.1.1.1 and 1.1.1.2 qualify; max is 99.0.
    assert val == 99.0


def test_full_wildcard_matches_everything():
    val, _, _ = find_max_value(
        "EC-.-.-.-",
        _brenda(kcat_rows=[
            ("1.1.1.1", "*", "yeast", 5.0),
            ("9.9.9.9", "*", "moo", 1000.0),
        ]),
    )
    assert val == 1000.0


def test_wildcard_searches_both_tables():
    val, _, param = find_max_value(
        "EC1.1.-.-",
        _brenda(
            kcat_rows=[("1.1.1.1", "*", "yeast", 5.0)],
            sa_rows=[("1.1.2.1", "ecoli", 100.0, 50.0)],
        ),
    )
    assert val == 100.0
    assert param == "SA*Mw"


# --------------------------------------------------------------------------- #
# Multi-row picking (max within a single matched set)
# --------------------------------------------------------------------------- #

def test_multiple_rows_for_same_ec_picks_max():
    val, org, _ = find_max_value(
        "EC1.1.1.1",
        _brenda(kcat_rows=[
            ("1.1.1.1", "*", "low_org", 1.0),
            ("1.1.1.1", "*", "high_org", 50.0),
            ("1.1.1.1", "*", "mid_org", 10.0),
        ]),
    )
    assert val == 50.0
    assert org == "high_org"


# --------------------------------------------------------------------------- #
# Empty BRENDA tables
# --------------------------------------------------------------------------- #

def test_empty_brenda_returns_zero():
    val, _, _ = find_max_value("EC1.1.1.1", _brenda())
    assert val == 0.0
