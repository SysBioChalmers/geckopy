"""Tests for load_brenda_data."""
import logging
from pathlib import Path

import pytest

from geckopy.databases import BrendaData, load_brenda_data


# Header lines that load_brenda_data is expected to skip.
_COMMENT_KCAT = "# BRENDA release test generated 2026-05-18 - CC BY 4.0 - kcat in 1/s\n"
_COMMENT_SA = (
    "# BRENDA release test generated 2026-05-18 - CC BY 4.0 - "
    "specific activity in umol/min/mg\n"
)
_COMMENT_MW = (
    "# BRENDA release test generated 2026-05-18 - CC BY 4.0 - "
    "molecular weight in g/mol\n"
)
_HEADER_KCAT = (
    "ec_code\tsubstrate\torganism\tkcat_max\tkcat_median\tn\treferences\n"
)
_HEADER_SA = (
    "ec_code\tsubstrate\torganism\tsa_max\tsa_median\tn\treferences\n"
)
_HEADER_MW = "ec_code\tsubstrate\torganism\tmw\tn\treferences\n"


def _write_brenda_files(
    folder: Path,
    *,
    kcat: str = "",
    sa: str = "",
    mw: str = "",
    with_headers: bool = True,
) -> None:
    """Create the three BRENDA TSVs with optional header lines."""
    kcat_text = (_COMMENT_KCAT + _HEADER_KCAT + kcat) if with_headers else kcat
    sa_text = (_COMMENT_SA + _HEADER_SA + sa) if with_headers else sa
    mw_text = (_COMMENT_MW + _HEADER_MW + mw) if with_headers else mw
    (folder / "kcat.tsv").write_text(kcat_text, encoding="utf-8")
    (folder / "sa.tsv").write_text(sa_text, encoding="utf-8")
    (folder / "mw.tsv").write_text(mw_text, encoding="utf-8")


# --------------------------------------------------------------------------- #
# Trivial cases
# --------------------------------------------------------------------------- #

def test_returns_brenda_data_with_split_views(tmp_path):
    _write_brenda_files(tmp_path)
    result = load_brenda_data(tmp_path)
    assert isinstance(result, BrendaData)
    expected_kcat = ["ec_code", "substrate", "organism", "kcat", "n"]
    expected_sa = ["ec_code", "organism", "kcat", "mw", "n"]
    assert list(result.kcat_max.columns) == expected_kcat
    assert list(result.kcat_median.columns) == expected_kcat
    assert list(result.sa_max.columns) == expected_sa
    assert list(result.sa_median.columns) == expected_sa


def test_kcat_for_picks_view_by_aggregation(tmp_path):
    _write_brenda_files(
        tmp_path,
        kcat="1.1.1.1\tm1\torg1\t10\t4\t3\t*\n",
    )
    result = load_brenda_data(tmp_path)
    assert result.kcat_for("max").iloc[0]["kcat"] == 10.0
    assert result.kcat_for("median").iloc[0]["kcat"] == 4.0


def test_kcat_for_invalid_aggregation_raises(tmp_path):
    _write_brenda_files(tmp_path)
    result = load_brenda_data(tmp_path)
    with pytest.raises(ValueError, match="aggregation must be"):
        result.kcat_for("mean")
    with pytest.raises(ValueError, match="aggregation must be"):
        result.sa_for("bogus")


def test_all_empty_files_yield_empty_dataframes(tmp_path):
    _write_brenda_files(tmp_path)
    result = load_brenda_data(tmp_path)
    assert result.kcat_max.empty
    assert result.kcat_median.empty
    assert result.sa_max.empty
    assert result.sa_median.empty


# --------------------------------------------------------------------------- #
# Missing files
# --------------------------------------------------------------------------- #

def test_missing_kcat_file_raises(tmp_path):
    (tmp_path / "sa.tsv").write_text("")
    (tmp_path / "mw.tsv").write_text("")
    with pytest.raises(FileNotFoundError, match="kcat.tsv"):
        load_brenda_data(tmp_path)


def test_missing_sa_file_raises(tmp_path):
    (tmp_path / "kcat.tsv").write_text("")
    (tmp_path / "mw.tsv").write_text("")
    with pytest.raises(FileNotFoundError, match="sa.tsv"):
        load_brenda_data(tmp_path)


def test_missing_mw_file_raises(tmp_path):
    (tmp_path / "kcat.tsv").write_text("")
    (tmp_path / "sa.tsv").write_text("")
    with pytest.raises(FileNotFoundError, match="mw.tsv"):
        load_brenda_data(tmp_path)


def test_missing_folder_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        load_brenda_data(tmp_path / "does_not_exist")


# --------------------------------------------------------------------------- #
# KCAT file parsing
# --------------------------------------------------------------------------- #

def test_kcat_basic_parse_max_and_median(tmp_path):
    _write_brenda_files(
        tmp_path,
        kcat="1.1.1.1\tm1\torg1\t10\t10\t1\tPMID:1\n",
    )
    result = load_brenda_data(tmp_path)
    assert len(result.kcat_max) == 1
    assert len(result.kcat_median) == 1
    row_max = result.kcat_max.iloc[0]
    row_med = result.kcat_median.iloc[0]
    assert row_max["ec_code"] == "1.1.1.1"
    assert row_max["substrate"] == "m1"
    assert row_max["organism"] == "org1"
    assert row_max["kcat"] == 10.0
    assert int(row_max["n"]) == 1
    assert row_med["kcat"] == 10.0


def test_kcat_no_scaling_applied(tmp_path):
    _write_brenda_files(
        tmp_path,
        kcat="1.1.1.1\tm1\torg1\t42.5\t42.5\t1\t*\n",
    )
    result = load_brenda_data(tmp_path)
    assert result.kcat_max.iloc[0]["kcat"] == pytest.approx(42.5)
    assert result.kcat_median.iloc[0]["kcat"] == pytest.approx(42.5)


def test_kcat_views_carry_different_values_from_same_row(tmp_path):
    """A wide row with distinct max and median lands in both views with
    its respective value."""
    _write_brenda_files(
        tmp_path,
        kcat=(
            "1.1.1.1\tm1\torg1\t10\t4\t3\t*\n"
            "2.7.7.7\tm2\torg2\t100\t40\t5\t*\n"
        ),
    )
    result = load_brenda_data(tmp_path)
    assert sorted(result.kcat_max["kcat"]) == [10.0, 100.0]
    assert sorted(result.kcat_median["kcat"]) == [4.0, 40.0]


def test_header_and_comment_lines_skipped(tmp_path):
    """Both the ``#`` release comment and the TSV column-header line are
    skipped."""
    _write_brenda_files(
        tmp_path,
        kcat="1.1.1.1\tethanol\tyeast\t23.5\t23.5\t1\tPMID:1\n",
    )
    result = load_brenda_data(tmp_path)
    assert len(result.kcat_max) == 1
    assert result.kcat_max.iloc[0]["ec_code"] == "1.1.1.1"


# --------------------------------------------------------------------------- #
# Malformed lines
# --------------------------------------------------------------------------- #

def test_malformed_too_few_columns_skipped_with_warning(tmp_path, caplog):
    _write_brenda_files(
        tmp_path,
        kcat=(
            "1.1.1.1\tm1\torg1\t1\t1\t1\t*\n"
            "1.1.1.1\tm2\torg1\n"
            "1.1.1.1\tm3\torg1\t3\t3\t1\t*\n"
        ),
    )
    with caplog.at_level(logging.WARNING):
        result = load_brenda_data(tmp_path)
    assert len(result.kcat_max) == 2
    assert "tab-delimited fields" in caplog.text


def test_malformed_non_numeric_value_skipped_with_warning(tmp_path, caplog):
    _write_brenda_files(
        tmp_path,
        kcat=(
            "1.1.1.1\tm1\torg1\t1\t1\t1\t*\n"
            "1.1.1.1\tm2\torg1\tnot_a_number\t2\t1\t*\n"
            "1.1.1.1\tm3\torg1\t3\t3\t1\t*\n"
        ),
    )
    with caplog.at_level(logging.WARNING):
        result = load_brenda_data(tmp_path)
    assert len(result.kcat_max) == 2
    assert "non-numeric" in caplog.text


# --------------------------------------------------------------------------- #
# SA + MW join
# --------------------------------------------------------------------------- #

def test_sa_mw_join_computes_kcat(tmp_path):
    """SA = 60 [umol/min/mg] = 1 [mmol/s/g]; MW = 1000 [g/mol] = 1 [g/mmol].
    kcat = SA * MW = 1 [1/s]."""
    _write_brenda_files(
        tmp_path,
        sa="1.1.1.1\t*\torg1\t60\t60\t1\t*\n",
        mw="1.1.1.1\t*\torg1\t1000\t1\t*\n",
    )
    result = load_brenda_data(tmp_path)
    assert len(result.sa_max) == 1
    row = result.sa_max.iloc[0]
    assert row["ec_code"] == "1.1.1.1"
    assert row["organism"] == "org1"
    assert row["kcat"] == pytest.approx(1.0)
    assert row["mw"] == pytest.approx(1.0)


def test_sa_views_have_separate_kcat_values_from_max_and_median_columns(tmp_path):
    _write_brenda_files(
        tmp_path,
        sa="1.1.1.1\t*\torg1\t120\t60\t3\t*\n",
        mw="1.1.1.1\t*\torg1\t1000\t1\t*\n",
    )
    result = load_brenda_data(tmp_path)
    assert result.sa_max.iloc[0]["kcat"] == pytest.approx(2.0)
    assert result.sa_median.iloc[0]["kcat"] == pytest.approx(1.0)


def test_sa_without_mw_match_dropped(tmp_path):
    _write_brenda_files(
        tmp_path,
        sa="1.1.1.1\t*\torgX\t60\t60\t1\t*\n",
        mw="1.1.1.1\t*\torgY\t1000\t1\t*\n",
    )
    result = load_brenda_data(tmp_path)
    assert result.sa_max.empty
    assert result.sa_median.empty


def test_sa_mw_join_case_insensitive_on_organism(tmp_path):
    _write_brenda_files(
        tmp_path,
        sa="1.1.1.1\t*\tORGANISM ONE\t60\t60\t1\t*\n",
        mw="1.1.1.1\t*\torganism one\t1000\t1\t*\n",
    )
    result = load_brenda_data(tmp_path)
    assert len(result.sa_max) == 1
    assert len(result.sa_median) == 1


def test_sa_first_mw_match_wins(tmp_path):
    _write_brenda_files(
        tmp_path,
        sa="1.1.1.1\t*\torg1\t60\t60\t1\t*\n",
        mw=(
            "1.1.1.1\t*\torg1\t1000\t1\t*\n"
            "1.1.1.1\t*\torg1\t9000\t1\t*\n"
        ),
    )
    result = load_brenda_data(tmp_path)
    assert len(result.sa_max) == 1
    assert result.sa_max.iloc[0]["mw"] == pytest.approx(1.0)


def test_empty_sa_yields_empty_sa_table_even_with_mw(tmp_path):
    _write_brenda_files(
        tmp_path,
        sa="",
        mw="1.1.1.1\t*\torg1\t1000\t1\t*\n",
    )
    result = load_brenda_data(tmp_path)
    assert result.sa_max.empty
    assert result.sa_median.empty


def test_empty_mw_yields_empty_sa_table_even_with_sa(tmp_path):
    _write_brenda_files(
        tmp_path,
        sa="1.1.1.1\t*\torg1\t60\t60\t1\t*\n",
        mw="",
    )
    result = load_brenda_data(tmp_path)
    assert result.sa_max.empty
    assert result.sa_median.empty


# --------------------------------------------------------------------------- #
# Real fixture from examples/ecTestGEM
# --------------------------------------------------------------------------- #

@pytest.fixture
def ec_test_gem_brenda_folder() -> Path:
    return Path(__file__).resolve().parent.parent / "examples" / "ecTestGEM" / "data"


def test_real_fixture_loads_kcat(ec_test_gem_brenda_folder):
    result = load_brenda_data(ec_test_gem_brenda_folder)
    # The fixture has 3 (ec, substrate, organism) triples each with a
    # single source measurement, so max == median per triple.
    assert len(result.kcat_max) == 3
    assert list(result.kcat_max["ec_code"]) == ["1.1.1.1", "1.1.1.1", "1.1.2.2"]
    assert list(result.kcat_max["substrate"]) == ["m1", "m2", "m1"]
    assert list(result.kcat_max["kcat"]) == [1.0, 10.0, 100.0]
    assert list(result.kcat_median["kcat"]) == [1.0, 10.0, 100.0]


def test_real_fixture_sa_table_empty_due_to_mismatched_organisms(
    ec_test_gem_brenda_folder,
):
    """The fixture has SA for `acetobacter pasteurianus` and MW for
    `testus testus` - different organisms, so no SA row joins."""
    result = load_brenda_data(ec_test_gem_brenda_folder)
    assert result.sa_max.empty
    assert result.sa_median.empty
