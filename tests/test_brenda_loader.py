"""Tests for load_brenda_data."""
import logging
from pathlib import Path

import pytest

from geckopy.databases import BrendaData, load_brenda_data


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #

def _write_brenda_files(
    folder: Path,
    *,
    kcat: str = "",
    sa: str = "",
    mw: str = "",
) -> None:
    """Create the three BRENDA files in the given folder with the given content."""
    (folder / "max_KCAT.txt").write_text(kcat, encoding="utf-8")
    (folder / "max_SA.txt").write_text(sa, encoding="utf-8")
    (folder / "max_MW.txt").write_text(mw, encoding="utf-8")


# --------------------------------------------------------------------------- #
# Trivial cases
# --------------------------------------------------------------------------- #

def test_returns_brenda_data_with_two_dataframes(tmp_path):
    _write_brenda_files(tmp_path)
    result = load_brenda_data(tmp_path)
    assert isinstance(result, BrendaData)
    assert list(result.kcat.columns) == ["ec_code", "substrate", "organism", "kcat"]
    assert list(result.sa.columns) == ["ec_code", "organism", "kcat", "mw"]


def test_all_empty_files_yield_empty_dataframes(tmp_path):
    _write_brenda_files(tmp_path)
    result = load_brenda_data(tmp_path)
    assert result.kcat.empty
    assert result.sa.empty


# --------------------------------------------------------------------------- #
# Missing files
# --------------------------------------------------------------------------- #

def test_missing_kcat_file_raises(tmp_path):
    (tmp_path / "max_SA.txt").write_text("")
    (tmp_path / "max_MW.txt").write_text("")
    with pytest.raises(FileNotFoundError, match="max_KCAT.txt"):
        load_brenda_data(tmp_path)


def test_missing_sa_file_raises(tmp_path):
    (tmp_path / "max_KCAT.txt").write_text("")
    (tmp_path / "max_MW.txt").write_text("")
    with pytest.raises(FileNotFoundError, match="max_SA.txt"):
        load_brenda_data(tmp_path)


def test_missing_mw_file_raises(tmp_path):
    (tmp_path / "max_KCAT.txt").write_text("")
    (tmp_path / "max_SA.txt").write_text("")
    with pytest.raises(FileNotFoundError, match="max_MW.txt"):
        load_brenda_data(tmp_path)


def test_missing_folder_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        load_brenda_data(tmp_path / "does_not_exist")


# --------------------------------------------------------------------------- #
# KCAT file parsing
# --------------------------------------------------------------------------- #

def test_kcat_basic_parse(tmp_path):
    _write_brenda_files(
        tmp_path,
        kcat="EC1.1.1.1\tm1\torg1//taxonomy//kegg\t10\t*\n",
    )
    result = load_brenda_data(tmp_path)
    assert len(result.kcat) == 1
    row = result.kcat.iloc[0]
    assert row["ec_code"] == "1.1.1.1"
    assert row["substrate"] == "m1"
    assert row["organism"] == "org1"
    assert row["kcat"] == 10.0


def test_kcat_no_scaling_applied(tmp_path):
    """KCAT values pass through unchanged."""
    _write_brenda_files(
        tmp_path,
        kcat="EC1.1.1.1\tm1\torg1\t42.5\t*\n",
    )
    result = load_brenda_data(tmp_path)
    assert result.kcat.iloc[0]["kcat"] == pytest.approx(42.5)


def test_kcat_multiple_rows(tmp_path):
    _write_brenda_files(
        tmp_path,
        kcat=(
            "EC1.1.1.1\tm1\torg1\t1\t*\n"
            "EC1.1.1.1\tm2\torg2\t2\t*\n"
            "EC2.7.7.7\tm3\torg1\t3\t*\n"
        ),
    )
    result = load_brenda_data(tmp_path)
    assert len(result.kcat) == 3
    assert list(result.kcat["ec_code"]) == ["1.1.1.1", "1.1.1.1", "2.7.7.7"]
    assert list(result.kcat["substrate"]) == ["m1", "m2", "m3"]
    assert list(result.kcat["kcat"]) == [1.0, 2.0, 3.0]


def test_kcat_organism_taxonomy_stripped(tmp_path):
    _write_brenda_files(
        tmp_path,
        kcat="EC1.1.1.1\tm1\thomo sapiens//eukaryota;animalia//hsa\t1\t*\n",
    )
    result = load_brenda_data(tmp_path)
    assert result.kcat.iloc[0]["organism"] == "homo sapiens"


def test_kcat_organism_without_double_slash(tmp_path):
    """Organism with no `//` is unchanged."""
    _write_brenda_files(
        tmp_path,
        kcat="EC1.1.1.1\tm1\tjust an organism\t1\t*\n",
    )
    result = load_brenda_data(tmp_path)
    assert result.kcat.iloc[0]["organism"] == "just an organism"


def test_fifth_column_dropped(tmp_path):
    """5th column (always `*` in real dumps) does not appear in output."""
    _write_brenda_files(
        tmp_path,
        kcat="EC1.1.1.1\tm1\torg1\t1\tsome_metadata\n",
    )
    result = load_brenda_data(tmp_path)
    assert list(result.kcat.columns) == ["ec_code", "substrate", "organism", "kcat"]


def test_blank_lines_skipped_silently(tmp_path, caplog):
    _write_brenda_files(
        tmp_path,
        kcat=(
            "EC1.1.1.1\tm1\torg1\t1\t*\n"
            "\n"
            "EC1.1.1.1\tm2\torg1\t2\t*\n"
        ),
    )
    with caplog.at_level(logging.WARNING):
        result = load_brenda_data(tmp_path)
    assert len(result.kcat) == 2
    assert "skipped" not in caplog.text


# --------------------------------------------------------------------------- #
# Malformed lines: warn and skip
# --------------------------------------------------------------------------- #

def test_malformed_too_few_columns_skipped_with_warning(tmp_path, caplog):
    _write_brenda_files(
        tmp_path,
        kcat=(
            "EC1.1.1.1\tm1\torg1\t1\t*\n"
            "EC1.1.1.1\tm2\torg1\n"  # only 3 columns
            "EC1.1.1.1\tm3\torg1\t3\t*\n"
        ),
    )
    with caplog.at_level(logging.WARNING):
        result = load_brenda_data(tmp_path)
    assert len(result.kcat) == 2
    assert "expected 5 tab-delimited fields" in caplog.text


def test_malformed_non_numeric_value_skipped_with_warning(tmp_path, caplog):
    _write_brenda_files(
        tmp_path,
        kcat=(
            "EC1.1.1.1\tm1\torg1\t1\t*\n"
            "EC1.1.1.1\tm2\torg1\tnot_a_number\t*\n"
            "EC1.1.1.1\tm3\torg1\t3\t*\n"
        ),
    )
    with caplog.at_level(logging.WARNING):
        result = load_brenda_data(tmp_path)
    assert len(result.kcat) == 2
    assert "is not numeric" in caplog.text


# --------------------------------------------------------------------------- #
# SA + MW join: unit conversions and lookup
# --------------------------------------------------------------------------- #

def test_sa_mw_join_computes_kcat(tmp_path):
    """SA = 60 [umol/min/mg] = 1 [mmol/s/g]; MW = 1000 [g/mol] = 1 [g/mmol].
    kcat = SA * MW = 1 * 1 = 1 [1/s]."""
    _write_brenda_files(
        tmp_path,
        sa="EC1.1.1.1\t*\torg1\t60\t*\n",
        mw="EC1.1.1.1\t*\torg1\t1000\t*\n",
    )
    result = load_brenda_data(tmp_path)
    assert len(result.sa) == 1
    row = result.sa.iloc[0]
    assert row["ec_code"] == "1.1.1.1"
    assert row["organism"] == "org1"
    assert row["kcat"] == pytest.approx(1.0)
    assert row["mw"] == pytest.approx(1.0)


def test_sa_scaling_factor_one_over_sixty(tmp_path):
    _write_brenda_files(
        tmp_path,
        sa="EC1.1.1.1\t*\torg1\t120\t*\n",  # 120 / 60 = 2
        mw="EC1.1.1.1\t*\torg1\t1000\t*\n",  # 1000 / 1000 = 1
    )
    result = load_brenda_data(tmp_path)
    # SA = 2 [mmol/s/g], MW = 1 [g/mmol] -> kcat = 2 [1/s]
    assert result.sa.iloc[0]["kcat"] == pytest.approx(2.0)
    assert result.sa.iloc[0]["mw"] == pytest.approx(1.0)


def test_mw_scaling_factor_one_over_thousand(tmp_path):
    _write_brenda_files(
        tmp_path,
        sa="EC1.1.1.1\t*\torg1\t60\t*\n",  # -> 1 mmol/s/g
        mw="EC1.1.1.1\t*\torg1\t50000\t*\n",  # 50000 / 1000 = 50 g/mmol
    )
    result = load_brenda_data(tmp_path)
    assert result.sa.iloc[0]["kcat"] == pytest.approx(50.0)
    assert result.sa.iloc[0]["mw"] == pytest.approx(50.0)


def test_sa_without_mw_match_dropped(tmp_path):
    _write_brenda_files(
        tmp_path,
        sa="EC1.1.1.1\t*\torgX\t60\t*\n",
        mw="EC1.1.1.1\t*\torgY\t1000\t*\n",  # different organism
    )
    result = load_brenda_data(tmp_path)
    assert result.sa.empty


def test_sa_without_mw_for_ec_dropped(tmp_path):
    _write_brenda_files(
        tmp_path,
        sa="EC1.1.1.1\t*\torg1\t60\t*\n",
        mw="EC2.7.7.7\t*\torg1\t1000\t*\n",  # different EC
    )
    result = load_brenda_data(tmp_path)
    assert result.sa.empty


def test_sa_mw_join_case_insensitive_on_ec(tmp_path):
    _write_brenda_files(
        tmp_path,
        sa="EC1.1.1.1\t*\torg1\t60\t*\n",
        mw="ec1.1.1.1\t*\torg1\t1000\t*\n",
    )
    result = load_brenda_data(tmp_path)
    assert len(result.sa) == 1


def test_sa_mw_join_case_insensitive_on_organism(tmp_path):
    _write_brenda_files(
        tmp_path,
        sa="EC1.1.1.1\t*\tORGANISM ONE\t60\t*\n",
        mw="EC1.1.1.1\t*\torganism one\t1000\t*\n",
    )
    result = load_brenda_data(tmp_path)
    assert len(result.sa) == 1


def test_sa_first_mw_match_wins(tmp_path):
    _write_brenda_files(
        tmp_path,
        sa="EC1.1.1.1\t*\torg1\t60\t*\n",
        mw=(
            "EC1.1.1.1\t*\torg1\t1000\t*\n"  # first match: 1 g/mmol
            "EC1.1.1.1\t*\torg1\t9000\t*\n"  # second match: ignored
        ),
    )
    result = load_brenda_data(tmp_path)
    assert len(result.sa) == 1
    assert result.sa.iloc[0]["mw"] == pytest.approx(1.0)


def test_sa_organism_taxonomy_stripped_before_join(tmp_path):
    """Both SA and MW have //-suffixed organism strings; the join uses
    the trimmed form."""
    _write_brenda_files(
        tmp_path,
        sa="EC1.1.1.1\t*\torg1//tax1//kegg1\t60\t*\n",
        mw="EC1.1.1.1\t*\torg1//tax2//kegg2\t1000\t*\n",
    )
    result = load_brenda_data(tmp_path)
    assert len(result.sa) == 1
    assert result.sa.iloc[0]["organism"] == "org1"


def test_multiple_sa_rows_each_join_independently(tmp_path):
    _write_brenda_files(
        tmp_path,
        sa=(
            "EC1.1.1.1\t*\torg1\t60\t*\n"
            "EC1.1.1.1\t*\torg2\t120\t*\n"
            "EC2.7.7.7\t*\torg1\t60\t*\n"
        ),
        mw=(
            "EC1.1.1.1\t*\torg1\t1000\t*\n"
            "EC1.1.1.1\t*\torg2\t2000\t*\n"
            "EC2.7.7.7\t*\torg1\t1000\t*\n"
        ),
    )
    result = load_brenda_data(tmp_path)
    assert len(result.sa) == 3
    # All three SA rows had a matching MW row.


def test_empty_sa_yields_empty_sa_table_even_with_mw(tmp_path):
    _write_brenda_files(
        tmp_path,
        sa="",
        mw="EC1.1.1.1\t*\torg1\t1000\t*\n",
    )
    result = load_brenda_data(tmp_path)
    assert result.sa.empty


def test_empty_mw_yields_empty_sa_table_even_with_sa(tmp_path):
    _write_brenda_files(
        tmp_path,
        sa="EC1.1.1.1\t*\torg1\t60\t*\n",
        mw="",
    )
    result = load_brenda_data(tmp_path)
    assert result.sa.empty


# --------------------------------------------------------------------------- #
# Real fixture from examples/ecTestGEM
# --------------------------------------------------------------------------- #

@pytest.fixture
def ec_test_gem_brenda_folder() -> Path:
    return Path(__file__).resolve().parent.parent / "examples" / "ecTestGEM" / "data"


def test_real_fixture_loads_kcat(ec_test_gem_brenda_folder):
    result = load_brenda_data(ec_test_gem_brenda_folder)
    assert len(result.kcat) == 3
    # Expected EC codes from the fixture (EC prefix stripped).
    assert list(result.kcat["ec_code"]) == ["1.1.1.1", "1.1.1.1", "1.1.2.2"]
    assert list(result.kcat["substrate"]) == ["m1", "m2", "m1"]
    assert list(result.kcat["kcat"]) == [1.0, 10.0, 100.0]


def test_real_fixture_sa_table_empty_due_to_mismatched_organisms(
    ec_test_gem_brenda_folder,
):
    """The fixture has SA for `acetobacter pasteurianus` and MW for
    `testus testus` — different organisms, so no SA row joins."""
    result = load_brenda_data(ec_test_gem_brenda_folder)
    assert result.sa.empty


def test_real_fixture_organism_taxonomy_stripped(ec_test_gem_brenda_folder):
    result = load_brenda_data(ec_test_gem_brenda_folder)
    # Organism strings in the fixture KCAT file are like
    # `testus testus//*//*`; after trimming the kcat rows have plain organism.
    organisms = list(result.kcat["organism"])
    assert all("//" not in o for o in organisms)
