"""Tests for load_dlkcat_ignore_lists."""

from geckopy.databases import DLKcatIgnoreLists, load_dlkcat_ignore_lists
from geckopy.databases.dlkcat_ignore_lists import _normalize


# --------------------------------------------------------------------------- #
# _normalize helper
# --------------------------------------------------------------------------- #

def test_normalize_lowercases_and_strips_special():
    assert _normalize("H2O") == "h2o"
    assert _normalize("H+") == "h"
    assert _normalize("Mg(2+)") == "mg2"
    assert _normalize("Cu2(+)") == "cu2"
    assert _normalize("ATP") == "atp"


def test_normalize_strips_whitespace():
    assert _normalize("hydrogen sulfide") == "hydrogensulfide"


def test_normalize_handles_empty():
    assert _normalize("") == ""
    assert _normalize("()") == ""


# --------------------------------------------------------------------------- #
# Defaults shipped with geckopy
# --------------------------------------------------------------------------- #

def test_load_defaults_returns_populated_lists():
    """When no project folder is given, the shipped defaults are loaded."""
    lists = load_dlkcat_ignore_lists()
    assert isinstance(lists, DLKcatIgnoreLists)
    assert len(lists.ignore_names) > 0
    assert len(lists.ignore_smiles) > 0
    assert len(lists.currency_pairs) > 0


def test_default_ignore_names_include_water():
    lists = load_dlkcat_ignore_lists()
    assert "h2o" in lists.ignore_names


def test_default_currency_pairs_include_atp_adp():
    lists = load_dlkcat_ignore_lists()
    assert ("atp", "adp") in lists.currency_pairs


def test_default_smiles_include_water_smiles():
    lists = load_dlkcat_ignore_lists()
    # "H2O\tO" -> SMILES "O".
    assert "O" in lists.ignore_smiles


# --------------------------------------------------------------------------- #
# Project-folder override
# --------------------------------------------------------------------------- #

def test_project_folder_overrides_ignore_file(tmp_path):
    (tmp_path / "DLKcatIgnoreMets.tsv").write_text(
        "custom_met\tCUSTOMSMILES\n", encoding="utf-8",
    )
    lists = load_dlkcat_ignore_lists(tmp_path)
    assert lists.ignore_names == ["custommet"]
    assert lists.ignore_smiles == ["CUSTOMSMILES"]
    # Currency pairs still come from defaults since only ignore was overridden.
    assert ("atp", "adp") in lists.currency_pairs


def test_project_folder_overrides_currency_file(tmp_path):
    (tmp_path / "DLKcatCurrencyMets.tsv").write_text(
        "first\tsecond\n", encoding="utf-8",
    )
    lists = load_dlkcat_ignore_lists(tmp_path)
    assert lists.currency_pairs == [("first", "second")]
    # Ignore lists still come from defaults.
    assert "h2o" in lists.ignore_names


def test_missing_project_folder_uses_defaults(tmp_path):
    # tmp_path exists but contains no .tsv files.
    lists = load_dlkcat_ignore_lists(tmp_path)
    assert "h2o" in lists.ignore_names
    assert ("atp", "adp") in lists.currency_pairs


def test_project_folder_none_uses_defaults():
    lists = load_dlkcat_ignore_lists(None)
    assert "h2o" in lists.ignore_names


# --------------------------------------------------------------------------- #
# Parsing edge cases
# --------------------------------------------------------------------------- #

def test_blank_lines_skipped(tmp_path):
    (tmp_path / "DLKcatIgnoreMets.tsv").write_text(
        "name1\tSMI1\n\nname2\tSMI2\n", encoding="utf-8",
    )
    lists = load_dlkcat_ignore_lists(tmp_path)
    assert lists.ignore_names == ["name1", "name2"]


def test_one_column_ignore_line_skipped(tmp_path):
    """The ignore file requires both columns; a one-column line is skipped."""
    (tmp_path / "DLKcatIgnoreMets.tsv").write_text(
        "name1\tSMI1\nbroken_line\nname2\tSMI2\n", encoding="utf-8",
    )
    lists = load_dlkcat_ignore_lists(tmp_path)
    assert lists.ignore_names == ["name1", "name2"]


def test_empty_smiles_field_omitted_from_ignore_smiles(tmp_path):
    """A row with empty SMILES still contributes to ignore_names but
    does not add an empty string to ignore_smiles."""
    (tmp_path / "DLKcatIgnoreMets.tsv").write_text(
        "name1\tSMI1\nname2\t\n", encoding="utf-8",
    )
    lists = load_dlkcat_ignore_lists(tmp_path)
    assert lists.ignore_names == ["name1", "name2"]
    assert lists.ignore_smiles == ["SMI1"]


def test_currency_pair_with_empty_field_skipped(tmp_path):
    (tmp_path / "DLKcatCurrencyMets.tsv").write_text(
        "a\tb\n\tc\nd\t\n", encoding="utf-8",
    )
    lists = load_dlkcat_ignore_lists(tmp_path)
    assert lists.currency_pairs == [("a", "b")]


def test_currency_pair_normalization(tmp_path):
    (tmp_path / "DLKcatCurrencyMets.tsv").write_text(
        "Mg(2+)\tH2O\n", encoding="utf-8",
    )
    lists = load_dlkcat_ignore_lists(tmp_path)
    assert lists.currency_pairs == [("mg2", "h2o")]
