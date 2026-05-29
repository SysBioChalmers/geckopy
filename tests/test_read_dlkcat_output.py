"""Tests for read_dlkcat_output."""
from pathlib import Path

import cobra
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from geckopy.ec_model import EcModel
from geckopy.ec_model.ec_data import EcData
from geckopy.gather_kcats import read_dlkcat_output


# --------------------------------------------------------------------------- #
# Fixture builders
# --------------------------------------------------------------------------- #

def _ec_model(
    ec_rxns: list[str],
    metabolite_names: list[str],
) -> EcModel:
    """Build an EcModel with the given ec.rxns and metabolites (by name)."""
    model = EcModel("test")
    for i, name in enumerate(metabolite_names):
        m = cobra.Metabolite(f"m{i}", compartment="c")
        m.name = name
        model.add_metabolites([m])

    n = len(ec_rxns)
    model.ec = EcData(
        rxns=list(ec_rxns),
        kcat=np.full(n, np.nan, dtype=float),
        source=[""] * n,
        notes=[""] * n,
        eccodes=[""] * n,
        rxn_enz_mat=sparse.csr_matrix((n, 0), dtype=float),
    )
    return model


def _write_dlkcat_output(tmp_path: Path, rows: list[tuple]) -> Path:
    """Write a DLKcat output TSV with the standard header.

    Each row is (rxn_id, gene, substrate, smiles, sequence, kcat).
    """
    p = tmp_path / "DLKcat.tsv"
    with open(p, "w", encoding="utf-8") as f:
        f.write("Reaction\tGene\tSubstrate\tSMILES\tSequence\tKcat\n")
        for row in rows:
            f.write("\t".join(str(c) for c in row) + "\n")
    return p


# --------------------------------------------------------------------------- #
# Trivial cases
# --------------------------------------------------------------------------- #

def test_missing_file_raises(tmp_path):
    model = _ec_model(["r1"], ["alpha"])
    with pytest.raises(FileNotFoundError, match="DLKcat output file"):
        read_dlkcat_output(model, tmp_path / "missing.tsv")


def test_all_na_kcats_raises(tmp_path):
    model = _ec_model(["r1"], ["alpha"])
    p = _write_dlkcat_output(tmp_path, [
        ("r1", "g1", "alpha", "CC", "MASEQ", "NA"),
        ("r1", "g1", "alpha", "CC", "MASEQ", "NA"),
    ])
    with pytest.raises(ValueError, match="no numeric kcat values"):
        read_dlkcat_output(model, p)


def test_empty_kcats_raises(tmp_path):
    model = _ec_model(["r1"], ["alpha"])
    p = _write_dlkcat_output(tmp_path, [
        ("r1", "g1", "alpha", "CC", "MASEQ", ""),
    ])
    with pytest.raises(ValueError, match="no numeric kcat values"):
        read_dlkcat_output(model, p)


# --------------------------------------------------------------------------- #
# Happy path
# --------------------------------------------------------------------------- #

def test_single_row_basic_parse(tmp_path):
    model = _ec_model(["r1"], ["alpha"])
    p = _write_dlkcat_output(tmp_path, [
        ("r1", "g1", "alpha", "CC", "MASEQ", "5.0"),
    ])
    df = read_dlkcat_output(model, p)
    assert len(df) == 1
    row = df.iloc[0]
    assert row["rxn_id"] == "r1"
    assert row["source"] == "DLKcat"
    assert row["substrates"] == ["alpha"]
    assert row["genes"] == ["g1"]
    assert row["kcat"] == pytest.approx(5.0)
    assert row["eccode"] == ""
    assert pd.isna(row["wildcard_level"])
    assert pd.isna(row["origin"])


def test_multi_row_parse(tmp_path):
    model = _ec_model(["r1", "r2", "r3"], ["alpha", "beta", "gamma"])
    p = _write_dlkcat_output(tmp_path, [
        ("r1", "g1", "alpha", "CC", "M1", "5.0"),
        ("r2", "g2", "beta", "CN", "M2", "7.5"),
        ("r3", "g3", "gamma", "CO", "M3", "9.0"),
    ])
    df = read_dlkcat_output(model, p)
    assert len(df) == 3
    assert list(df["rxn_id"]) == ["r1", "r2", "r3"]
    assert list(df["kcat"]) == [5.0, 7.5, 9.0]


def test_substrates_and_genes_are_single_element_lists(tmp_path):
    model = _ec_model(["r1"], ["alpha"])
    p = _write_dlkcat_output(tmp_path, [
        ("r1", "g1", "alpha", "CC", "MASEQ", "5.0"),
    ])
    df = read_dlkcat_output(model, p)
    assert isinstance(df.iloc[0]["substrates"], list)
    assert isinstance(df.iloc[0]["genes"], list)
    assert len(df.iloc[0]["substrates"]) == 1
    assert len(df.iloc[0]["genes"]) == 1


def test_schema_matches_fuzzy_output(tmp_path):
    """The output column set must match what fuzzy_kcat_matching produces
    so they can be merged downstream."""
    model = _ec_model(["r1"], ["alpha"])
    p = _write_dlkcat_output(tmp_path, [
        ("r1", "g1", "alpha", "CC", "MASEQ", "5.0"),
    ])
    df = read_dlkcat_output(model, p)
    assert list(df.columns) == [
        "rxn_id", "source", "eccode", "substrates", "genes", "kcat",
        "wildcard_level", "origin",
    ]


# --------------------------------------------------------------------------- #
# Filtering: rows with non-numeric kcat are dropped silently
# --------------------------------------------------------------------------- #

def test_mixed_valid_and_na_kept_only_valid(tmp_path, caplog):
    import logging
    model = _ec_model(["r1", "r2", "r3"], ["alpha", "beta", "gamma"])
    p = _write_dlkcat_output(tmp_path, [
        ("r1", "g1", "alpha", "CC", "M1", "5.0"),
        ("r2", "g2", "beta", "CN", "M2", "NA"),
        ("r3", "g3", "gamma", "CO", "M3", "9.0"),
    ])
    with caplog.at_level(logging.INFO):
        df = read_dlkcat_output(model, p)
    assert len(df) == 2
    assert list(df["rxn_id"]) == ["r1", "r3"]
    assert "dropped 1" in caplog.text


# --------------------------------------------------------------------------- #
# Case-insensitive substrate matching (geckopy divergence from MATLAB)
# --------------------------------------------------------------------------- #

def test_substrate_match_is_case_insensitive(tmp_path):
    model = _ec_model(["r1"], ["Alpha"])
    p = _write_dlkcat_output(tmp_path, [
        ("r1", "g1", "ALPHA", "CC", "MASEQ", "5.0"),
    ])
    df = read_dlkcat_output(model, p)
    assert len(df) == 1
    # The output preserves the file's casing for the substrate.
    assert df.iloc[0]["substrates"] == ["ALPHA"]


def test_unknown_substrate_raises(tmp_path):
    model = _ec_model(["r1"], ["alpha"])
    p = _write_dlkcat_output(tmp_path, [
        ("r1", "g1", "completely_unknown", "CC", "M", "5.0"),
    ])
    with pytest.raises(ValueError, match="substrate name"):
        read_dlkcat_output(model, p)


def test_unknown_substrate_listed_in_error(tmp_path):
    model = _ec_model(["r1"], ["alpha"])
    p = _write_dlkcat_output(tmp_path, [
        ("r1", "g1", "ghost", "CC", "M", "5.0"),
        ("r1", "g1", "phantom", "CC", "M", "5.0"),
    ])
    with pytest.raises(ValueError) as exc_info:
        read_dlkcat_output(model, p)
    msg = str(exc_info.value)
    assert "ghost" in msg
    assert "phantom" in msg


# --------------------------------------------------------------------------- #
# Reaction ID validation
# --------------------------------------------------------------------------- #

def test_unknown_reaction_raises(tmp_path):
    model = _ec_model(["r1"], ["alpha"])
    p = _write_dlkcat_output(tmp_path, [
        ("r_nonexistent", "g1", "alpha", "CC", "M", "5.0"),
    ])
    with pytest.raises(ValueError, match="reaction ID"):
        read_dlkcat_output(model, p)


def test_reaction_match_is_case_sensitive(tmp_path):
    """Reaction IDs are case-sensitive (matching cobrapy convention)."""
    model = _ec_model(["r1"], ["alpha"])
    p = _write_dlkcat_output(tmp_path, [
        ("R1", "g1", "alpha", "CC", "M", "5.0"),
    ])
    with pytest.raises(ValueError, match="reaction ID"):
        read_dlkcat_output(model, p)


# --------------------------------------------------------------------------- #
# Validation runs on the full file, not just the kept-rows subset
# --------------------------------------------------------------------------- #

def test_unknown_substrate_in_na_row_still_raises(tmp_path):
    """A bad substrate in a row that would be filtered out anyway still
    raises, because the validation runs against the full file."""
    model = _ec_model(["r1"], ["alpha"])
    p = _write_dlkcat_output(tmp_path, [
        ("r1", "g1", "alpha", "CC", "M", "5.0"),
        ("r1", "g1", "unknown_thing", "CN", "M", "NA"),
    ])
    with pytest.raises(ValueError, match="substrate name"):
        read_dlkcat_output(model, p)


# --------------------------------------------------------------------------- #
# Integration with apply_kcat_list
# --------------------------------------------------------------------------- #

def test_output_can_feed_apply_kcat_list(tmp_path):
    """End-to-end: read_dlkcat_output -> apply_kcat_list flow."""
    from geckopy.gather_kcats import apply_kcat_list

    model = _ec_model(["r1", "r2"], ["alpha", "beta"])
    p = _write_dlkcat_output(tmp_path, [
        ("r1", "g1", "alpha", "CC", "M", "5.0"),
        ("r2", "g2", "beta", "CN", "M", "7.0"),
    ])
    df = read_dlkcat_output(model, p)
    updated = apply_kcat_list(model, df)
    assert sorted(updated) == ["r1", "r2"]
    # apply_kcat_list lowercases the source token (no fuzzy metadata).
    assert model.ec.source == ["dlkcat", "dlkcat"]
    np.testing.assert_array_equal(model.ec.kcat, [5.0, 7.0])
