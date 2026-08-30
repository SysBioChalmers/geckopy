"""Tests for get_reactions_from_enzyme."""
from pathlib import Path

import cobra
import pandas as pd
import pytest

from geckopy import EcModel, ModelAdapter, make_ec_model
from geckopy.ec_model.pipeline import (
    get_reactions_from_enzyme,
    set_kcat_for_reactions,
)

EXAMPLE_DIR = Path(__file__).parents[1] / "examples" / "ecTestGEM"


_ECTESTGEM_CACHE: EcModel | None = None


def _ectestgem_ec_model() -> EcModel:
    """Cached build of the ecTestGEM ecModel; deep-copied per call."""
    import copy as _copy
    global _ECTESTGEM_CACHE
    if _ECTESTGEM_CACHE is None:
        adapter = ModelAdapter.from_folder(EXAMPLE_DIR)
        cobra_model = cobra.io.read_sbml_model(str(adapter.params.conv_gem))
        _ECTESTGEM_CACHE = make_ec_model(cobra_model, adapter)
    return _copy.deepcopy(_ECTESTGEM_CACHE)


def test_returns_dataframe_with_expected_columns():
    ec_model = _ectestgem_ec_model()
    df = get_reactions_from_enzyme(ec_model, "P4")
    assert isinstance(df, pd.DataFrame)
    assert list(df.columns) == ["rxn_id", "kcat", "name", "gpr"]


def test_p4_catalyzes_only_r3():
    ec_model = _ectestgem_ec_model()
    df = get_reactions_from_enzyme(ec_model, "P4")
    assert df["rxn_id"].tolist() == ["R3"]
    assert df["gpr"].tolist() == ["G4"]


def test_p1_catalyzes_both_r2_expansions():
    ec_model = _ectestgem_ec_model()
    df = get_reactions_from_enzyme(ec_model, "P1")
    assert sorted(df["rxn_id"].tolist()) == ["R2_EXP_1", "R2_REV_EXP_1"]


def test_p3_catalyzes_both_r2_isozyme_branches():
    ec_model = _ectestgem_ec_model()
    df = get_reactions_from_enzyme(ec_model, "P3")
    assert sorted(df["rxn_id"].tolist()) == ["R2_EXP_2", "R2_REV_EXP_2"]


def test_kcat_initially_zero():
    ec_model = _ectestgem_ec_model()
    df = get_reactions_from_enzyme(ec_model, "P4")
    assert (df["kcat"] == 0).all()


def test_kcat_reflects_set_value():
    ec_model = _ectestgem_ec_model()
    set_kcat_for_reactions(ec_model, ["R3"], 10.0, apply=False)
    df = get_reactions_from_enzyme(ec_model, "P4")
    assert df["kcat"].tolist() == [10.0]


def test_unknown_protein_raises():
    ec_model = _ectestgem_ec_model()
    with pytest.raises(ValueError, match="not found in ec.enzymes"):
        get_reactions_from_enzyme(ec_model, "P_UNKNOWN")


def test_case_sensitive_match():
    ec_model = _ectestgem_ec_model()
    with pytest.raises(ValueError):
        get_reactions_from_enzyme(ec_model, "p4")


def test_rxn_names_populated_from_cobra_model():
    ec_model = _ectestgem_ec_model()
    df = get_reactions_from_enzyme(ec_model, "P4")
    assert df["name"].tolist() == ["R3"]


def test_dataframe_supports_natural_pandas_idioms():
    """The returned DataFrame composes with ordinary pandas operations
    (e.g. sort_values), not just positional indexing."""
    ec_model = _ectestgem_ec_model()
    set_kcat_for_reactions(ec_model, ["R2_EXP_1"], 50.0, apply=False)
    set_kcat_for_reactions(ec_model, ["R2_REV_EXP_1"], 25.0, apply=False)

    df = get_reactions_from_enzyme(ec_model, "P1")
    sorted_df = df.sort_values("kcat", ascending=False)
    assert sorted_df["rxn_id"].iloc[0] == "R2_EXP_1"
    assert sorted_df["rxn_id"].iloc[1] == "R2_REV_EXP_1"
