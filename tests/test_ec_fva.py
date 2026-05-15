"""Tests for ec_fva."""
import cobra
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from geckopy.ec_model import EcModel
from geckopy.ec_model.ec_data import EcData
from geckopy.utilities import ec_fva


# --------------------------------------------------------------------------- #
# Fixture builders
# --------------------------------------------------------------------------- #

def _build_simple_ec(*, ub_R: float = 5.0) -> tuple[EcModel, cobra.Model]:
    """ec_model with a single enzyme-constrained R; matching conv model
    has just R (no _REV, no _EXP)."""
    ec = EcModel("ec")
    A_e = cobra.Metabolite("A_e", compartment="e")
    A_c = cobra.Metabolite("A_c", compartment="c")
    B_c = cobra.Metabolite("B_c", compartment="c")
    pool = cobra.Metabolite("prot_pool", compartment="c")
    prot_E = cobra.Metabolite("prot_E", compartment="c")
    ec.add_metabolites([A_e, A_c, B_c, pool, prot_E])

    EX_A = cobra.Reaction("EX_A")
    EX_A.add_metabolites({A_e: -1.0})
    EX_A.lower_bound = -1000.0; EX_A.upper_bound = 0.0

    TR_A = cobra.Reaction("TR_A")
    TR_A.add_metabolites({A_e: -1.0, A_c: 1.0})
    TR_A.lower_bound = 0.0; TR_A.upper_bound = 1000.0

    R = cobra.Reaction("R")
    R.add_metabolites({A_c: -1.0, B_c: 1.0, prot_E: -1.0})
    R.lower_bound = 0.0; R.upper_bound = ub_R

    SK_B = cobra.Reaction("SK_B")
    SK_B.add_metabolites({B_c: -1.0})
    SK_B.lower_bound = 0.0; SK_B.upper_bound = 1000.0

    pool_ex = cobra.Reaction("prot_pool_exchange")
    pool_ex.add_metabolites({pool: 1.0})
    pool_ex.lower_bound = 0.0; pool_ex.upper_bound = 1000.0

    usage = cobra.Reaction("usage_prot_E")
    usage.add_metabolites({pool: -1.0, prot_E: 1.0})
    usage.lower_bound = 0.0; usage.upper_bound = 5.0

    ec.add_reactions([EX_A, TR_A, R, SK_B, pool_ex, usage])
    ec.objective = "SK_B"
    ec.ec = EcData(
        rxns=[], kcat=np.empty(0), source=[], notes=[], eccodes=[],
        rxn_enz_mat=sparse.csr_matrix((0, 0)),
    )

    # Conv model matches ec but without prot machinery and usage rxn.
    conv = cobra.Model("conv")
    A_e2 = cobra.Metabolite("A_e", compartment="e")
    A_c2 = cobra.Metabolite("A_c", compartment="c")
    B_c2 = cobra.Metabolite("B_c", compartment="c")
    conv.add_metabolites([A_e2, A_c2, B_c2])
    EX_A2 = cobra.Reaction("EX_A")
    EX_A2.add_metabolites({A_e2: -1.0})
    EX_A2.lower_bound = -1000.0; EX_A2.upper_bound = 0.0
    TR_A2 = cobra.Reaction("TR_A")
    TR_A2.add_metabolites({A_e2: -1.0, A_c2: 1.0})
    TR_A2.lower_bound = 0.0; TR_A2.upper_bound = 1000.0
    R2 = cobra.Reaction("R")
    R2.add_metabolites({A_c2: -1.0, B_c2: 1.0})
    R2.lower_bound = 0.0; R2.upper_bound = ub_R
    SK_B2 = cobra.Reaction("SK_B")
    SK_B2.add_metabolites({B_c2: -1.0})
    SK_B2.lower_bound = 0.0; SK_B2.upper_bound = 1000.0
    conv.add_reactions([EX_A2, TR_A2, R2, SK_B2])
    return ec, conv


def _build_with_rev() -> tuple[EcModel, cobra.Model]:
    """ec model with R and R_REV; conv model has just R."""
    ec, conv = _build_simple_ec()

    A_c = ec.metabolites.get_by_id("A_c")
    B_c = ec.metabolites.get_by_id("B_c")
    prot_E = ec.metabolites.get_by_id("prot_E")
    R_REV = cobra.Reaction("R_REV")
    R_REV.add_metabolites({B_c: -1.0, A_c: 1.0, prot_E: -1.0})
    R_REV.lower_bound = 0.0; R_REV.upper_bound = 5.0
    ec.add_reactions([R_REV])
    return ec, conv


# --------------------------------------------------------------------------- #
# Trivial cases
# --------------------------------------------------------------------------- #

def test_returns_dataframe_with_correct_index():
    ec, conv = _build_simple_ec()
    fva = ec_fva(ec, conv)
    assert isinstance(fva, pd.DataFrame)
    assert list(fva.index) == [r.id for r in conv.reactions]
    assert list(fva.columns) == ["min_flux", "max_flux"]


def test_empty_ec_model_returns_empty_df():
    ec = EcModel("empty")
    ec.ec = EcData(
        rxns=[], kcat=np.empty(0), source=[], notes=[], eccodes=[],
        rxn_enz_mat=sparse.csr_matrix((0, 0)),
    )
    conv = cobra.Model("empty_conv")
    fva = ec_fva(ec, conv)
    assert fva.empty
    assert list(fva.columns) == ["min_flux", "max_flux"]


# --------------------------------------------------------------------------- #
# Single-rxn FVA
# --------------------------------------------------------------------------- #

def test_min_le_max_for_every_reaction():
    """The basic FVA invariant."""
    ec, conv = _build_simple_ec()
    fva = ec_fva(ec, conv)
    assert (fva["min_flux"] <= fva["max_flux"] + 1e-9).all()


def test_R_max_capped_by_enzyme_ub():
    """R has enzyme stoich -1, and usage_prot_E has ub=5; max R = 5."""
    ec, conv = _build_simple_ec(ub_R=10.0)
    fva = ec_fva(ec, conv)
    assert fva.loc["R", "max_flux"] == pytest.approx(5.0, rel=1e-6)


def test_irreversible_rxn_min_is_zero():
    """R is forward-only (lb=0); min_flux should be 0."""
    ec, conv = _build_simple_ec()
    fva = ec_fva(ec, conv)
    assert fva.loc["R", "min_flux"] == pytest.approx(0.0, abs=1e-9)


def test_R_with_REV_yields_valid_fva():
    """When the ec model has both R and R_REV, FVA should still
    produce valid (min <= max) results in the conv space. (This
    topology doesn't allow R_REV to carry flux because B has no
    source other than R, but the function should still complete
    cleanly.)"""
    ec, conv = _build_with_rev()
    fva = ec_fva(ec, conv)
    assert (fva["min_flux"] <= fva["max_flux"] + 1e-9).all()
    # R is in the result; mapped flux is finite.
    assert "R" in fva.index
    assert np.isfinite(fva.loc["R", "min_flux"])
    assert np.isfinite(fva.loc["R", "max_flux"])


# --------------------------------------------------------------------------- #
# Output shape
# --------------------------------------------------------------------------- #

def test_index_is_named_rxn_id():
    ec, conv = _build_simple_ec()
    fva = ec_fva(ec, conv)
    assert fva.index.name == "rxn_id"


def test_includes_every_conv_reaction():
    ec, conv = _build_simple_ec()
    fva = ec_fva(ec, conv)
    assert set(fva.index) == {r.id for r in conv.reactions}


# --------------------------------------------------------------------------- #
# Bounds are respected
# --------------------------------------------------------------------------- #

def test_max_flux_within_rxn_upper_bound():
    ec, conv = _build_simple_ec()
    fva = ec_fva(ec, conv)
    for rxn_id in fva.index:
        rxn = conv.reactions.get_by_id(rxn_id)
        assert fva.loc[rxn_id, "max_flux"] <= rxn.upper_bound + 1e-9


def test_min_flux_within_rxn_lower_bound():
    ec, conv = _build_simple_ec()
    fva = ec_fva(ec, conv)
    for rxn_id in fva.index:
        rxn = conv.reactions.get_by_id(rxn_id)
        assert fva.loc[rxn_id, "min_flux"] >= rxn.lower_bound - 1e-9
