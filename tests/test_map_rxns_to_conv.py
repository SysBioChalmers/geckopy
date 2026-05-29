"""Tests for map_rxns_to_conv."""
import cobra
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from geckopy.ec_model import EcModel
from geckopy.ec_model.ec_data import EcData
from geckopy.utilities import MapRxnsResult, map_rxns_to_conv


# --------------------------------------------------------------------------- #
# Fixture builders
# --------------------------------------------------------------------------- #

def _conv_model(rxn_ids: list[str]) -> cobra.Model:
    """Conventional model with named reactions, no metabolites needed
    for mapping purposes."""
    model = cobra.Model("conv")
    rxns = [cobra.Reaction(rid) for rid in rxn_ids]
    model.add_reactions(rxns)
    return model


def _ec_model(rxn_ids: list[str]) -> EcModel:
    """ec model with the given reaction IDs (no stoichiometry needed)."""
    model = EcModel("ec")
    rxns = [cobra.Reaction(rid) for rid in rxn_ids]
    model.add_reactions(rxns)
    model.ec = EcData(
        rxns=[], kcat=np.empty(0), source=[], notes=[], eccodes=[],
        rxn_enz_mat=sparse.csr_matrix((0, 0)),
    )
    return model


# --------------------------------------------------------------------------- #
# Pre-flight
# --------------------------------------------------------------------------- #

def test_empty_flux_vect_raises():
    ec = _ec_model(["R"])
    conv = _conv_model(["R"])
    with pytest.raises(ValueError, match="empty"):
        map_rxns_to_conv(ec, conv, np.array([]))


def test_wrong_length_1d_raises():
    ec = _ec_model(["R"])
    conv = _conv_model(["R"])
    with pytest.raises(ValueError, match="length"):
        map_rxns_to_conv(ec, conv, np.array([1.0, 2.0]))


def test_wrong_axis0_2d_raises():
    ec = _ec_model(["R"])
    conv = _conv_model(["R"])
    with pytest.raises(ValueError, match="axis-0 length"):
        map_rxns_to_conv(ec, conv, np.array([[1.0], [2.0]]))


def test_3d_input_raises():
    ec = _ec_model(["R"])
    conv = _conv_model(["R"])
    with pytest.raises(ValueError, match="must be 1-D or 2-D"):
        map_rxns_to_conv(ec, conv, np.zeros((2, 2, 2)))


def test_missing_conventional_rxn_raises():
    ec = _ec_model(["R1"])
    conv = _conv_model(["R1", "R_missing"])
    with pytest.raises(ValueError, match="not found"):
        map_rxns_to_conv(ec, conv, np.array([1.0]))


# --------------------------------------------------------------------------- #
# Trivial cases
# --------------------------------------------------------------------------- #

def test_one_to_one_mapping_passes_through():
    ec = _ec_model(["R1", "R2"])
    conv = _conv_model(["R1", "R2"])
    result = map_rxns_to_conv(ec, conv, np.array([5.0, 7.0]))
    assert isinstance(result, MapRxnsResult)
    np.testing.assert_array_equal(result.mapped_flux, [5.0, 7.0])
    assert result.usage_enz == []


def test_reorder_to_conv_order():
    """ec_model rxn order != conv order; output follows conv order."""
    ec = _ec_model(["R1", "R2"])
    conv = _conv_model(["R2", "R1"])
    result = map_rxns_to_conv(ec, conv, np.array([5.0, 7.0]))
    np.testing.assert_array_equal(result.mapped_flux, [7.0, 5.0])


# --------------------------------------------------------------------------- #
# _REV reactions are negated
# --------------------------------------------------------------------------- #

def test_reverse_rxn_flux_negated_and_summed():
    """R has flux 5 forward; R_REV has flux 3 reverse-direction. Net
    flux = 5 - 3 = 2."""
    ec = _ec_model(["R", "R_REV"])
    conv = _conv_model(["R"])
    result = map_rxns_to_conv(ec, conv, np.array([5.0, 3.0]))
    assert result.mapped_flux[0] == pytest.approx(2.0)


def test_only_reverse_rxn_carries_flux():
    """Only R_REV carries flux 5 -> mapped flux = -5."""
    ec = _ec_model(["R", "R_REV"])
    conv = _conv_model(["R"])
    result = map_rxns_to_conv(ec, conv, np.array([0.0, 5.0]))
    assert result.mapped_flux[0] == pytest.approx(-5.0)


# --------------------------------------------------------------------------- #
# _EXP_<N> isozyme split combination
# --------------------------------------------------------------------------- #

def test_isozyme_expansion_summed_to_base_rxn():
    ec = _ec_model(["R_EXP_1", "R_EXP_2"])
    conv = _conv_model(["R"])
    result = map_rxns_to_conv(ec, conv, np.array([3.0, 4.0]))
    assert result.mapped_flux[0] == pytest.approx(7.0)


def test_combined_rev_and_exp_handled():
    """All four variants of R combined: R, R_REV, R_EXP_1, R_REV_EXP_1."""
    ec = _ec_model([
        "R", "R_REV", "R_EXP_1", "R_REV_EXP_1",
    ])
    conv = _conv_model(["R"])
    # Forward 10 + 4 = 14; reverse 1 + 2 = 3; net = 11.
    result = map_rxns_to_conv(ec, conv, np.array([10.0, 1.0, 4.0, 2.0]))
    assert result.mapped_flux[0] == pytest.approx(11.0)


# --------------------------------------------------------------------------- #
# Enzyme-usage extraction
# --------------------------------------------------------------------------- #

def test_usage_prot_rxns_extracted():
    ec = _ec_model(["R", "usage_prot_E"])
    conv = _conv_model(["R"])
    result = map_rxns_to_conv(ec, conv, np.array([5.0, 3.0]))
    assert result.usage_enz == ["E"]
    np.testing.assert_array_equal(result.enz_usage_flux, [3.0])


def test_prot_pool_exchange_labeled_pool():
    ec = _ec_model(["R", "prot_pool_exchange"])
    conv = _conv_model(["R"])
    result = map_rxns_to_conv(ec, conv, np.array([5.0, 7.0]))
    assert result.usage_enz == ["pool"]
    np.testing.assert_array_equal(result.enz_usage_flux, [7.0])


def test_usage_rxns_excluded_from_mapped_flux():
    """usage rxns appear only in enz_usage_flux, NOT in mapped_flux."""
    ec = _ec_model(["R", "usage_prot_E", "prot_pool_exchange"])
    conv = _conv_model(["R"])
    result = map_rxns_to_conv(ec, conv, np.array([5.0, 3.0, 7.0]))
    assert len(result.mapped_flux) == 1
    assert result.mapped_flux[0] == pytest.approx(5.0)
    assert sorted(result.usage_enz) == ["E", "pool"]


# --------------------------------------------------------------------------- #
# 1-D vs 2-D input
# --------------------------------------------------------------------------- #

def test_1d_input_yields_1d_output():
    ec = _ec_model(["R"])
    conv = _conv_model(["R"])
    result = map_rxns_to_conv(ec, conv, np.array([5.0]))
    assert result.mapped_flux.ndim == 1


def test_2d_input_yields_2d_output():
    ec = _ec_model(["R", "R_REV"])
    conv = _conv_model(["R"])
    flux = np.array([[5.0, 10.0], [1.0, 2.0]])
    result = map_rxns_to_conv(ec, conv, flux)
    assert result.mapped_flux.shape == (1, 2)
    np.testing.assert_array_equal(result.mapped_flux, [[4.0, 8.0]])


def test_pd_series_input_supported():
    ec = _ec_model(["R", "usage_prot_E"])
    conv = _conv_model(["R"])
    series = pd.Series({"R": 5.0, "usage_prot_E": 3.0})
    result = map_rxns_to_conv(ec, conv, series)
    assert result.mapped_flux[0] == pytest.approx(5.0)
    assert result.enz_usage_flux[0] == pytest.approx(3.0)


def test_pd_series_missing_rxn_defaults_to_zero():
    ec = _ec_model(["R1", "R2"])
    conv = _conv_model(["R1", "R2"])
    series = pd.Series({"R1": 5.0})  # R2 missing
    result = map_rxns_to_conv(ec, conv, series)
    np.testing.assert_array_equal(result.mapped_flux, [5.0, 0.0])
