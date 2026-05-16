"""Tests for EcData validation and construction."""
import numpy as np
import pytest
from scipy import sparse

from geckopy import EcData


def test_default_ec_data_is_empty_and_valid():
    ec = EcData()
    assert ec.n_rxns == 0
    assert ec.n_enzymes == 0
    ec.validate()


def test_empty_for_reactions_shape():
    ec = EcData.empty_for_reactions(n_rxns=3, n_enzymes=2)
    assert ec.n_rxns == 3
    assert ec.n_enzymes == 2
    assert ec.kcat.shape == (3,)
    assert ec.mw.shape == (2,)
    assert ec.rxn_enz_mat.shape == (3, 2)
    assert np.isnan(ec.kcat).all()
    assert np.isnan(ec.mw).all()
    assert np.isnan(ec.concs).all()
    assert ec.rxns == ["", "", ""]
    ec.validate()


def test_validate_rejects_mismatched_reaction_fields():
    ec = EcData(
        rxns=["r1", "r2"],
        kcat=np.array([1.0]),  # wrong length
        source=["", ""],
        notes=["", ""],
        eccodes=["", ""],
    )
    with pytest.raises(ValueError, match="ec.kcat"):
        ec.validate()


def test_validate_rejects_mismatched_enzyme_fields():
    ec = EcData(
        genes=["g1", "g2"],
        enzymes=["P1", "P2"],
        mw=np.array([30.0, 40.0, 50.0]),  # wrong length
        sequence=["", ""],
        concs=np.array([np.nan, np.nan]),
    )
    with pytest.raises(ValueError, match="ec.mw"):
        ec.validate()


def test_validate_rejects_mismatched_coupling_matrix():
    ec = EcData.empty_for_reactions(n_rxns=3, n_enzymes=2)
    ec.rxn_enz_mat = sparse.csr_matrix((3, 5), dtype=float)  # wrong shape
    with pytest.raises(ValueError, match="rxn_enz_mat"):
        ec.validate()


def test_gecko_light_flag_roundtrips():
    ec = EcData(gecko_light=True)
    assert ec.gecko_light is True
