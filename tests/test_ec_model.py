"""Tests for the EcModel class."""
import cobra
import pytest

from geckopy import EcData, EcModel


def _tiny_cobra_model() -> cobra.Model:
    model = cobra.Model("tiny")
    m_a = cobra.Metabolite("A", compartment="c")
    m_b = cobra.Metabolite("B", compartment="c")
    r1 = cobra.Reaction("r1")
    r1.add_metabolites({m_a: -1, m_b: 1})
    r1.lower_bound, r1.upper_bound = 0, 1000
    model.add_reactions([r1])
    return model


def test_ec_model_has_empty_ec_by_default():
    ec_model = EcModel("test")
    assert isinstance(ec_model.ec, EcData)
    assert ec_model.ec.n_rxns == 0
    assert ec_model.adapter is None
    assert ec_model.ec.gecko_light is False


def test_gecko_light_flag_propagates():
    ec_model = EcModel("test", gecko_light=True)
    assert ec_model.ec.gecko_light is True


def test_from_cobra_preserves_reactions():
    base = _tiny_cobra_model()
    ec_model = EcModel.from_cobra(base, adapter=None)
    assert "r1" in {r.id for r in ec_model.reactions}
    assert ec_model.adapter is None


def test_validate_ec_catches_unknown_reaction_id():
    base = _tiny_cobra_model()
    ec_model = EcModel.from_cobra(base, adapter=None)
    ec_model.ec = EcData.empty_for_reactions(n_rxns=1)
    ec_model.ec.rxns = ["not_in_model"]
    with pytest.raises(ValueError, match="not present in the model"):
        ec_model.validate_ec()


def test_validate_ec_accepts_consistent_state():
    base = _tiny_cobra_model()
    ec_model = EcModel.from_cobra(base, adapter=None)
    ec_model.ec = EcData.empty_for_reactions(n_rxns=1)
    ec_model.ec.rxns = ["r1"]
    ec_model.validate_ec()  # should not raise
