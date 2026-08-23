"""Tests for the EcModel class."""
import cobra
import numpy as np
import pytest
from scipy import sparse

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
    ec_model.ec = EcData.empty(n_rxns=1)
    ec_model.ec.rxns = ["not_in_model"]
    with pytest.raises(ValueError, match="not present in the model"):
        ec_model.validate_ec()


def test_validate_ec_accepts_consistent_state():
    base = _tiny_cobra_model()
    ec_model = EcModel.from_cobra(base, adapter=None)
    ec_model.ec = EcData.empty(n_rxns=1)
    ec_model.ec.rxns = ["r1"]
    ec_model.validate_ec()  # should not raise


# --------------------------------------------------------------------------- #
# copy(): the ec substructure must be cloned, not shared
#
# cobra.Model.copy carries non-network attributes over by reference, so
# without EcModel.copy every one of these would write through to the
# original. MATLAB GECKO has value semantics throughout.
# --------------------------------------------------------------------------- #

def _ec_model_with_data() -> EcModel:
    """A tiny EcModel with every ec field populated."""
    ec_model = EcModel.from_cobra(_tiny_cobra_model(), adapter=None)
    ec_model.ec = EcData.empty(n_rxns=1, n_enzymes=1)
    ec_model.ec.rxns = ["r1"]
    ec_model.ec.kcat[0] = 10.0
    ec_model.ec.source = ["brenda"]
    ec_model.ec.notes = ["a note"]
    ec_model.ec.eccodes = ["1.1.1.1"]
    ec_model.ec.genes = ["G1"]
    ec_model.ec.enzymes = ["P1"]
    ec_model.ec.mw[0] = 10000.0
    ec_model.ec.sequence = ["MRAL"]
    ec_model.ec.concs[0] = 0.5
    ec_model.ec.rxn_enz_mat = sparse.csr_matrix(np.array([[1.0]]))
    return ec_model


def test_copy_returns_an_ec_model():
    assert isinstance(_ec_model_with_data().copy(), EcModel)


def test_copy_gives_an_independent_ec_object():
    original = _ec_model_with_data()
    assert original.copy().ec is not original.ec


def test_copy_preserves_ec_contents():
    original = _ec_model_with_data()
    copied = original.copy()
    assert copied.ec.rxns == original.ec.rxns
    assert copied.ec.source == original.ec.source
    assert copied.ec.notes == original.ec.notes
    assert copied.ec.eccodes == original.ec.eccodes
    assert copied.ec.genes == original.ec.genes
    assert copied.ec.enzymes == original.ec.enzymes
    assert copied.ec.sequence == original.ec.sequence
    assert copied.ec.gecko_light == original.ec.gecko_light
    np.testing.assert_array_equal(copied.ec.kcat, original.ec.kcat)
    np.testing.assert_array_equal(copied.ec.mw, original.ec.mw)
    np.testing.assert_array_equal(copied.ec.concs, original.ec.concs)
    np.testing.assert_array_equal(
        copied.ec.rxn_enz_mat.toarray(), original.ec.rxn_enz_mat.toarray()
    )


@pytest.mark.parametrize("field_name", ["kcat", "mw", "concs"])
def test_mutating_a_numpy_field_on_the_copy_leaves_the_original(field_name):
    original = _ec_model_with_data()
    before = getattr(original.ec, field_name)[0]
    copied = original.copy()
    getattr(copied.ec, field_name)[0] = 99.0
    assert getattr(original.ec, field_name)[0] == before
    assert getattr(copied.ec, field_name)[0] == 99.0


@pytest.mark.parametrize(
    "field_name", ["rxns", "source", "notes", "eccodes", "genes", "enzymes", "sequence"]
)
def test_mutating_a_list_field_on_the_copy_leaves_the_original(field_name):
    original = _ec_model_with_data()
    before = list(getattr(original.ec, field_name))
    copied = original.copy()
    getattr(copied.ec, field_name)[0] = "CHANGED"
    assert getattr(original.ec, field_name) == before
    assert getattr(copied.ec, field_name)[0] == "CHANGED"


def test_appending_to_a_list_field_on_the_copy_leaves_the_original():
    original = _ec_model_with_data()
    copied = original.copy()
    copied.ec.rxns.append("r2")
    assert original.ec.rxns == ["r1"]


def test_mutating_rxn_enz_mat_on_the_copy_leaves_the_original():
    original = _ec_model_with_data()
    copied = original.copy()
    copied.ec.rxn_enz_mat[0, 0] = 42.0
    assert original.ec.rxn_enz_mat[0, 0] == 1.0


def test_replacing_ec_wholesale_on_the_copy_leaves_the_original():
    original = _ec_model_with_data()
    copied = original.copy()
    copied.ec = EcData.empty(n_rxns=0)
    assert original.ec.rxns == ["r1"]


def test_copy_rebinds_the_enzymes_view_to_the_copy():
    """The lazy EnzymeView holds the model it reads through; a copy whose
    view still pointed at the original would read and write the original's
    ec data."""
    original = _ec_model_with_data()
    copied = original.copy()
    copied.ec.enzymes = ["P2"]
    assert [e.id for e in copied.enzymes] == ["P2"]
    assert [e.id for e in original.enzymes] == ["P1"]


def test_copy_shares_the_adapter_reference():
    """The adapter is project configuration, not model state."""
    original = EcModel.from_cobra(_tiny_cobra_model(), adapter=None)
    sentinel = object()
    original.adapter = sentinel
    assert original.copy().adapter is sentinel


def test_copy_still_copies_the_network_like_cobra():
    original = _ec_model_with_data()
    copied = original.copy()
    assert {r.id for r in copied.reactions} == {r.id for r in original.reactions}
    assert copied.reactions.get_by_id("r1") is not original.reactions.get_by_id("r1")
    copied.reactions.get_by_id("r1").upper_bound = 7.0
    assert original.reactions.get_by_id("r1").upper_bound == 1000.0
