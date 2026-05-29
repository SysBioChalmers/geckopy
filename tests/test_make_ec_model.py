"""Tests for the top-level make_ec_model orchestrator."""
import logging
from pathlib import Path

import cobra
import numpy as np
import pytest

from geckopy import EcModel, ModelAdapter, make_ec_model
from geckopy.databases import UniprotDB, load_uniprot_tsv

EXAMPLE_DIR = Path(__file__).parents[1] / "examples" / "ecTestGEM"


# --------------------------------------------------------------------------- #
# End-to-end on the ecTestGEM fixture
# --------------------------------------------------------------------------- #

def _load_fresh_ectestgem() -> tuple[cobra.Model, ModelAdapter]:
    adapter = ModelAdapter.from_folder(EXAMPLE_DIR)
    model = cobra.io.read_sbml_model(str(adapter.params.conv_gem))
    return model, adapter


def test_make_ec_model_returns_ec_model_instance():
    model, adapter = _load_fresh_ectestgem()
    ec_model = make_ec_model(model, adapter)
    assert isinstance(ec_model, EcModel)


def test_make_ec_model_no_unmatched_genes_in_ectestgem(caplog):
    """All 5 genes in ecTestGEM match UniProt, so no warning is logged."""
    model, adapter = _load_fresh_ectestgem()
    with caplog.at_level(logging.WARNING, logger="geckopy.ec_model.make_ec_model"):
        make_ec_model(model, adapter)
    assert "not found in UniProt" not in caplog.text


def test_make_ec_model_populates_ec_rxns():
    model, adapter = _load_fresh_ectestgem()
    ec_model = make_ec_model(model, adapter)
    assert sorted(ec_model.ec.rxns) == sorted([
        "R2_EXP_1", "R2_EXP_2",
        "R2_REV_EXP_1", "R2_REV_EXP_2",
        "R3", "R5",
    ])


def test_make_ec_model_populates_ec_genes_alphabetically():
    model, adapter = _load_fresh_ectestgem()
    ec_model = make_ec_model(model, adapter)
    assert ec_model.ec.genes == ["G1", "G2", "G3", "G4", "G5"]
    assert ec_model.ec.enzymes == ["P1", "P2", "P3", "P4", "P5"]
    np.testing.assert_array_almost_equal(
        ec_model.ec.mw, [10000.0, 20000.0, 30000.0, 40000.0, 50000.0]
    )


def test_make_ec_model_kcat_initialized_nan():
    model, adapter = _load_fresh_ectestgem()
    ec_model = make_ec_model(model, adapter)
    assert np.isnan(ec_model.ec.kcat).all()


def test_make_ec_model_adds_prot_metabolites():
    model, adapter = _load_fresh_ectestgem()
    ec_model = make_ec_model(model, adapter)
    met_ids = {m.id for m in ec_model.metabolites}
    assert "prot_pool" in met_ids
    for p in ["prot_P1", "prot_P2", "prot_P3", "prot_P4", "prot_P5"]:
        assert p in met_ids


def test_make_ec_model_adds_usage_and_exchange_reactions():
    model, adapter = _load_fresh_ectestgem()
    ec_model = make_ec_model(model, adapter)
    rxn_ids = {r.id for r in ec_model.reactions}
    assert "prot_pool_exchange" in rxn_ids
    for u in [f"usage_prot_P{i}" for i in range(1, 6)]:
        assert u in rxn_ids


def test_make_ec_model_coupling_matrix_shape_and_values():
    model, adapter = _load_fresh_ectestgem()
    ec_model = make_ec_model(model, adapter)
    mat = ec_model.ec.rxn_enz_mat.toarray()
    assert mat.shape == (6, 5)
    assert set(np.unique(mat)).issubset({0.0, 1.0})
    assert int(mat.sum()) == 8


def test_make_ec_model_validates_ec_consistency():
    model, adapter = _load_fresh_ectestgem()
    ec_model = make_ec_model(model, adapter)
    ec_model.validate_ec()


def test_make_ec_model_attaches_adapter():
    model, adapter = _load_fresh_ectestgem()
    ec_model = make_ec_model(model, adapter)
    assert ec_model.adapter is adapter


# --------------------------------------------------------------------------- #
# Branching and error handling
# --------------------------------------------------------------------------- #

def test_make_ec_model_accepts_preloaded_uniprot_db():
    model, adapter = _load_fresh_ectestgem()
    db = load_uniprot_tsv(EXAMPLE_DIR / "data" / "uniprot.tsv")
    ec_model = make_ec_model(model, adapter, uniprot_db=db)
    assert ec_model.ec.n_enzymes == 5


def test_make_ec_model_gecko_light_raises():
    model, adapter = _load_fresh_ectestgem()
    with pytest.raises(NotImplementedError, match="gecko_light"):
        make_ec_model(model, adapter, gecko_light=True)


def test_make_ec_model_missing_uniprot_tsv_raises(tmp_path):
    """If the adapter points at a folder without uniprot.tsv, and no
    explicit uniprot_db is given, raise FileNotFoundError."""
    toml = (
        'conv_gem = "dummy.xml"\n'
        'org_name = "test"\n'
        'enzyme_comp = "c"\n'
    )
    (tmp_path / "model_adapter.toml").write_text(toml)
    adapter = ModelAdapter.from_folder(tmp_path)

    model = cobra.Model("m")
    m = cobra.Metabolite("A", compartment="c")
    model.add_metabolites([m])
    r = cobra.Reaction("r1")
    r.add_metabolites({m: -1})
    r.gene_reaction_rule = "g1"
    model.add_reactions([r])

    with pytest.raises(FileNotFoundError, match="uniprot.tsv"):
        make_ec_model(model, adapter)


def test_make_ec_model_refuses_rerun_on_populated_ecmodel():
    """Calling make_ec_model twice on the same model should raise the
    second time, because the model already has a populated ec."""
    model, adapter = _load_fresh_ectestgem()
    ec_model = make_ec_model(model, adapter)

    with pytest.raises(ValueError, match="already has a populated ec"):
        make_ec_model(ec_model, adapter)


def test_make_ec_model_unmatched_genes_logged_and_annotated(tmp_path, caplog):
    """A model whose gene is absent from the UniProt DB should produce
    a logged warning and per-reaction note, but no return value."""
    toml = (
        'conv_gem = "dummy.xml"\n'
        'org_name = "test"\n'
        'enzyme_comp = "c"\n'
    )
    (tmp_path / "model_adapter.toml").write_text(toml)
    adapter = ModelAdapter.from_folder(tmp_path)

    model = cobra.Model("m")
    m_a = cobra.Metabolite("A", compartment="c")
    m_b = cobra.Metabolite("B", compartment="c")
    model.add_metabolites([m_a, m_b])
    r = cobra.Reaction("r1")
    r.add_metabolites({m_a: -1, m_b: 1})
    r.lower_bound, r.upper_bound = 0.0, 1000.0
    r.gene_reaction_rule = "g_known and g_missing"
    model.add_reactions([r])

    db = UniprotDB(
        ids=["P1"],
        genes=["g_known"],
        eccodes=[""],
        mw=np.array([10.0]),
        sequences=["S"],
    )

    with caplog.at_level(logging.WARNING, logger="geckopy.ec_model.make_ec_model"):
        ec_model = make_ec_model(model, adapter, uniprot_db=db)

    assert "g_missing" in caplog.text
    warned = ec_model.reactions.get_by_id("r1")
    assert "geckopy_warning" in warned.notes
    assert "g_missing" in warned.notes["geckopy_warning"]
