"""Tests for the preprocessing stages of make_ec_model.

These tests build small synthetic cobra models so every branch of the
ported function is exercised, independent of whether the ecTestGEM
example happens to trigger that branch.
"""
from pathlib import Path

import cobra

from geckopy import ModelAdapter
from geckopy.ec_model.pipeline import remove_pseudoreaction_gprs


def _minimal_adapter(tmp_path: Path) -> ModelAdapter:
    """Adapter rooted in tmp_path. No real SBML needed for these tests."""
    (tmp_path / "model_adapter.toml").write_text(
        'conv_gem = "dummy.xml"\norg_name = "test"\n'
    )
    return ModelAdapter.from_folder(tmp_path)


def _build_model(reactions: list[tuple[str, str, str]]) -> cobra.Model:
    """Build a cobra.Model from (rxn_id, rxn_name, gpr) tuples."""
    model = cobra.Model("test")
    m = cobra.Metabolite("M_dummy", compartment="c")
    for rxn_id, rxn_name, gpr in reactions:
        rxn = cobra.Reaction(rxn_id, name=rxn_name)
        rxn.add_metabolites({m: -1})
        rxn.gene_reaction_rule = gpr
        model.add_reactions([rxn])
    return model


def test_clears_gprs_for_reactions_with_pseudoreaction_in_name(tmp_path):
    adapter = _minimal_adapter(tmp_path)
    model = _build_model([
        ("r1", "real reaction", "g1"),
        ("r2", "biomass pseudoreaction", "g2"),
        ("r3", "maintenance pseudoreaction", "g3 and g4"),
    ])

    cleared = remove_pseudoreaction_gprs(model, adapter)

    assert cleared == ["r2", "r3"]
    assert model.reactions.get_by_id("r1").gene_reaction_rule == "g1"
    assert model.reactions.get_by_id("r2").gene_reaction_rule == ""
    assert model.reactions.get_by_id("r3").gene_reaction_rule == ""


def test_case_sensitive_match_like_matlab(tmp_path):
    adapter = _minimal_adapter(tmp_path)
    model = _build_model([
        ("r1", "Pseudoreaction biomass", "g1"),  # capital P: no match
        ("r2", "biomass pseudoreaction", "g2"),  # match
    ])

    cleared = remove_pseudoreaction_gprs(model, adapter)

    assert cleared == ["r2"]
    assert model.reactions.get_by_id("r1").gene_reaction_rule == "g1"


def test_clears_gprs_listed_in_pseudoRxns_tsv(tmp_path):
    adapter = _minimal_adapter(tmp_path)
    (tmp_path / "data").mkdir()
    (tmp_path / "data" / "pseudoRxns.tsv").write_text(
        "r1\tgrowth\nr3\tmaintenance\n"
    )
    model = _build_model([
        ("r1", "normal name", "g1"),
        ("r2", "normal name", "g2"),
        ("r3", "normal name", "g3"),
    ])

    cleared = remove_pseudoreaction_gprs(model, adapter)

    assert cleared == ["r1", "r3"]
    assert model.reactions.get_by_id("r2").gene_reaction_rule == "g2"


def test_combines_name_match_and_tsv(tmp_path):
    adapter = _minimal_adapter(tmp_path)
    (tmp_path / "data").mkdir()
    (tmp_path / "data" / "pseudoRxns.tsv").write_text("r1\n")
    model = _build_model([
        ("r1", "from tsv", "g1"),
        ("r2", "biomass pseudoreaction", "g2"),
        ("r3", "regular", "g3"),
    ])

    cleared = remove_pseudoreaction_gprs(model, adapter)

    assert cleared == ["r1", "r2"]
    assert model.reactions.get_by_id("r3").gene_reaction_rule == "g3"


def test_ignores_tsv_ids_not_present_in_model(tmp_path):
    adapter = _minimal_adapter(tmp_path)
    (tmp_path / "data").mkdir()
    (tmp_path / "data" / "pseudoRxns.tsv").write_text("nonexistent\nr1\n")
    model = _build_model([("r1", "regular", "g1")])

    cleared = remove_pseudoreaction_gprs(model, adapter)
    assert cleared == ["r1"]


def test_no_pseudoreactions_leaves_model_unchanged(tmp_path):
    adapter = _minimal_adapter(tmp_path)
    model = _build_model([
        ("r1", "regular", "g1"),
        ("r2", "another regular", "g2"),
    ])

    cleared = remove_pseudoreaction_gprs(model, adapter)
    assert cleared == []
    assert model.reactions.get_by_id("r1").gene_reaction_rule == "g1"
    assert model.reactions.get_by_id("r2").gene_reaction_rule == "g2"


def test_clearing_gpr_also_clears_reaction_genes(tmp_path):
    """cobrapy detail worth pinning: setting gene_reaction_rule=''
    also removes the reaction's gene associations. The ported function
    relies on this so callers do not need a separate gene cleanup step."""
    adapter = _minimal_adapter(tmp_path)
    model = _build_model([("r1", "a pseudoreaction", "g1 and g2")])

    assert len(model.reactions.get_by_id("r1").genes) == 2
    remove_pseudoreaction_gprs(model, adapter)
    assert len(model.reactions.get_by_id("r1").genes) == 0