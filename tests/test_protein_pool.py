"""Tests for stages 9, 10, 11, 12 of make_ec_model (protein pool machinery)."""
from pathlib import Path

import cobra
import numpy as np
import pytest

from geckopy import EcModel, ModelAdapter
from geckopy.databases import UniprotDB, load_uniprot_tsv
from geckopy.ec_model.pipeline import (
    add_protein_pool_exchange_reaction,
    add_protein_pool_pseudometabolite,
    add_protein_pseudometabolites,
    add_protein_usage_reactions,
    allocate_ec_for_catalyzed_reactions,
    build_rxn_enzyme_coupling,
    convert_to_irreversible,
    expand_model,
    invert_backwards_only_reactions,
    populate_enzyme_data,
    remove_pseudoreaction_gprs,
)

EXAMPLE_DIR = Path(__file__).parents[1] / "examples" / "ecTestGEM"


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #

def _minimal_adapter(tmp_path: Path, enzyme_comp: str = "c") -> ModelAdapter:
    (tmp_path / "model_adapter.toml").write_text(
        f'conv_gem = "dummy.xml"\norg_name = "test"\nenzyme_comp = "{enzyme_comp}"\n'
    )
    return ModelAdapter.from_folder(tmp_path)


def _tiny_ec_model_with_two_enzymes(tmp_path: Path) -> EcModel:
    """Build an EcModel manually with enough state for stages 9-12."""
    adapter = _minimal_adapter(tmp_path)
    base = cobra.Model("tiny")
    m_a = cobra.Metabolite("A", compartment="c")
    m_b = cobra.Metabolite("B", compartment="c")
    base.add_metabolites([m_a, m_b])

    r1 = cobra.Reaction("r1")
    r1.add_metabolites({m_a: -1, m_b: 1})
    r1.lower_bound, r1.upper_bound = 0.0, 1000.0
    r1.gene_reaction_rule = "g1 and g2"
    base.add_reactions([r1])

    ec_model = EcModel.from_cobra(base, adapter=adapter)
    allocate_ec_for_catalyzed_reactions(ec_model)

    db = UniprotDB(
        ids=["P1", "P2"],
        genes=["g1", "g2"],
        eccodes=["", ""],
        mw=np.array([10.0, 20.0]),
        sequences=["S1", "S2"],
    )
    populate_enzyme_data(ec_model, db)
    build_rxn_enzyme_coupling(ec_model)
    return ec_model


# --------------------------------------------------------------------------- #
# Compartment resolution
# --------------------------------------------------------------------------- #

def test_enzyme_comp_matches_by_name_when_names_differ_from_ids(tmp_path):
    """Realistic yeast-like case: compartment name "cytoplasm", ID "c"."""
    adapter = _minimal_adapter(tmp_path, enzyme_comp="cytoplasm")
    base = cobra.Model("test")
    base.compartments = {"c": "cytoplasm", "e": "extracellular"}

    m = cobra.Metabolite("A", compartment="c")
    base.add_metabolites([m])
    r = cobra.Reaction("r1")
    r.add_metabolites({m: -1})
    r.gene_reaction_rule = "g1"
    base.add_reactions([r])

    ec_model = EcModel.from_cobra(base, adapter=adapter)
    allocate_ec_for_catalyzed_reactions(ec_model)
    db = UniprotDB(
        ids=["P1"], genes=["g1"], eccodes=[""],
        mw=np.array([10.0]), sequences=["S"],
    )
    populate_enzyme_data(ec_model, db)

    added = add_protein_pseudometabolites(ec_model)
    assert added == ["prot_P1"]
    assert ec_model.metabolites.get_by_id("prot_P1").compartment == "c"


def test_enzyme_comp_not_found_raises(tmp_path):
    adapter = _minimal_adapter(tmp_path, enzyme_comp="nucleus")
    ec_model = _tiny_ec_model_with_two_enzymes(tmp_path)
    ec_model.adapter = adapter  # override compartment target

    with pytest.raises(ValueError, match="enzyme compartment|enzyme_comp"):
        add_protein_pseudometabolites(ec_model)


# --------------------------------------------------------------------------- #
# Stage 9: add_protein_pseudometabolites
# --------------------------------------------------------------------------- #

def test_stage9_adds_one_met_per_enzyme(tmp_path):
    ec_model = _tiny_ec_model_with_two_enzymes(tmp_path)
    added = add_protein_pseudometabolites(ec_model)

    assert added == ["prot_P1", "prot_P2"]
    assert "prot_P1" in {m.id for m in ec_model.metabolites}
    assert "prot_P2" in {m.id for m in ec_model.metabolites}


def test_stage9_sets_compartment(tmp_path):
    ec_model = _tiny_ec_model_with_two_enzymes(tmp_path)
    add_protein_pseudometabolites(ec_model)
    assert ec_model.metabolites.get_by_id("prot_P1").compartment == "c"


def test_stage9_sets_sbo_annotation(tmp_path):
    ec_model = _tiny_ec_model_with_two_enzymes(tmp_path)
    add_protein_pseudometabolites(ec_model)
    met = ec_model.metabolites.get_by_id("prot_P1")
    assert met.annotation.get("sbo") == "SBO:0000252"


def test_stage9_is_idempotent(tmp_path):
    ec_model = _tiny_ec_model_with_two_enzymes(tmp_path)
    add_protein_pseudometabolites(ec_model)
    second = add_protein_pseudometabolites(ec_model)
    assert second == []
    met_ids = {m.id for m in ec_model.metabolites}
    assert sum(1 for i in met_ids if i.startswith("prot_")) == 2


def test_stage9_empty_enzymes_adds_nothing(tmp_path):
    """Stage 9 should no-op when ec.enzymes is empty, without needing
    to resolve a compartment (so it works even on an empty model)."""
    adapter = _minimal_adapter(tmp_path)
    base = cobra.Model("empty")
    ec_model = EcModel.from_cobra(base, adapter=adapter)
    added = add_protein_pseudometabolites(ec_model)
    assert added == []
    assert "prot_" not in {m.id for m in ec_model.metabolites}


def test_stages_9_through_12_on_empty_ec_but_nonempty_model(tmp_path):
    """Stages 10 and 12 still run (pool is always added) even when there
    are no enzymes to add pseudometabolites for, as long as the model
    has a valid enzyme_comp."""
    adapter = _minimal_adapter(tmp_path)
    base = cobra.Model("model_without_genes")
    m = cobra.Metabolite("A", compartment="c")
    base.add_metabolites([m])
    # cobrapy populates model.compartments from metabolites, so "c" exists now.

    ec_model = EcModel.from_cobra(base, adapter=adapter)
    # ec.enzymes is empty (no stage-7 call).

    assert add_protein_pseudometabolites(ec_model) == []
    add_protein_pool_pseudometabolite(ec_model)
    assert "prot_pool" in {mm.id for mm in ec_model.metabolites}
    assert add_protein_usage_reactions(ec_model) == []  # no enzymes, no usages
    add_protein_pool_exchange_reaction(ec_model)
    assert "prot_pool_exchange" in {r.id for r in ec_model.reactions}


def test_stage9_requires_adapter():
    ec_model = EcModel("test")
    ec_model.adapter = None
    with pytest.raises(ValueError, match="adapter"):
        add_protein_pseudometabolites(ec_model)


# --------------------------------------------------------------------------- #
# Stage 10: add_protein_pool_pseudometabolite
# --------------------------------------------------------------------------- #

def test_stage10_adds_prot_pool(tmp_path):
    ec_model = _tiny_ec_model_with_two_enzymes(tmp_path)
    add_protein_pool_pseudometabolite(ec_model)
    assert "prot_pool" in {m.id for m in ec_model.metabolites}
    assert ec_model.metabolites.get_by_id("prot_pool").compartment == "c"


def test_stage10_is_idempotent(tmp_path):
    ec_model = _tiny_ec_model_with_two_enzymes(tmp_path)
    add_protein_pool_pseudometabolite(ec_model)
    add_protein_pool_pseudometabolite(ec_model)
    assert sum(
        1 for m in ec_model.metabolites if m.id == "prot_pool"
    ) == 1


# --------------------------------------------------------------------------- #
# Stage 11: add_protein_usage_reactions
# --------------------------------------------------------------------------- #

def test_stage11_adds_one_usage_reaction_per_enzyme(tmp_path):
    ec_model = _tiny_ec_model_with_two_enzymes(tmp_path)
    add_protein_pseudometabolites(ec_model)
    add_protein_pool_pseudometabolite(ec_model)

    added = add_protein_usage_reactions(ec_model)
    assert added == ["usage_prot_P1", "usage_prot_P2"]


def test_stage11_stoichiometry_forward_direction(tmp_path):
    ec_model = _tiny_ec_model_with_two_enzymes(tmp_path)
    add_protein_pseudometabolites(ec_model)
    add_protein_pool_pseudometabolite(ec_model)
    add_protein_usage_reactions(ec_model)

    r = ec_model.reactions.get_by_id("usage_prot_P1")
    stoich = {m.id: c for m, c in r.metabolites.items()}
    assert stoich == {"prot_pool": -1.0, "prot_P1": 1.0}
    assert r.lower_bound == 0.0
    assert r.upper_bound == 1000.0


def test_stage11_sets_gpr_to_single_gene(tmp_path):
    ec_model = _tiny_ec_model_with_two_enzymes(tmp_path)
    add_protein_pseudometabolites(ec_model)
    add_protein_pool_pseudometabolite(ec_model)
    add_protein_usage_reactions(ec_model)

    r1 = ec_model.reactions.get_by_id("usage_prot_P1")
    r2 = ec_model.reactions.get_by_id("usage_prot_P2")
    assert r1.gene_reaction_rule == "g1"
    assert r2.gene_reaction_rule == "g2"


def test_stage11_sets_subsystem(tmp_path):
    ec_model = _tiny_ec_model_with_two_enzymes(tmp_path)
    add_protein_pseudometabolites(ec_model)
    add_protein_pool_pseudometabolite(ec_model)
    add_protein_usage_reactions(ec_model)
    assert ec_model.reactions.get_by_id(
        "usage_prot_P1"
    ).subsystem == "Protein usage"


def test_stage11_is_idempotent(tmp_path):
    ec_model = _tiny_ec_model_with_two_enzymes(tmp_path)
    add_protein_pseudometabolites(ec_model)
    add_protein_pool_pseudometabolite(ec_model)
    add_protein_usage_reactions(ec_model)
    second = add_protein_usage_reactions(ec_model)
    assert second == []


# --------------------------------------------------------------------------- #
# Stage 12: add_protein_pool_exchange_reaction
# --------------------------------------------------------------------------- #

def test_stage12_adds_pool_exchange(tmp_path):
    ec_model = _tiny_ec_model_with_two_enzymes(tmp_path)
    add_protein_pool_pseudometabolite(ec_model)
    add_protein_pool_exchange_reaction(ec_model)

    r = ec_model.reactions.get_by_id("prot_pool_exchange")
    stoich = {m.id: c for m, c in r.metabolites.items()}
    assert stoich == {"prot_pool": 1.0}
    assert r.lower_bound == 0.0
    assert r.upper_bound == 1000.0
    assert r.subsystem == "Protein usage"


def test_stage12_is_idempotent(tmp_path):
    ec_model = _tiny_ec_model_with_two_enzymes(tmp_path)
    add_protein_pool_pseudometabolite(ec_model)
    add_protein_pool_exchange_reaction(ec_model)
    add_protein_pool_exchange_reaction(ec_model)
    count = sum(
        1 for r in ec_model.reactions if r.id == "prot_pool_exchange"
    )
    assert count == 1


# --------------------------------------------------------------------------- #
# End-to-end: stages 1-12 on ecTestGEM fixture
# --------------------------------------------------------------------------- #

def _run_all_stages_on_ectestgem() -> EcModel:
    adapter = ModelAdapter.from_folder(EXAMPLE_DIR)
    cobra_model = cobra.io.read_sbml_model(str(adapter.params.conv_gem))

    remove_pseudoreaction_gprs(cobra_model, adapter)
    invert_backwards_only_reactions(cobra_model)
    convert_to_irreversible(cobra_model)
    expand_model(cobra_model)

    ec_model = EcModel.from_cobra(cobra_model, adapter=adapter)
    allocate_ec_for_catalyzed_reactions(ec_model)
    uniprot_db = load_uniprot_tsv(EXAMPLE_DIR / "data" / "uniprot.tsv")
    populate_enzyme_data(ec_model, uniprot_db)
    build_rxn_enzyme_coupling(ec_model)

    add_protein_pseudometabolites(ec_model)
    add_protein_pool_pseudometabolite(ec_model)
    add_protein_usage_reactions(ec_model)
    add_protein_pool_exchange_reaction(ec_model)

    return ec_model


def test_ectestgem_has_expected_prot_metabolites():
    ec_model = _run_all_stages_on_ectestgem()
    met_ids = {m.id for m in ec_model.metabolites}
    # 5 prot_Pi + 1 prot_pool + 4 original metabolites = 10
    assert "prot_pool" in met_ids
    for pid in ["prot_P1", "prot_P2", "prot_P3", "prot_P4", "prot_P5"]:
        assert pid in met_ids


def test_ectestgem_has_expected_usage_reactions():
    ec_model = _run_all_stages_on_ectestgem()
    rxn_ids = {r.id for r in ec_model.reactions}
    for uid in [
        "usage_prot_P1", "usage_prot_P2", "usage_prot_P3",
        "usage_prot_P4", "usage_prot_P5",
    ]:
        assert uid in rxn_ids
    assert "prot_pool_exchange" in rxn_ids


def test_ectestgem_usage_reactions_have_correct_stoichiometry():
    ec_model = _run_all_stages_on_ectestgem()
    for enzyme in ["P1", "P2", "P3", "P4", "P5"]:
        r = ec_model.reactions.get_by_id(f"usage_prot_{enzyme}")
        stoich = {m.id: c for m, c in r.metabolites.items()}
        assert stoich == {"prot_pool": -1.0, f"prot_{enzyme}": 1.0}


def test_ectestgem_pool_exchange_stoichiometry():
    ec_model = _run_all_stages_on_ectestgem()
    r = ec_model.reactions.get_by_id("prot_pool_exchange")
    stoich = {m.id: c for m, c in r.metabolites.items()}
    assert stoich == {"prot_pool": 1.0}


def test_ectestgem_usage_gpr_matches_enzyme():
    ec_model = _run_all_stages_on_ectestgem()
    # P1<->G1, P2<->G2, ..., P5<->G5 based on the uniprot.tsv fixture.
    expected = {
        "usage_prot_P1": "G1",
        "usage_prot_P2": "G2",
        "usage_prot_P3": "G3",
        "usage_prot_P4": "G4",
        "usage_prot_P5": "G5",
    }
    for rxn_id, gpr in expected.items():
        assert ec_model.reactions.get_by_id(rxn_id).gene_reaction_rule == gpr


def test_ectestgem_validate_ec_after_all_stages():
    ec_model = _run_all_stages_on_ectestgem()
    ec_model.validate_ec()


# --------------------------------------------------------------------------- #
# set_prot_pool_size
# --------------------------------------------------------------------------- #

from geckopy.ec_model.pipeline import set_prot_pool_size


def _tiny_ec_model_with_pool(tmp_path: Path) -> EcModel:
    """Tiny EcModel with the full stages 1-12 completed, so pool_exchange exists."""
    ec_model = _tiny_ec_model_with_two_enzymes(tmp_path)
    add_protein_pseudometabolites(ec_model)
    add_protein_pool_pseudometabolite(ec_model)
    add_protein_usage_reactions(ec_model)
    add_protein_pool_exchange_reaction(ec_model)
    return ec_model


def test_set_prot_pool_size_default_uses_adapter_params(tmp_path):
    ec_model = _tiny_ec_model_with_pool(tmp_path)
    # adapter defaults: p_tot=0.5, f=0.5, sigma=0.5 -> 0.5*0.5*0.5*1000 = 125.0
    bound = set_prot_pool_size(ec_model)
    assert bound == 125.0
    assert ec_model.reactions.get_by_id(
        "prot_pool_exchange"
    ).upper_bound == 125.0


def test_set_prot_pool_size_explicit_overrides_adapter(tmp_path):
    ec_model = _tiny_ec_model_with_pool(tmp_path)
    bound = set_prot_pool_size(ec_model, p_tot=0.4, f=0.3, sigma=0.2)
    assert bound == pytest.approx(0.4 * 0.3 * 0.2 * 1000.0)
    assert ec_model.reactions.get_by_id(
        "prot_pool_exchange"
    ).upper_bound == pytest.approx(24.0)


def test_set_prot_pool_size_partial_override(tmp_path):
    """Override only p_tot; f and sigma come from adapter (both 0.5)."""
    ec_model = _tiny_ec_model_with_pool(tmp_path)
    bound = set_prot_pool_size(ec_model, p_tot=1.0)
    assert bound == 1.0 * 0.5 * 0.5 * 1000.0  # = 250.0


def test_set_prot_pool_size_raises_without_pool_exchange(tmp_path):
    """If prot_pool_exchange is absent, the function must raise."""
    ec_model = _tiny_ec_model_with_two_enzymes(tmp_path)
    # Intentionally do not add the pool machinery.
    with pytest.raises(ValueError, match="prot_pool_exchange"):
        set_prot_pool_size(ec_model)


def test_set_prot_pool_size_raises_without_adapter_when_args_missing(tmp_path):
    """If the adapter is absent and an arg is missing, raise."""
    ec_model = _tiny_ec_model_with_pool(tmp_path)
    ec_model.adapter = None
    with pytest.raises(ValueError, match="No ModelAdapter available"):
        set_prot_pool_size(ec_model)  # all three args missing


def test_set_prot_pool_size_works_without_adapter_when_all_args_given(tmp_path):
    """No adapter is fine if the user supplies all three values."""
    ec_model = _tiny_ec_model_with_pool(tmp_path)
    ec_model.adapter = None
    bound = set_prot_pool_size(ec_model, p_tot=0.5, f=0.5, sigma=0.5)
    assert bound == 125.0


def test_set_prot_pool_size_is_idempotent(tmp_path):
    """Calling twice with the same args yields the same bound."""
    ec_model = _tiny_ec_model_with_pool(tmp_path)
    bound1 = set_prot_pool_size(ec_model, p_tot=0.3, f=0.4, sigma=0.5)
    bound2 = set_prot_pool_size(ec_model, p_tot=0.3, f=0.4, sigma=0.5)
    assert bound1 == bound2
