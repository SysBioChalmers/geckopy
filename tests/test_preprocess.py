"""Tests for the preprocessing stages of make_ec_model.

These tests build small synthetic cobra models so every branch of the
ported function is exercised, independent of whether the ecTestGEM
example happens to trigger that branch.
"""
from pathlib import Path

import cobra

from geckopy import ModelAdapter
from geckopy.ec_model.pipeline import (
    convert_to_irreversible,
    invert_backwards_only_reactions,
    remove_pseudoreaction_gprs,
)


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #

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


def _build_model_with_bounds(
    reactions: list[tuple[str, dict[str, float], float, float]],
) -> cobra.Model:
    """Build from (rxn_id, {met_id: coef}, lb, ub) tuples."""
    model = cobra.Model("test")
    mets: dict[str, cobra.Metabolite] = {}
    for _, stoich, _, _ in reactions:
        for met_id in stoich:
            if met_id not in mets:
                mets[met_id] = cobra.Metabolite(met_id, compartment="c")

    for rxn_id, stoich, lb, ub in reactions:
        rxn = cobra.Reaction(rxn_id)
        rxn.lower_bound = lb
        rxn.upper_bound = ub
        rxn.add_metabolites({mets[m]: c for m, c in stoich.items()})
        model.add_reactions([rxn])
    return model


# --------------------------------------------------------------------------- #
# Stage 1: remove_pseudoreaction_gprs
# --------------------------------------------------------------------------- #

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


# --------------------------------------------------------------------------- #
# Stage 2: invert_backwards_only_reactions
# --------------------------------------------------------------------------- #

def test_inverts_single_backwards_only_reaction():
    model = _build_model_with_bounds([
        ("r1", {"A": -1.0, "B": 1.0}, -1000.0, 0.0),
    ])

    inverted = invert_backwards_only_reactions(model)

    r1 = model.reactions.get_by_id("r1")
    assert inverted == ["r1"]
    assert r1.lower_bound == 0.0
    assert r1.upper_bound == 1000.0
    # Stoichiometry was flipped.
    stoich = {m.id: c for m, c in r1.metabolites.items()}
    assert stoich == {"A": 1.0, "B": -1.0}


def test_does_not_touch_forward_only_reactions():
    model = _build_model_with_bounds([
        ("r1", {"A": -1.0, "B": 1.0}, 0.0, 1000.0),
    ])

    inverted = invert_backwards_only_reactions(model)

    r1 = model.reactions.get_by_id("r1")
    assert inverted == []
    assert r1.lower_bound == 0.0
    assert r1.upper_bound == 1000.0
    stoich = {m.id: c for m, c in r1.metabolites.items()}
    assert stoich == {"A": -1.0, "B": 1.0}


def test_does_not_touch_reversible_reactions():
    model = _build_model_with_bounds([
        ("r1", {"A": -1.0, "B": 1.0}, -1000.0, 1000.0),
    ])

    inverted = invert_backwards_only_reactions(model)
    assert inverted == []
    r1 = model.reactions.get_by_id("r1")
    assert r1.lower_bound == -1000.0
    assert r1.upper_bound == 1000.0


def test_does_not_touch_blocked_reactions():
    """A reaction with lb == 0 and ub == 0 should not be inverted."""
    model = _build_model_with_bounds([
        ("r1", {"A": -1.0, "B": 1.0}, 0.0, 0.0),
    ])

    inverted = invert_backwards_only_reactions(model)
    assert inverted == []


def test_does_not_touch_negative_only_reactions_with_nonzero_ub():
    """lb < 0 and ub < 0 (both negative) should not be inverted -- only
    lb < 0 with ub == 0 exactly triggers inversion."""
    model = _build_model_with_bounds([
        ("r1", {"A": -1.0, "B": 1.0}, -1000.0, -100.0),
    ])

    inverted = invert_backwards_only_reactions(model)
    assert inverted == []
    r1 = model.reactions.get_by_id("r1")
    assert r1.lower_bound == -1000.0
    assert r1.upper_bound == -100.0


def test_inverts_multiple_reactions_independently():
    model = _build_model_with_bounds([
        ("r1", {"A": -1.0, "B": 1.0}, -500.0, 0.0),
        ("r2", {"B": -2.0, "C": 3.0}, 0.0, 1000.0),     # not touched
        ("r3", {"C": -1.0, "D": 1.0}, -200.0, 0.0),
    ])

    inverted = invert_backwards_only_reactions(model)

    assert inverted == ["r1", "r3"]

    r1 = model.reactions.get_by_id("r1")
    assert r1.bounds == (0.0, 500.0)
    assert {m.id: c for m, c in r1.metabolites.items()} == {"A": 1.0, "B": -1.0}

    r2 = model.reactions.get_by_id("r2")
    assert r2.bounds == (0.0, 1000.0)
    assert {m.id: c for m, c in r2.metabolites.items()} == {"B": -2.0, "C": 3.0}

    r3 = model.reactions.get_by_id("r3")
    assert r3.bounds == (0.0, 200.0)
    assert {m.id: c for m, c in r3.metabolites.items()} == {"C": 1.0, "D": -1.0}


def test_inversion_preserves_gpr_and_name():
    model = _build_model_with_bounds([
        ("r1", {"A": -1.0, "B": 1.0}, -500.0, 0.0),
    ])
    r1 = model.reactions.get_by_id("r1")
    r1.name = "backwards reaction"
    r1.gene_reaction_rule = "g1 and g2"

    invert_backwards_only_reactions(model)

    r1 = model.reactions.get_by_id("r1")
    assert r1.name == "backwards reaction"
    assert r1.gene_reaction_rule == "g1 and g2"
    assert {g.id for g in r1.genes} == {"g1", "g2"}


def test_returns_empty_list_when_nothing_to_invert():
    model = _build_model_with_bounds([
        ("r1", {"A": -1.0, "B": 1.0}, 0.0, 1000.0),
        ("r2", {"B": -1.0, "C": 1.0}, -1000.0, 1000.0),
    ])
    assert invert_backwards_only_reactions(model) == []


# --------------------------------------------------------------------------- #
# Stage 4: convert_to_irreversible
#
# Stage 3 (set `model.rev` based on bounds) is absorbed into this stage,
# since in MATLAB the `rev` field is only ever consumed by convertToIrrev
# and never read again. Exchange reactions are excluded here to match
# MATLAB behavior.
# --------------------------------------------------------------------------- #

def test_splits_single_reversible_non_exchange():
    model = _build_model_with_bounds([
        ("r1", {"A": -1.0, "B": 1.0}, -500.0, 1000.0),
    ])

    added = convert_to_irreversible(model)
    assert added == ["r1_REV"]

    fwd = model.reactions.get_by_id("r1")
    rev = model.reactions.get_by_id("r1_REV")

    assert fwd.bounds == (0.0, 1000.0)
    assert {m.id: c for m, c in fwd.metabolites.items()} == {"A": -1.0, "B": 1.0}

    assert rev.bounds == (0.0, 500.0)
    assert {m.id: c for m, c in rev.metabolites.items()} == {"A": 1.0, "B": -1.0}


def test_does_not_split_forward_only_reaction():
    model = _build_model_with_bounds([
        ("r1", {"A": -1.0, "B": 1.0}, 0.0, 1000.0),
    ])
    added = convert_to_irreversible(model)
    assert added == []
    assert "r1_REV" not in {r.id for r in model.reactions}


def test_does_not_split_exchange_reaction_even_if_reversible():
    """Exchange reactions (one metabolite) are excluded from the
    irreversibility step, regardless of bounds."""
    model = _build_model_with_bounds([
        ("EX_A", {"A": -1.0}, -1000.0, 1000.0),
    ])
    added = convert_to_irreversible(model)
    assert added == []
    ex = model.reactions.get_by_id("EX_A")
    assert ex.bounds == (-1000.0, 1000.0)


def test_splits_multiple_mixed_reactions():
    model = _build_model_with_bounds([
        ("r1", {"A": -1.0, "B": 1.0}, -500.0, 1000.0),   # split
        ("r2", {"B": -2.0, "C": 3.0}, 0.0, 1000.0),      # forward only
        ("EX_A", {"A": -1.0}, -1000.0, 1000.0),          # exchange
        ("r3", {"C": -1.0, "D": 1.0}, -200.0, 200.0),    # split
    ])

    added = convert_to_irreversible(model)
    assert added == ["r1_REV", "r3_REV"]

    assert model.reactions.get_by_id("r1").bounds == (0.0, 1000.0)
    assert model.reactions.get_by_id("r1_REV").bounds == (0.0, 500.0)
    assert model.reactions.get_by_id("r2").bounds == (0.0, 1000.0)
    assert model.reactions.get_by_id("EX_A").bounds == (-1000.0, 1000.0)
    assert model.reactions.get_by_id("r3").bounds == (0.0, 200.0)
    assert model.reactions.get_by_id("r3_REV").bounds == (0.0, 200.0)


def test_reverse_reaction_inherits_gpr():
    model = _build_model_with_bounds([
        ("r1", {"A": -1.0, "B": 1.0}, -500.0, 1000.0),
    ])
    model.reactions.get_by_id("r1").gene_reaction_rule = "g1 and g2"

    convert_to_irreversible(model)

    rev = model.reactions.get_by_id("r1_REV")
    assert rev.gene_reaction_rule == "g1 and g2"
    assert {g.id for g in rev.genes} == {"g1", "g2"}


def test_forward_reaction_lb_is_clamped_to_zero():
    """After splitting, the original (forward) reaction's lower bound
    is clamped to 0."""
    model = _build_model_with_bounds([
        ("r1", {"A": -1.0, "B": 1.0}, -500.0, 1000.0),
    ])
    convert_to_irreversible(model)
    assert model.reactions.get_by_id("r1").lower_bound == 0.0


def test_no_reverse_reaction_has_negative_bound():
    """Check the invariant that after conversion, no reaction in the
    model is allowed to carry negative flux (except exchanges)."""
    model = _build_model_with_bounds([
        ("r1", {"A": -1.0, "B": 1.0}, -500.0, 1000.0),
        ("r2", {"B": -1.0, "C": 1.0}, -1000.0, 0.0),      # blocked reverse
        ("EX_A", {"A": -1.0}, -1000.0, 1000.0),
    ])
    convert_to_irreversible(model)
    for rxn in model.reactions:
        if rxn.boundary:
            continue
        assert rxn.lower_bound >= 0, f"{rxn.id} still has lb < 0"


def test_returns_empty_list_when_nothing_to_split():
    model = _build_model_with_bounds([
        ("r1", {"A": -1.0, "B": 1.0}, 0.0, 1000.0),
        ("EX_A", {"A": -1.0}, -1000.0, 1000.0),
    ])
    assert convert_to_irreversible(model) == []


def test_conversion_is_idempotent_after_first_pass():
    """Running convert_to_irreversible twice should not create
    `_REV_REV` reactions, because the first pass already clamped
    all non-exchange lb to 0."""
    model = _build_model_with_bounds([
        ("r1", {"A": -1.0, "B": 1.0}, -500.0, 1000.0),
    ])
    convert_to_irreversible(model)
    second = convert_to_irreversible(model)
    assert second == []
    assert "r1_REV_REV" not in {r.id for r in model.reactions}
