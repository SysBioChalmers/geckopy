"""Tests for set_kcat_for_reactions."""
from pathlib import Path

import cobra
import numpy as np
import pytest

from geckopy import EcModel, ModelAdapter, make_ec_model
from geckopy.ec_model.pipeline import (
    apply_kcat_constraints,
    set_kcat_for_reactions,
)

EXAMPLE_DIR = Path(__file__).parents[1] / "examples" / "ecTestGEM"


# --------------------------------------------------------------------------- #
# Fixtures
# --------------------------------------------------------------------------- #

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


def _kcat_at(model: EcModel, rxn_id: str) -> float:
    return float(model.ec.kcat[model.ec.rxns.index(rxn_id)])


def _source_at(model: EcModel, rxn_id: str) -> str:
    return model.ec.source[model.ec.rxns.index(rxn_id)]


def _get_s_coef(rxn: cobra.Reaction, met_id: str) -> float:
    for m, c in rxn.metabolites.items():
        if m.id == met_id:
            return c
    return 0.0


# --------------------------------------------------------------------------- #
# Single suffixed ID
# --------------------------------------------------------------------------- #

def test_set_single_suffixed_rxn():
    ec_model = _ectestgem_ec_model()
    updated = set_kcat_for_reactions(ec_model, ["R3"], 10.0, apply=False)
    assert updated == ["R3"]
    assert _kcat_at(ec_model, "R3") == 10.0
    assert _source_at(ec_model, "R3") == "manual"


def test_set_single_kcat_scalar_only_form():
    ec_model = _ectestgem_ec_model()
    updated = set_kcat_for_reactions(ec_model, ["R3"], 25.0, apply=False)
    assert updated == ["R3"]
    assert _kcat_at(ec_model, "R3") == 25.0


# --------------------------------------------------------------------------- #
# Un-suffixed ID expands to all isozymes
# --------------------------------------------------------------------------- #

def test_unsuffixed_rxn_expands_to_all_isozymes():
    """R2 expands to R2_EXP_1 and R2_EXP_2 (and R2_REV_EXP_1 / R2_REV_EXP_2
    if those count). Let's check what the actual matches are first."""
    ec_model = _ectestgem_ec_model()

    # Check what ec.rxns contains starting with R2.
    r2_entries = [r for r in ec_model.ec.rxns if r.startswith("R2")]
    # Expected: R2_EXP_1, R2_EXP_2, R2_REV_EXP_1, R2_REV_EXP_2.
    # All have the base name R2 after stripping _EXP_<n>.
    assert sorted(r2_entries) == sorted([
        "R2_EXP_1", "R2_EXP_2", "R2_REV_EXP_1", "R2_REV_EXP_2"
    ])

    # The base name "R2" should match all four, because stripping
    # _EXP_<n> from R2_EXP_1 gives "R2" but stripping from R2_REV_EXP_1
    # gives "R2_REV". So actually only R2_EXP_1 and R2_EXP_2 match "R2".
    updated = set_kcat_for_reactions(ec_model, ["R2"], 5.0, apply=False)
    assert sorted(updated) == ["R2_EXP_1", "R2_EXP_2"]
    assert _kcat_at(ec_model, "R2_EXP_1") == 5.0
    assert _kcat_at(ec_model, "R2_EXP_2") == 5.0
    # R2_REV_EXP_1 should NOT have been updated.
    assert np.isnan(_kcat_at(ec_model, "R2_REV_EXP_1"))


def test_unsuffixed_rxn_with_scalar_kcat_broadcasts():
    ec_model = _ectestgem_ec_model()
    set_kcat_for_reactions(ec_model, ["R2"], 7.5, apply=False)
    assert _kcat_at(ec_model, "R2_EXP_1") == 7.5
    assert _kcat_at(ec_model, "R2_EXP_2") == 7.5


def test_unsuffixed_rxn_with_list_kcat_forbidden():
    """Strict rule: un-suffixed ID expanding to multiple matches with a
    length-N kcat list is forbidden."""
    ec_model = _ectestgem_ec_model()
    with pytest.raises(ValueError, match="must be a scalar"):
        set_kcat_for_reactions(ec_model, ["R2"], [1.0, 2.0], apply=False)


def test_set_specific_isozyme_with_suffix():
    """Passing R2_EXP_1 should match only that one, not its sibling."""
    ec_model = _ectestgem_ec_model()
    updated = set_kcat_for_reactions(
        ec_model, ["R2_EXP_1"], 99.0, apply=False
    )
    assert updated == ["R2_EXP_1"]
    assert _kcat_at(ec_model, "R2_EXP_1") == 99.0
    assert np.isnan(_kcat_at(ec_model, "R2_EXP_2"))


# --------------------------------------------------------------------------- #
# Multiple input IDs
# --------------------------------------------------------------------------- #

def test_multiple_ids_with_scalar_kcat():
    """Scalar kcat broadcasts to all matched reactions across all inputs."""
    ec_model = _ectestgem_ec_model()
    updated = set_kcat_for_reactions(
        ec_model, ["R3", "R5"], 12.0, apply=False
    )
    assert sorted(updated) == ["R3", "R5"]
    assert _kcat_at(ec_model, "R3") == 12.0
    assert _kcat_at(ec_model, "R5") == 12.0


def test_multiple_ids_with_per_input_kcat_list():
    """A list whose length equals len(rxn_ids) gives one value per input,
    broadcasting if that input expanded to multiple matches."""
    ec_model = _ectestgem_ec_model()
    set_kcat_for_reactions(
        ec_model,
        ["R3", "R5"],
        [10.0, 20.0],
        apply=False,
    )
    assert _kcat_at(ec_model, "R3") == 10.0
    assert _kcat_at(ec_model, "R5") == 20.0


def test_multiple_ids_with_per_match_kcat_list_when_all_unique():
    """When every input matches exactly one reaction, a length-N kcat
    where N == total_matches is also accepted."""
    ec_model = _ectestgem_ec_model()
    # Both R3 and R5 are 1:1, so per-input == per-match.
    set_kcat_for_reactions(
        ec_model,
        ["R3", "R5"],
        [10.0, 20.0],
        apply=False,
    )
    assert _kcat_at(ec_model, "R3") == 10.0
    assert _kcat_at(ec_model, "R5") == 20.0


def test_kcat_list_wrong_length_raises():
    ec_model = _ectestgem_ec_model()
    with pytest.raises(ValueError, match="kcat has length 3"):
        set_kcat_for_reactions(
            ec_model, ["R3", "R5"], [1.0, 2.0, 3.0], apply=False,
        )


# --------------------------------------------------------------------------- #
# Unknown ID, empty input
# --------------------------------------------------------------------------- #

def test_unknown_rxn_id_raises():
    ec_model = _ectestgem_ec_model()
    with pytest.raises(ValueError, match="matched no entries"):
        set_kcat_for_reactions(ec_model, ["nonsense"], 1.0)


def test_unknown_with_known_still_raises():
    """If any one ID is unknown, the whole call must fail."""
    ec_model = _ectestgem_ec_model()
    with pytest.raises(ValueError, match="matched no entries"):
        set_kcat_for_reactions(ec_model, ["R3", "nonsense"], 1.0)


def test_empty_rxn_ids_returns_empty():
    ec_model = _ectestgem_ec_model()
    assert set_kcat_for_reactions(ec_model, [], 1.0) == []


# --------------------------------------------------------------------------- #
# apply=True propagates to S matrix
# --------------------------------------------------------------------------- #

def test_apply_true_writes_to_s_matrix():
    ec_model = _ectestgem_ec_model()
    set_kcat_for_reactions(ec_model, ["R3"], 10.0, apply=True)

    r3 = ec_model.reactions.get_by_id("R3")
    coef = _get_s_coef(r3, "prot_P4")
    expected = -(40000.0 / (10.0 * 3600.0))
    assert coef == pytest.approx(expected)


def test_apply_false_leaves_s_matrix_unchanged():
    ec_model = _ectestgem_ec_model()
    set_kcat_for_reactions(ec_model, ["R3"], 10.0, apply=False)

    r3 = ec_model.reactions.get_by_id("R3")
    assert _get_s_coef(r3, "prot_P4") == 0.0


def test_apply_default_is_true():
    """The default value of apply should be True."""
    ec_model = _ectestgem_ec_model()
    set_kcat_for_reactions(ec_model, ["R3"], 10.0)  # no apply arg
    r3 = ec_model.reactions.get_by_id("R3")
    assert _get_s_coef(r3, "prot_P4") != 0.0


def test_apply_only_updates_passed_reactions():
    """apply_kcat_constraints is called with update_rxns=updated_ids,
    not all reactions, so unrelated reactions retain whatever they had."""
    ec_model = _ectestgem_ec_model()
    # Pre-set R5 and apply manually.
    ec_model.ec.kcat[ec_model.ec.rxns.index("R5")] = 50.0
    apply_kcat_constraints(ec_model)
    r5_pre = _get_s_coef(ec_model.reactions.get_by_id("R5"), "prot_P5")
    assert r5_pre != 0.0

    # Now use set_kcat_for_reactions to update R3. R5 must be untouched.
    set_kcat_for_reactions(ec_model, ["R3"], 10.0)
    r5_post = _get_s_coef(ec_model.reactions.get_by_id("R5"), "prot_P5")
    assert r5_post == r5_pre


# --------------------------------------------------------------------------- #
# Source string
# --------------------------------------------------------------------------- #

def test_source_string_is_manual():
    ec_model = _ectestgem_ec_model()
    set_kcat_for_reactions(ec_model, ["R3"], 10.0, apply=False)
    assert _source_at(ec_model, "R3") == "manual"
