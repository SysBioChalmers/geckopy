"""Tests for get_ec_from_gem."""
import logging

import cobra
import numpy as np
import pytest
from scipy import sparse

from geckopy.ec_model import EcModel
from geckopy.ec_model.ec_data import EcData
from geckopy.get_enzyme_data import get_ec_from_gem
from geckopy.get_enzyme_data.ec_from_gem import _normalize_annotation


# --------------------------------------------------------------------------- #
# Test fixture builders
# --------------------------------------------------------------------------- #

def _make_ec_model(
    reactions: list[tuple[str, str | list[str] | None]],
    *,
    gecko_light: bool = False,
    ec_rxn_prefix: str = "",
) -> EcModel:
    """Build an EcModel with reactions carrying optional `ec-code` annotations.

    Each entry of ``reactions`` is ``(rxn_id, ec_value)`` where ``ec_value``
    is the value to set as ``reaction.annotation["ec-code"]`` (or ``None``
    to leave the annotation unset). One generic metabolite "A" is added
    so the reactions are well-formed.

    The returned EcModel has ``ec.rxns`` initialised in the same order
    as ``reactions``, with each ``ec.rxns[i]`` set to
    ``ec_rxn_prefix + reactions[i][0]`` so that gecko-light's 4-char
    prefix can be exercised.
    """
    model = EcModel("test", gecko_light=gecko_light)
    met = cobra.Metabolite("A", compartment="c")
    model.add_metabolites([met])

    for rxn_id, ec_value in reactions:
        rxn = cobra.Reaction(rxn_id)
        rxn.lower_bound = 0.0
        rxn.upper_bound = 1000.0
        rxn.add_metabolites({met: 1.0})
        if ec_value is not None:
            rxn.annotation["ec-code"] = ec_value
        model.add_reactions([rxn])

    n = len(reactions)
    model.ec = EcData(
        gecko_light=gecko_light,
        rxns=[ec_rxn_prefix + r for r, _ in reactions],
        kcat=np.full(n, np.nan, dtype=float),
        source=[""] * n,
        notes=[""] * n,
        eccodes=[""] * n,
        rxn_enz_mat=sparse.csr_matrix((n, 0), dtype=float),
    )
    return model


# --------------------------------------------------------------------------- #
# Trivial cases
# --------------------------------------------------------------------------- #

def test_empty_ec_rxns_is_a_noop():
    model = _make_ec_model([])
    get_ec_from_gem(model)
    assert model.ec.eccodes == []


def test_single_valid_ec_populated():
    model = _make_ec_model([("r1", "1.2.3.4")])
    get_ec_from_gem(model)
    assert model.ec.eccodes == ["1.2.3.4"]


def test_multiple_reactions_each_populated():
    model = _make_ec_model([
        ("r1", "1.2.3.4"),
        ("r2", "5.6.7.8"),
    ])
    get_ec_from_gem(model)
    assert model.ec.eccodes == ["1.2.3.4", "5.6.7.8"]


def test_missing_annotation_leaves_empty():
    model = _make_ec_model([("r1", None)])
    get_ec_from_gem(model)
    assert model.ec.eccodes == [""]


def test_empty_string_annotation_leaves_empty():
    model = _make_ec_model([("r1", "")])
    get_ec_from_gem(model)
    assert model.ec.eccodes == [""]


# --------------------------------------------------------------------------- #
# Annotation shape tolerance (cobrapy lets ec-code be str or list)
# --------------------------------------------------------------------------- #

def test_annotation_as_list_is_joined_with_semicolons():
    model = _make_ec_model([("r1", ["1.2.3.4", "5.6.7.8"])])
    get_ec_from_gem(model)
    assert model.ec.eccodes == ["1.2.3.4;5.6.7.8"]


def test_annotation_as_single_element_list():
    model = _make_ec_model([("r1", ["1.2.3.4"])])
    get_ec_from_gem(model)
    assert model.ec.eccodes == ["1.2.3.4"]


def test_annotation_as_list_with_empty_strings_filtered():
    model = _make_ec_model([("r1", ["1.2.3.4", "", "5.6.7.8"])])
    get_ec_from_gem(model)
    assert model.ec.eccodes == ["1.2.3.4;5.6.7.8"]


def test_annotation_string_with_semicolons_preserved():
    model = _make_ec_model([("r1", "1.2.3.4;5.6.7.8")])
    get_ec_from_gem(model)
    assert model.ec.eccodes == ["1.2.3.4;5.6.7.8"]


# --------------------------------------------------------------------------- #
# Wildcard placeholders
# --------------------------------------------------------------------------- #

def test_dash_in_last_level_accepted():
    model = _make_ec_model([("r1", "1.2.3.-")])
    get_ec_from_gem(model)
    assert model.ec.eccodes == ["1.2.3.-"]


def test_full_dash_cascade_accepted():
    model = _make_ec_model([("r1", "1.-.-.-")])
    get_ec_from_gem(model)
    assert model.ec.eccodes == ["1.-.-.-"]


def test_multi_ec_with_dashes():
    model = _make_ec_model([("r1", "1.2.3.-;5.6.-.-")])
    get_ec_from_gem(model)
    assert model.ec.eccodes == ["1.2.3.-;5.6.-.-"]


def test_multi_digit_levels():
    model = _make_ec_model([("r1", "3.4.21.1;2.7.11.1")])
    get_ec_from_gem(model)
    assert model.ec.eccodes == ["3.4.21.1;2.7.11.1"]


# --------------------------------------------------------------------------- #
# Validation: invalid strings cleared and warned
# --------------------------------------------------------------------------- #

def test_three_level_token_rejected(caplog):
    model = _make_ec_model([("r1", "1.2.3")])
    with caplog.at_level(logging.WARNING):
        get_ec_from_gem(model)
    assert model.ec.eccodes == [""]
    assert "skipped" in caplog.text
    assert "'1.2.3'" in caplog.text


def test_underscore_separator_rejected(caplog):
    model = _make_ec_model([("r1", "1_2_3_4")])
    with caplog.at_level(logging.WARNING):
        get_ec_from_gem(model)
    assert model.ec.eccodes == [""]
    assert "skipped" in caplog.text


def test_pipe_separator_for_multi_ec_rejected(caplog):
    """Per the MATLAB docstring: only `;` is a valid multi-EC separator."""
    model = _make_ec_model([("r1", "1.2.3.4|5.6.7.8")])
    with caplog.at_level(logging.WARNING):
        get_ec_from_gem(model)
    assert model.ec.eccodes == [""]
    assert "skipped" in caplog.text


def test_partial_invalid_in_multi_ec_rejects_whole_string(caplog):
    """If any token in the `;`-joined string is malformed, the whole entry is rejected."""
    model = _make_ec_model([("r1", "1.2.3.4;junk")])
    with caplog.at_level(logging.WARNING):
        get_ec_from_gem(model)
    assert model.ec.eccodes == [""]
    assert "skipped" in caplog.text


def test_mixed_valid_and_invalid_across_reactions(caplog):
    model = _make_ec_model([
        ("r1", "1.2.3.4"),
        ("r2", "junk"),
        ("r3", "5.6.7.-"),
    ])
    with caplog.at_level(logging.WARNING):
        get_ec_from_gem(model)
    assert model.ec.eccodes == ["1.2.3.4", "", "5.6.7.-"]
    assert "skipped 1" in caplog.text
    assert "'junk'" in caplog.text


def test_warning_lists_all_invalid_strings(caplog):
    model = _make_ec_model([
        ("r1", "junk1"),
        ("r2", "1.2.3.4"),
        ("r3", "junk2"),
    ])
    with caplog.at_level(logging.WARNING):
        get_ec_from_gem(model)
    assert "junk1" in caplog.text
    assert "junk2" in caplog.text


def test_no_warning_when_all_valid(caplog):
    model = _make_ec_model([
        ("r1", "1.2.3.4"),
        ("r2", "5.6.7.-"),
    ])
    with caplog.at_level(logging.WARNING):
        get_ec_from_gem(model)
    assert "skipped" not in caplog.text


# --------------------------------------------------------------------------- #
# ec_rxns parameter (cobrapy-style list of IDs)
# --------------------------------------------------------------------------- #

def test_ec_rxns_subset_updates_only_specified():
    model = _make_ec_model([
        ("r1", "1.2.3.4"),
        ("r2", "5.6.7.8"),
        ("r3", "9.9.9.9"),
    ])
    # Pre-populate r2's eccode with a sentinel; it should be left alone.
    model.ec.eccodes = ["", "preexisting", ""]

    get_ec_from_gem(model, ec_rxns=["r1", "r3"])

    assert model.ec.eccodes == ["1.2.3.4", "preexisting", "9.9.9.9"]


def test_ec_rxns_empty_is_a_noop():
    model = _make_ec_model([("r1", "1.2.3.4")])
    model.ec.eccodes = ["preexisting"]
    get_ec_from_gem(model, ec_rxns=[])
    assert model.ec.eccodes == ["preexisting"]


def test_ec_rxns_unknown_id_raises():
    model = _make_ec_model([("r1", "1.2.3.4")])
    with pytest.raises(ValueError, match="not present in model.ec.rxns"):
        get_ec_from_gem(model, ec_rxns=["nonexistent"])


def test_ec_rxns_accepts_iterables_not_just_lists():
    model = _make_ec_model([
        ("r1", "1.2.3.4"),
        ("r2", "5.6.7.8"),
    ])
    get_ec_from_gem(model, ec_rxns=(rid for rid in ["r1"]))
    assert model.ec.eccodes == ["1.2.3.4", ""]


# --------------------------------------------------------------------------- #
# gecko_light: 4-char prefix on ec.rxns
# --------------------------------------------------------------------------- #

def test_gecko_light_strips_4char_prefix_to_find_reaction():
    model = _make_ec_model(
        [("r1", "1.2.3.4"), ("r2", "5.6.7.8")],
        gecko_light=True,
        ec_rxn_prefix="001_",
    )
    # ec.rxns now reads ["001_r1", "001_r2"], model.reactions has r1, r2.
    get_ec_from_gem(model)
    assert model.ec.eccodes == ["1.2.3.4", "5.6.7.8"]


def test_gecko_light_unknown_underlying_reaction_raises():
    model = _make_ec_model(
        [("r1", "1.2.3.4")],
        gecko_light=True,
        ec_rxn_prefix="001_",
    )
    # Corrupt ec.rxns to point at a nonexistent stripped name.
    model.ec.rxns = ["XXXX_does_not_exist"]
    with pytest.raises(KeyError):
        get_ec_from_gem(model)


# --------------------------------------------------------------------------- #
# Integration with expand_model: _EXP_N reactions inherit ec-code
# --------------------------------------------------------------------------- #

def test_inherits_ec_code_from_expanded_reactions():
    """End-to-end: a reaction with isozymes is expanded by expand_model
    and its ec-code annotation should propagate to the _EXP_N copies,
    which get_ec_from_gem then picks up."""
    from geckopy.ec_model.pipeline import expand_model
    from geckopy.ec_model.pipeline.populate_ec import (
        allocate_ec_for_catalyzed_reactions,
    )

    model = EcModel("test")
    met = cobra.Metabolite("A", compartment="c")
    model.add_metabolites([met])

    rxn = cobra.Reaction("r1")
    rxn.lower_bound = 0.0
    rxn.upper_bound = 1000.0
    rxn.add_metabolites({met: 1.0})
    rxn.gene_reaction_rule = "g1 or g2"
    rxn.annotation["ec-code"] = "1.2.3.4"
    model.add_reactions([rxn])

    expand_model(model)
    allocate_ec_for_catalyzed_reactions(model)
    get_ec_from_gem(model)

    # ec.rxns should be the two expanded reactions, both inheriting the
    # parent's EC code.
    assert sorted(model.ec.rxns) == ["r1_EXP_1", "r1_EXP_2"]
    # ec.eccodes is parallel to ec.rxns; both entries should be "1.2.3.4".
    eccode_by_rxn = dict(zip(model.ec.rxns, model.ec.eccodes))
    assert eccode_by_rxn == {"r1_EXP_1": "1.2.3.4", "r1_EXP_2": "1.2.3.4"}


# --------------------------------------------------------------------------- #
# _normalize_annotation helper (covered indirectly above; spot-check)
# --------------------------------------------------------------------------- #

def test_normalize_handles_none():
    assert _normalize_annotation(None) == ""


def test_normalize_handles_empty_string():
    assert _normalize_annotation("") == ""


def test_normalize_handles_string():
    assert _normalize_annotation("1.2.3.4") == "1.2.3.4"


def test_normalize_handles_list():
    assert _normalize_annotation(["1.2.3.4", "5.6.7.8"]) == "1.2.3.4;5.6.7.8"


def test_normalize_handles_tuple():
    assert _normalize_annotation(("1.2.3.4", "5.6.7.8")) == "1.2.3.4;5.6.7.8"


def test_normalize_filters_empty_list_entries():
    assert _normalize_annotation(["1.2.3.4", "", None]) == "1.2.3.4"
