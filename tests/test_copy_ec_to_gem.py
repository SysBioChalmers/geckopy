"""Tests for copy_ec_to_gem."""
import cobra
import numpy as np
from scipy import sparse

from geckopy.ec_model import EcModel
from geckopy.ec_model.ec_data import EcData
from geckopy.get_enzyme_data import copy_ec_to_gem, fill_eccodes_from_gem


# --------------------------------------------------------------------------- #
# Fixture builders
# --------------------------------------------------------------------------- #

def _make_ec_model(
    rxn_specs: list[tuple[str, object]],
    *,
    gecko_light: bool = False,
    ec_rxn_prefix: str = "",
    ec_eccodes: list[str] | None = None,
) -> EcModel:
    """Build an EcModel with reactions carrying the given initial
    ``ec-code`` annotation values.

    ``rxn_specs`` is a list of ``(rxn_id, initial_ec_annotation)``;
    ``initial_ec_annotation`` may be a string, list, or None (no
    annotation). The model's ``ec.rxns`` is set to
    ``ec_rxn_prefix + rxn_id`` for each rxn so gecko-light's prefix
    handling can be exercised.
    """
    model = EcModel("test", gecko_light=gecko_light)
    met = cobra.Metabolite("A", compartment="c")
    model.add_metabolites([met])

    for rxn_id, initial in rxn_specs:
        rxn = cobra.Reaction(rxn_id)
        rxn.lower_bound = 0.0
        rxn.upper_bound = 1000.0
        rxn.add_metabolites({met: 1.0})
        if initial is not None:
            rxn.annotation["ec-code"] = initial
        model.add_reactions([rxn])

    n = len(rxn_specs)
    if ec_eccodes is None:
        ec_eccodes = [""] * n

    model.ec = EcData(
        gecko_light=gecko_light,
        rxns=[ec_rxn_prefix + r for r, _ in rxn_specs],
        kcat=np.full(n, np.nan, dtype=float),
        source=[""] * n,
        notes=[""] * n,
        eccodes=list(ec_eccodes),
        rxn_enz_mat=sparse.csr_matrix((n, 0), dtype=float),
    )
    return model


# --------------------------------------------------------------------------- #
# Trivial cases
# --------------------------------------------------------------------------- #

def test_empty_ec_rxns_is_a_noop():
    model = _make_ec_model([])
    copy_ec_to_gem(model)
    # Nothing to assert; just shouldn't crash.


def test_empty_eccodes_does_not_pollute_existing_annotation():
    """An empty ec.eccodes entry must not clobber a populated annotation,
    even with overwrite=True."""
    model = _make_ec_model(
        [("r1", "9.9.9.9")], ec_eccodes=[""],
    )
    copy_ec_to_gem(model, overwrite=True)
    rxn = model.reactions.get_by_id("r1")
    assert rxn.annotation["ec-code"] == "9.9.9.9"


def test_empty_eccodes_does_not_create_annotation():
    model = _make_ec_model([("r1", None)], ec_eccodes=[""])
    copy_ec_to_gem(model)
    rxn = model.reactions.get_by_id("r1")
    assert "ec-code" not in rxn.annotation


def test_single_value_written_as_list():
    """Output value is a list of EC tokens, even for a single token."""
    model = _make_ec_model([("r1", None)], ec_eccodes=["1.1.1.1"])
    copy_ec_to_gem(model)
    assert model.reactions.get_by_id("r1").annotation["ec-code"] == ["1.1.1.1"]


def test_multi_ec_split_on_semicolon():
    model = _make_ec_model(
        [("r1", None)], ec_eccodes=["1.1.1.1;2.2.2.2;3.3.3.3"],
    )
    copy_ec_to_gem(model)
    assert model.reactions.get_by_id("r1").annotation["ec-code"] == [
        "1.1.1.1", "2.2.2.2", "3.3.3.3",
    ]


# --------------------------------------------------------------------------- #
# overwrite semantics
# --------------------------------------------------------------------------- #

def test_overwrite_false_fills_only_empty_annotations():
    model = _make_ec_model(
        [
            ("r1", None),                     # missing -> filled
            ("r2", ""),                       # empty string -> filled
            ("r3", []),                       # empty list -> filled
            ("r4", "9.9.9.9"),                # populated -> preserved
            ("r5", ["8.8.8.8"]),              # populated list -> preserved
        ],
        ec_eccodes=["1.1.1.1", "1.1.1.1", "1.1.1.1", "1.1.1.1", "1.1.1.1"],
    )
    copy_ec_to_gem(model, overwrite=False)
    assert model.reactions.get_by_id("r1").annotation["ec-code"] == ["1.1.1.1"]
    assert model.reactions.get_by_id("r2").annotation["ec-code"] == ["1.1.1.1"]
    assert model.reactions.get_by_id("r3").annotation["ec-code"] == ["1.1.1.1"]
    assert model.reactions.get_by_id("r4").annotation["ec-code"] == "9.9.9.9"
    assert model.reactions.get_by_id("r5").annotation["ec-code"] == ["8.8.8.8"]


def test_overwrite_true_replaces_populated_annotations():
    model = _make_ec_model(
        [("r1", "9.9.9.9"), ("r2", ["8.8.8.8"])],
        ec_eccodes=["1.1.1.1", "2.2.2.2"],
    )
    copy_ec_to_gem(model, overwrite=True)
    assert model.reactions.get_by_id("r1").annotation["ec-code"] == ["1.1.1.1"]
    assert model.reactions.get_by_id("r2").annotation["ec-code"] == ["2.2.2.2"]


def test_overwrite_default_is_false():
    """Default is `overwrite=False`: calling without the kwarg leaves an
    existing annotation untouched."""
    model = _make_ec_model(
        [("r1", "9.9.9.9")], ec_eccodes=["1.1.1.1"],
    )
    copy_ec_to_gem(model)  # no overwrite kwarg
    assert model.reactions.get_by_id("r1").annotation["ec-code"] == "9.9.9.9"


# --------------------------------------------------------------------------- #
# Integrity edge cases
# --------------------------------------------------------------------------- #

def test_reaction_in_model_but_not_in_ec_rxns_left_untouched():
    """`model.ec.rxns` is a subset of `model.reactions`. Reactions
    outside that subset must not gain annotations."""
    model = _make_ec_model(
        [("r1", None), ("r2", None)],
        ec_eccodes=["1.1.1.1", "2.2.2.2"],
    )
    # Drop r2 from ec.rxns to simulate a reaction that's in the model
    # but not catalysed.
    model.ec.rxns = ["r1"]
    model.ec.eccodes = ["1.1.1.1"]
    copy_ec_to_gem(model)
    assert model.reactions.get_by_id("r1").annotation["ec-code"] == ["1.1.1.1"]
    assert "ec-code" not in model.reactions.get_by_id("r2").annotation


def test_ec_rxn_not_in_model_silently_skipped():
    """If `ec.rxns` references a reaction that's no longer in the
    cobra model, that entry is silently skipped (not an error)."""
    model = _make_ec_model([("r1", None)], ec_eccodes=["1.1.1.1"])
    # Append a phantom ec.rxn that doesn't match any cobra reaction.
    model.ec.rxns = ["r1", "phantom"]
    model.ec.eccodes = ["1.1.1.1", "9.9.9.9"]
    copy_ec_to_gem(model)
    assert model.reactions.get_by_id("r1").annotation["ec-code"] == ["1.1.1.1"]


# --------------------------------------------------------------------------- #
# gecko_light prefix handling
# --------------------------------------------------------------------------- #

def test_gecko_light_strips_4char_prefix_to_find_reaction():
    model = _make_ec_model(
        [("r1", None), ("r2", None)],
        gecko_light=True,
        ec_rxn_prefix="001_",
        ec_eccodes=["1.1.1.1", "2.2.2.2"],
    )
    # ec.rxns is now ["001_r1", "001_r2"]; cobra reactions are r1, r2.
    copy_ec_to_gem(model)
    assert model.reactions.get_by_id("r1").annotation["ec-code"] == ["1.1.1.1"]
    assert model.reactions.get_by_id("r2").annotation["ec-code"] == ["2.2.2.2"]


# --------------------------------------------------------------------------- #
# Round-trip with fill_eccodes_from_gem
# --------------------------------------------------------------------------- #

def test_round_trip_with_fill_eccodes_from_gem_is_stable():
    """copy_ec_to_gem (writes list) followed by fill_eccodes_from_gem (reads
    list and joins with `;`) yields the same ec.eccodes string."""
    model = _make_ec_model(
        [("r1", None), ("r2", None)],
        ec_eccodes=["1.1.1.1;2.2.2.2", "3.3.3.3"],
    )
    original = list(model.ec.eccodes)
    copy_ec_to_gem(model)
    # Wipe ec.eccodes; expect fill_eccodes_from_gem to reconstitute it.
    model.ec.eccodes = ["", ""]
    fill_eccodes_from_gem(model)
    assert model.ec.eccodes == original
