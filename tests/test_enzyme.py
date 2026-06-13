"""Tests for the Enzyme / Kcats / EnzymeView proxy classes."""
from pathlib import Path

import cobra
import pytest

from geckopy import EcModel, ModelAdapter, make_ec_model
from geckopy.ec_model import Enzyme
from geckopy.ec_model.pipeline import apply_kcat_constraints

EXAMPLE_DIR = Path(__file__).parents[1] / "examples" / "ecTestGEM"


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


def _get_s_coef(rxn: cobra.Reaction, met_id: str) -> float:
    for m, c in rxn.metabolites.items():
        if m.id == met_id:
            return c
    return 0.0


# --------------------------------------------------------------------------- #
# Identity / live proxy
# --------------------------------------------------------------------------- #

def test_get_by_id_returns_live_proxy():
    model = _ectestgem_ec_model()
    enz = model.enzymes.get_by_id("P4")
    assert isinstance(enz, Enzyme)
    assert enz.id == "P4"

    idx = model.ec.enzymes.index("P4")
    model.ec.concs[idx] = 0.012
    assert enz.concentration == pytest.approx(0.012)


def test_unknown_uniprot_raises():
    model = _ectestgem_ec_model()
    with pytest.raises(KeyError):
        model.enzymes.get_by_id("P_ghost")


def test_iter_and_len():
    model = _ectestgem_ec_model()
    assert len(model.enzymes) == len(model.ec.enzymes)
    ids = [enz.id for enz in model.enzymes]
    assert ids == model.ec.enzymes


def test_contains():
    model = _ectestgem_ec_model()
    assert "P4" in model.enzymes
    assert "P_ghost" not in model.enzymes


# --------------------------------------------------------------------------- #
# Setters
# --------------------------------------------------------------------------- #

def test_concentration_setter_updates_usage_bound():
    model = _ectestgem_ec_model()
    enz = model.enzymes.get_by_id("P4")
    enz.concentration = 1e-3
    assert enz.usage_reaction.upper_bound == pytest.approx(1e-3)


def test_concentration_nan_returns_to_pool():
    model = _ectestgem_ec_model()
    enz = model.enzymes.get_by_id("P4")
    enz.concentration = 1e-3
    enz.concentration = float("nan")
    assert enz.usage_reaction.upper_bound == pytest.approx(1000.0)


def test_mw_setter_updates_all_catalyzed_reactions():
    """Set MW on P4, assert the coefficient on R3 (which is catalysed
    by P4) is recomputed using the new MW."""
    model = _ectestgem_ec_model()
    r3_idx = model.ec.rxns.index("R3")
    model.ec.kcat[r3_idx] = 10.0
    apply_kcat_constraints(model)

    enz = model.enzymes.get_by_id("P4")
    enz.mw = 50000.0
    coef = _get_s_coef(model.reactions.get_by_id("R3"), "prot_P4")
    expected = -(1.0 * 50000.0 / (10.0 * 3600.0))
    assert coef == pytest.approx(expected)


# --------------------------------------------------------------------------- #
# Kcats
# --------------------------------------------------------------------------- #

def test_kcats_get_returns_ec_kcat():
    model = _ectestgem_ec_model()
    r3_idx = model.ec.rxns.index("R3")
    model.ec.kcat[r3_idx] = 7.5
    enz = model.enzymes.get_by_id("P4")
    assert enz.kcats["R3"] == pytest.approx(7.5)


def test_kcats_set_single_enzyme_reaction():
    model = _ectestgem_ec_model()
    enz = model.enzymes.get_by_id("P4")
    enz.kcats["R3"] = 25.0
    expected = -(1.0 * 40000.0 / (25.0 * 3600.0))
    coef = _get_s_coef(model.reactions.get_by_id("R3"), "prot_P4")
    assert coef == pytest.approx(expected)


def test_kcats_set_complex_raises():
    """R2_EXP_1 is catalysed by P1+P2; per-enzyme kcat is ambiguous."""
    model = _ectestgem_ec_model()
    enz = model.enzymes.get_by_id("P1")
    with pytest.raises(ValueError, match="catalysed by"):
        enz.kcats["R2_EXP_1"] = 100.0


def test_kcats_unknown_rxn_raises():
    model = _ectestgem_ec_model()
    enz = model.enzymes.get_by_id("P4")
    with pytest.raises(KeyError):
        enz.kcats["R_ghost"]


def test_kcats_iter_lists_catalysed_reactions():
    model = _ectestgem_ec_model()
    enz = model.enzymes.get_by_id("P4")
    expected_ids = {
        model.ec.rxns[i]
        for i in model.ec.rxn_enz_mat.tocsc().getcol(enz.index).nonzero()[0]
    }
    assert set(enz.kcats) == expected_ids


# --------------------------------------------------------------------------- #
# Solver-side reads
# --------------------------------------------------------------------------- #

def test_flux_after_optimize():
    """Set a kcat so apply_kcat_constraints writes a real coefficient,
    pick the model's only exchange as objective, then solve."""
    model = _ectestgem_ec_model()
    r3_idx = model.ec.rxns.index("R3")
    model.ec.kcat[r3_idx] = 10.0
    apply_kcat_constraints(model)
    # Pick an objective that consumes P4 via R3.
    model.objective = "R3"
    sol = model.optimize()
    assert sol.status == "optimal"
    enz = model.enzymes.get_by_id("P4")
    expected = model.reactions.get_by_id("usage_prot_P4").flux
    assert enz.flux == pytest.approx(expected)


def test_reactions_matches_rxn_enz_mat():
    model = _ectestgem_ec_model()
    enz = model.enzymes.get_by_id("P4")
    col = model.ec.rxn_enz_mat.tocsc().getcol(enz.index)
    expected_ids = {model.ec.rxns[i] for i in col.nonzero()[0]}
    actual_ids = {r.id for r in enz.reactions}
    assert actual_ids == expected_ids


# --------------------------------------------------------------------------- #
# Identity properties
# --------------------------------------------------------------------------- #

def test_prot_metabolite_and_usage_reaction_ids():
    model = _ectestgem_ec_model()
    enz = model.enzymes.get_by_id("P4")
    assert enz.prot_metabolite_id == "prot_P4"
    assert enz.usage_reaction_id == "usage_prot_P4"
    assert enz.prot_metabolite.id == "prot_P4"
    assert enz.usage_reaction.id == "usage_prot_P4"


def test_gene_property():
    model = _ectestgem_ec_model()
    enz = model.enzymes.get_by_id("P4")
    assert enz.gene == "G4"


# --------------------------------------------------------------------------- #
# Query
# --------------------------------------------------------------------------- #

def test_query_by_gene_prefix():
    model = _ectestgem_ec_model()
    results = model.enzymes.query(lambda g: g.startswith("G1"), "gene")
    assert {e.gene for e in results} == {"G1"}


def test_query_by_mw():
    model = _ectestgem_ec_model()
    results = model.enzymes.query(lambda mw: mw >= 40000, "mw")
    assert all(e.mw >= 40000 for e in results)


# --------------------------------------------------------------------------- #
# gecko-light guards
# --------------------------------------------------------------------------- #

def test_gecko_light_blocks_concentration_setter():
    model = _ectestgem_ec_model()
    model.ec.gecko_light = True
    enz = model.enzymes.get_by_id("P4")
    with pytest.raises(NotImplementedError):
        enz.concentration = 1e-3


def test_gecko_light_blocks_kcats_setter():
    model = _ectestgem_ec_model()
    model.ec.gecko_light = True
    enz = model.enzymes.get_by_id("P4")
    with pytest.raises(NotImplementedError):
        enz.kcats["R3"] = 10.0


def test_gecko_light_allows_read_only_metadata():
    """Read-only metadata (mw, gene, sequence, concentration, kcats[r])
    works on light models because it reads from `model.ec.*` arrays
    that exist in both layouts. Only the per-enzyme prot/usage
    machinery is missing."""
    model = _ectestgem_ec_model()
    model.ec.gecko_light = True
    enz = model.enzymes.get_by_id("P4")
    # All of these should work without raising.
    assert enz.id == "P4"
    assert enz.gene == "G4"
    assert isinstance(enz.mw, float)
    assert enz.sequence == "MDFM"


def test_gecko_light_blocks_full_model_only_attrs():
    """Attributes that require the per-enzyme `prot_<id>` met and
    `usage_prot_<id>` reaction (which only the full layout has)
    raise a clear NotImplementedError on a light model."""
    model = _ectestgem_ec_model()
    model.ec.gecko_light = True
    enz = model.enzymes.get_by_id("P4")
    with pytest.raises(NotImplementedError, match="gecko-light"):
        _ = enz.prot_metabolite
    with pytest.raises(NotImplementedError, match="gecko-light"):
        _ = enz.usage_reaction
    with pytest.raises(NotImplementedError, match="gecko-light"):
        _ = enz.shadow_price


# --------------------------------------------------------------------------- #
# repr smoke
# --------------------------------------------------------------------------- #

def test_repr_html_smoke():
    """The HTML repr falls back to '-' when no flux is cached."""
    model = _ectestgem_ec_model()
    enz = model.enzymes.get_by_id("P4")
    html = enz._repr_html_()
    assert "P4" in html
    assert "G4" in html
    # No solve has happened -> flux/cap_usage fields are filled with '-'.
    assert ">-<" in html


def test_repr_smoke():
    model = _ectestgem_ec_model()
    enz = model.enzymes.get_by_id("P4")
    s = repr(enz)
    assert "P4" in s
    assert "G4" in s


# --------------------------------------------------------------------------- #
# Real gecko-light builds (not a flipped flag — actual ec.rxns shape)
# --------------------------------------------------------------------------- #

_ECTESTGEM_LIGHT_CACHE: EcModel | None = None


def _ectestgem_light_ec_model() -> EcModel:
    """Cached gecko-light build of ecTestGEM; deep-copied per call."""
    import copy as _copy
    global _ECTESTGEM_LIGHT_CACHE
    if _ECTESTGEM_LIGHT_CACHE is None:
        adapter = ModelAdapter.from_folder(EXAMPLE_DIR)
        cobra_model = cobra.io.read_sbml_model(str(adapter.params.conv_gem))
        _ECTESTGEM_LIGHT_CACHE = make_ec_model(
            cobra_model, adapter, gecko_light=True,
        )
    return _copy.deepcopy(_ECTESTGEM_LIGHT_CACHE)


def test_light_reactions_returns_cobra_reactions():
    """On a real light build, ec.rxns ids are prefixed (``001_R3``) but
    cobra has the unprefixed reaction. ``Enzyme.reactions`` must strip
    the prefix so it returns the real cobra reaction, not an empty set.
    """
    model = _ectestgem_light_ec_model()
    enz = model.enzymes.get_by_id("P4")  # G4 catalyses R3
    rxn_ids = {r.id for r in enz.reactions}
    assert "R3" in rxn_ids
    # And the entries are real cobra.Reaction instances.
    assert all(isinstance(r, cobra.Reaction) for r in enz.reactions)


def test_light_reactions_dedupes_isozyme_rows():
    """G1 is part of R2's complex isozyme (``001_R2`` AND ``001_R2_REV``).
    Both ec rows point to the same pair of cobra reactions, so the
    returned set must contain each cobra reaction once."""
    model = _ectestgem_light_ec_model()
    enz = model.enzymes.get_by_id("P1")
    cobra_ids = [r.id for r in enz.reactions]
    assert len(cobra_ids) == len(set(cobra_ids)), (
        f"Duplicate cobra reactions in Enzyme.reactions: {cobra_ids}"
    )
    assert set(cobra_ids) == {"R2", "R2_REV"}


def test_light_mw_setter_reapplies_constraints():
    """Changing MW on a light build must re-apply prot_pool coefficients
    for every reaction the enzyme participates in. Uses ec.rxns ids
    (the ``###_`` prefixed form) internally."""
    model = _ectestgem_light_ec_model()
    # Set a kcat on 001_R3 (the single isozyme of R3 = G4) and apply.
    r3_idx = model.ec.rxns.index("001_R3")
    model.ec.kcat[r3_idx] = 10.0
    apply_kcat_constraints(model)

    enz = model.enzymes.get_by_id("P4")
    enz.mw = 50000.0

    coef = _get_s_coef(model.reactions.get_by_id("R3"), "prot_pool")
    # cheapest isozyme cost: -MW_sum / (kcat * 3600), MW_sum = 50000.
    expected = -(50000.0 / (10.0 * 3600.0))
    assert coef == pytest.approx(expected)


def test_light_repr_html_does_not_crash():
    """The HTML repr accesses flux/cap_usage; on light those raise
    NotImplementedError, which the repr now catches alongside
    RuntimeError. Verify it returns a string without crashing."""
    model = _ectestgem_light_ec_model()
    enz = model.enzymes.get_by_id("P4")
    html = enz._repr_html_()
    assert "P4" in html
    assert ">-<" in html  # flux + cap_usage shown as '-'


def test_light_kcats_keys_are_prefixed_ec_rxn_ids():
    """Light kcats are keyed by ec.rxns ids (which carry the ``###_``
    prefix in light), not by cobra reaction ids. Different isozymes
    can hold different kcat values, so the prefix has to stay."""
    model = _ectestgem_light_ec_model()
    enz = model.enzymes.get_by_id("P4")  # only catalyses one row: 001_R3
    assert list(enz.kcats) == ["001_R3"]
    model.ec.kcat[model.ec.rxns.index("001_R3")] = 7.5
    assert enz.kcats["001_R3"] == pytest.approx(7.5)
