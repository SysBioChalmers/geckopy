"""Tests for apply_kcat_constraints."""
from pathlib import Path

import cobra
import pytest

from geckopy import EcModel, ModelAdapter, make_ec_model
from geckopy.ec_model.pipeline import apply_kcat_constraints

EXAMPLE_DIR = Path(__file__).parents[1] / "examples" / "ecTestGEM"


# --------------------------------------------------------------------------- #
# Fixtures
# --------------------------------------------------------------------------- #

_ECTESTGEM_CACHE: EcModel | None = None


def _ectestgem_ec_model() -> EcModel:
    """Load the ecTestGEM fixture through make_ec_model.

    The build is cached at module scope; each call returns a fresh
    deep copy so tests can mutate freely without leaking into each
    other. Building from scratch is ~3 s; deep-copy is <0.1 s.
    """
    import copy as _copy
    global _ECTESTGEM_CACHE
    if _ECTESTGEM_CACHE is None:
        adapter = ModelAdapter.from_folder(EXAMPLE_DIR)
        cobra_model = cobra.io.read_sbml_model(str(adapter.params.conv_gem))
        _ECTESTGEM_CACHE = make_ec_model(cobra_model, adapter)
    return _copy.deepcopy(_ECTESTGEM_CACHE)


def _get_s_coef(rxn: cobra.Reaction, met_id: str) -> float:
    """Get the stoichiometric coefficient of met_id in rxn, or 0 if absent."""
    for m, c in rxn.metabolites.items():
        if m.id == met_id:
            return c
    return 0.0


# --------------------------------------------------------------------------- #
# Basic behavior
# --------------------------------------------------------------------------- #

def test_writes_expected_coefficient_single_enzyme():
    """For a reaction with one enzyme (R3 / G4 / P4 in ecTestGEM), the
    coefficient should be -(1 * MW / (kcat * 3600)).
    R3 has one enzyme P4 with MW=40000 Da. Set kcat=10 /s."""
    ec_model = _ectestgem_ec_model()
    r3_idx = ec_model.ec.rxns.index("R3")
    ec_model.ec.kcat[r3_idx] = 10.0

    apply_kcat_constraints(ec_model)

    r3 = ec_model.reactions.get_by_id("R3")
    coef = _get_s_coef(r3, "prot_P4")
    expected = -(1.0 * 40000.0 / (10.0 * 3600.0))
    assert coef == pytest.approx(expected)


def test_writes_expected_coefficients_multi_subunit_reaction():
    """R2_EXP_1 has GPR 'G1 and G2' -> two enzymes P1 (MW=10000) and
    P2 (MW=20000), each with subunit count 1. Set kcat=100 /s."""
    ec_model = _ectestgem_ec_model()
    r_idx = ec_model.ec.rxns.index("R2_EXP_1")
    ec_model.ec.kcat[r_idx] = 100.0

    apply_kcat_constraints(ec_model)

    r = ec_model.reactions.get_by_id("R2_EXP_1")
    coef_p1 = _get_s_coef(r, "prot_P1")
    coef_p2 = _get_s_coef(r, "prot_P2")
    assert coef_p1 == pytest.approx(-(10000.0 / (100.0 * 3600.0)))
    assert coef_p2 == pytest.approx(-(20000.0 / (100.0 * 3600.0)))


def test_kcat_nan_produces_zero_coefficient():
    """When kcat is NaN, no coefficient is written. The all-NaN warning
    is expected: a freshly built ecModel has all-NaN kcats."""
    ec_model = _ectestgem_ec_model()
    # All kcats are initialized to NaN by make_ec_model.
    with pytest.warns(UserWarning, match="no valid entries"):
        apply_kcat_constraints(ec_model)
    r3 = ec_model.reactions.get_by_id("R3")
    assert _get_s_coef(r3, "prot_P4") == 0.0


def test_kcat_zero_produces_zero_coefficient():
    """All-zero/NaN selection triggers the documented warning."""
    ec_model = _ectestgem_ec_model()
    r3_idx = ec_model.ec.rxns.index("R3")
    ec_model.ec.kcat[r3_idx] = 0.0

    with pytest.warns(UserWarning, match="no valid entries"):
        apply_kcat_constraints(ec_model)

    r3 = ec_model.reactions.get_by_id("R3")
    assert _get_s_coef(r3, "prot_P4") == 0.0


def test_all_nan_kcats_warns_and_leaves_model_unchanged():
    ec_model = _ectestgem_ec_model()

    with pytest.warns(UserWarning, match="no valid entries"):
        apply_kcat_constraints(ec_model)

    # No reaction should gain any prot_ coefficient (except prot_pool,
    # which is handled by usage reactions, not apply_kcat).
    for rxn in ec_model.reactions:
        for met, coef in rxn.metabolites.items():
            if met.id.startswith("prot_") and met.id != "prot_pool":
                # Either 0 (if apply_kcat was about to write) or the
                # +1/-1 written by stage 11 for usage reactions.
                assert coef in (0.0, 1.0, -1.0), (
                    f"unexpected prot_ coef {coef} for {met.id} in {rxn.id}"
                )


def test_idempotent_when_applied_twice():
    ec_model = _ectestgem_ec_model()
    r3_idx = ec_model.ec.rxns.index("R3")
    ec_model.ec.kcat[r3_idx] = 10.0

    apply_kcat_constraints(ec_model)
    coef_first = _get_s_coef(ec_model.reactions.get_by_id("R3"), "prot_P4")

    apply_kcat_constraints(ec_model)
    coef_second = _get_s_coef(ec_model.reactions.get_by_id("R3"), "prot_P4")

    assert coef_first == coef_second


def test_changing_kcat_and_reapplying_overwrites():
    """Applying a kcat, then changing it and applying again, produces
    the new coefficient, not the sum or average."""
    ec_model = _ectestgem_ec_model()
    r3_idx = ec_model.ec.rxns.index("R3")
    ec_model.ec.kcat[r3_idx] = 10.0
    apply_kcat_constraints(ec_model)

    ec_model.ec.kcat[r3_idx] = 20.0
    apply_kcat_constraints(ec_model)

    expected = -(40000.0 / (20.0 * 3600.0))
    coef = _get_s_coef(ec_model.reactions.get_by_id("R3"), "prot_P4")
    assert coef == pytest.approx(expected)


# --------------------------------------------------------------------------- #
# update_rxns parameter
# --------------------------------------------------------------------------- #

def test_update_rxns_applies_only_to_specified():
    """Setting update_rxns=['R3'] should update R3 but not R5."""
    ec_model = _ectestgem_ec_model()
    ec_model.ec.kcat[ec_model.ec.rxns.index("R3")] = 10.0
    ec_model.ec.kcat[ec_model.ec.rxns.index("R5")] = 10.0

    apply_kcat_constraints(ec_model, update_rxns=["R3"])

    r3_coef = _get_s_coef(ec_model.reactions.get_by_id("R3"), "prot_P4")
    r5_coef = _get_s_coef(ec_model.reactions.get_by_id("R5"), "prot_P5")
    assert r3_coef != 0.0   # updated
    assert r5_coef == 0.0   # not updated


def test_update_rxns_unknown_id_raises():
    ec_model = _ectestgem_ec_model()
    with pytest.raises(ValueError, match="not present in ec.rxns"):
        apply_kcat_constraints(ec_model, update_rxns=["not_a_rxn"])


def test_update_rxns_empty_list_is_noop():
    ec_model = _ectestgem_ec_model()
    ec_model.ec.kcat[ec_model.ec.rxns.index("R3")] = 10.0

    apply_kcat_constraints(ec_model, update_rxns=[])

    r3_coef = _get_s_coef(ec_model.reactions.get_by_id("R3"), "prot_P4")
    assert r3_coef == 0.0


def test_setting_kcat_to_nan_and_reapplying_clears_old_constraint():
    """Apply a kcat, then set it back to NaN and re-apply: the previous
    coefficient must be cleared, even though the warning about all-NaN
    selection still fires."""
    ec_model = _ectestgem_ec_model()
    r3_idx = ec_model.ec.rxns.index("R3")
    ec_model.ec.kcat[r3_idx] = 10.0
    apply_kcat_constraints(ec_model)
    assert _get_s_coef(
        ec_model.reactions.get_by_id("R3"), "prot_P4"
    ) != 0.0

    ec_model.ec.kcat[r3_idx] = float("nan")
    with pytest.warns(UserWarning):
        apply_kcat_constraints(ec_model, update_rxns=["R3"])

    # Old coefficient cleared.
    assert _get_s_coef(
        ec_model.reactions.get_by_id("R3"), "prot_P4"
    ) == 0.0


def test_setting_kcat_to_zero_and_reapplying_clears_old_constraint():
    """Same as the NaN case but with 0.0 instead. Both should clear."""
    ec_model = _ectestgem_ec_model()
    r3_idx = ec_model.ec.rxns.index("R3")
    ec_model.ec.kcat[r3_idx] = 10.0
    apply_kcat_constraints(ec_model)
    assert _get_s_coef(
        ec_model.reactions.get_by_id("R3"), "prot_P4"
    ) != 0.0

    ec_model.ec.kcat[r3_idx] = 0.0
    with pytest.warns(UserWarning):
        apply_kcat_constraints(ec_model, update_rxns=["R3"])

    assert _get_s_coef(
        ec_model.reactions.get_by_id("R3"), "prot_P4"
    ) == 0.0


# --------------------------------------------------------------------------- #
# Gecko-light not yet implemented
# --------------------------------------------------------------------------- #

def test_gecko_light_raises():
    ec_model = _ectestgem_ec_model()
    ec_model.ec.gecko_light = True
    with pytest.raises(NotImplementedError, match="gecko-light"):
        apply_kcat_constraints(ec_model)


# --------------------------------------------------------------------------- #
# Multi-subunit (rxn_enz_mat entry > 1)
# --------------------------------------------------------------------------- #

def test_subunit_multiplier_applied():
    """If rxn_enz_mat[r, e] = 3, the coefficient magnitude triples."""
    ec_model = _ectestgem_ec_model()

    r3_idx = ec_model.ec.rxns.index("R3")
    p4_idx = ec_model.ec.genes.index("G4")

    # Edit rxn_enz_mat to set R3 <-> P4 at 3 subunits.
    mat = ec_model.ec.rxn_enz_mat.tolil()
    mat[r3_idx, p4_idx] = 3.0
    ec_model.ec.rxn_enz_mat = mat.tocsr()

    ec_model.ec.kcat[r3_idx] = 10.0
    apply_kcat_constraints(ec_model)

    r3 = ec_model.reactions.get_by_id("R3")
    coef = _get_s_coef(r3, "prot_P4")
    expected = -(3.0 * 40000.0 / (10.0 * 3600.0))
    assert coef == pytest.approx(expected)


# --------------------------------------------------------------------------- #
# End-to-end: apply to all ecTestGEM reactions with a uniform kcat
# --------------------------------------------------------------------------- #

def test_apply_to_all_ectestgem_reactions_with_uniform_kcat():
    """Set every catalyzed reaction's kcat to 100 /s and apply. Every
    (rxn, enzyme) pair listed in rxn_enz_mat should yield a negative
    coefficient on the corresponding prot_ metabolite."""
    ec_model = _ectestgem_ec_model()
    ec_model.ec.kcat[:] = 100.0

    apply_kcat_constraints(ec_model)

    mat = ec_model.ec.rxn_enz_mat.tocsr()
    for i, rxn_id in enumerate(ec_model.ec.rxns):
        row = mat.getrow(i)
        rxn = ec_model.reactions.get_by_id(rxn_id)
        for enz_idx in row.indices:
            enzyme = ec_model.ec.enzymes[enz_idx]
            coef = _get_s_coef(rxn, f"prot_{enzyme}")
            assert coef < 0.0, (
                f"expected negative prot_{enzyme} coefficient for {rxn_id}"
            )


def test_model_is_still_valid_cobra_model_after_apply():
    """Sanity: after applying kcats, the model should still solve
    a simple FBA (not necessarily to a meaningful objective)."""
    ec_model = _ectestgem_ec_model()
    ec_model.ec.kcat[:] = 100.0
    apply_kcat_constraints(ec_model)

    # Just verify it optimizes without error; optimum may be 0.
    sol = ec_model.optimize()
    assert sol.status in ("optimal", "infeasible")
