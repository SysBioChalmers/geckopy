"""Tests for get_conc_control_coeffs."""
import cobra
import numpy as np
import pytest
from scipy import sparse

from geckopy.ec_model import EcModel
from geckopy.ec_model.ec_data import EcData
from geckopy.limit_proteins import get_conc_control_coeffs


# --------------------------------------------------------------------------- #
# Tiny enzyme-constrained model fixture
# --------------------------------------------------------------------------- #

def _build_toy_ec_model(
    *,
    enzyme_ub: float = 5.0,
    enzyme_coeff_in_R: float = -1.0,
) -> EcModel:
    """Build a 1-substrate-1-enzyme model with a measurable growth rate.

    Topology:
        EX_A:               A_e <->                  (substrate exchange)
        TR_A:               A_e -> A_c               (transport)
        R_AB:    A_c + (-coeff)*prot_E -> B_c        (enzyme-constrained)
        SK_B:               B_c ->                   (sink, biomass-like)
        prot_pool_exchange: prot_pool ->             (pool source)
        usage_prot_E:       prot_pool -> prot_E      (usage; ub = enzyme_ub)

    The objective is SK_B. Maximum SK_B flux equals max usage flux,
    which equals enzyme_ub when enzyme_coeff_in_R = -1 (one unit of
    enzyme per unit of R flux).
    """
    model = EcModel("toy")

    A_e = cobra.Metabolite("A_e", compartment="e")
    A_c = cobra.Metabolite("A_c", compartment="c")
    B_c = cobra.Metabolite("B_c", compartment="c")
    pool = cobra.Metabolite("prot_pool", compartment="c")
    prot_E = cobra.Metabolite("prot_E", compartment="c")
    model.add_metabolites([A_e, A_c, B_c, pool, prot_E])

    EX_A = cobra.Reaction("EX_A")
    EX_A.add_metabolites({A_e: -1.0})
    EX_A.lower_bound = -1000.0
    EX_A.upper_bound = 0.0

    TR_A = cobra.Reaction("TR_A")
    TR_A.add_metabolites({A_e: -1.0, A_c: 1.0})
    TR_A.lower_bound = 0.0
    TR_A.upper_bound = 1000.0

    R_AB = cobra.Reaction("R_AB")
    R_AB.add_metabolites({A_c: -1.0, B_c: 1.0, prot_E: enzyme_coeff_in_R})
    R_AB.lower_bound = 0.0
    R_AB.upper_bound = 1000.0

    SK_B = cobra.Reaction("SK_B")
    SK_B.add_metabolites({B_c: -1.0})
    SK_B.lower_bound = 0.0
    SK_B.upper_bound = 1000.0

    pool_ex = cobra.Reaction("prot_pool_exchange")
    pool_ex.add_metabolites({pool: 1.0})
    pool_ex.lower_bound = 0.0
    pool_ex.upper_bound = 1000.0

    usage = cobra.Reaction("usage_prot_E")
    usage.add_metabolites({pool: -1.0, prot_E: 1.0})
    usage.lower_bound = 0.0
    usage.upper_bound = enzyme_ub

    model.add_reactions([EX_A, TR_A, R_AB, SK_B, pool_ex, usage])
    model.objective = "SK_B"

    model.ec = EcData(
        rxns=["R_AB"],
        kcat=np.array([1.0]),
        source=[""],
        notes=[""],
        eccodes=[""],
        genes=["g_E"],
        enzymes=["E"],
        mw=np.array([100.0]),
        sequence=[""],
        concs=np.array([float(enzyme_ub)]),
        rxn_enz_mat=sparse.csr_matrix([[1.0]]),
    )
    return model


# --------------------------------------------------------------------------- #
# Trivial cases
# --------------------------------------------------------------------------- #

def test_empty_proteins_list_returns_empty():
    model = _build_toy_ec_model()
    enz, coeffs = get_conc_control_coeffs(model, proteins=[])
    assert enz.size == 0
    assert coeffs.size == 0


def test_default_proteins_uses_model_ec_enzymes():
    model = _build_toy_ec_model()
    enz, coeffs = get_conc_control_coeffs(model)  # default proteins
    assert len(enz) == 1
    assert len(coeffs) == 1


def test_protein_not_in_model_silently_skipped():
    model = _build_toy_ec_model()
    enz, coeffs = get_conc_control_coeffs(
        model, proteins=["E", "ghost"],
    )
    assert enz[0]
    assert not enz[1]
    assert coeffs[1] == 0.0


# --------------------------------------------------------------------------- #
# Sensitivity analysis: limiting enzyme
# --------------------------------------------------------------------------- #

def test_limiting_enzyme_marked_and_has_positive_coeff():
    """The single enzyme E directly limits R_AB; relaxing it must
    yield a positive control coefficient."""
    model = _build_toy_ec_model(enzyme_ub=5.0)
    enz, coeffs = get_conc_control_coeffs(model)
    assert enz[0]
    assert coeffs[0] > 0


def test_control_coeff_matches_growth_per_unit_increase():
    """With enzyme_coeff_in_R = -1, growth = enzyme ub. Doubling ub
    from 5 to 10 increases growth by 5; coefficient = 5/(10-5) = 1."""
    model = _build_toy_ec_model(enzyme_ub=5.0)
    _, coeffs = get_conc_control_coeffs(
        model, fold_change=2.0,
    )
    assert coeffs[0] == pytest.approx(1.0, rel=1e-6)


def test_control_coeff_scales_with_enzyme_coefficient():
    """If R_AB consumes 2 units of enzyme per unit flux, doubling the
    enzyme ub from 5 to 10 increases growth from 2.5 to 5 (delta=2.5),
    and the coeff is 2.5/(10-5) = 0.5."""
    model = _build_toy_ec_model(enzyme_ub=5.0, enzyme_coeff_in_R=-2.0)
    _, coeffs = get_conc_control_coeffs(model, fold_change=2.0)
    assert coeffs[0] == pytest.approx(0.5, rel=1e-6)


# --------------------------------------------------------------------------- #
# limit parameter
# --------------------------------------------------------------------------- #

def test_default_limit_zero_includes_active_enzyme():
    """Default limit=0 -> any protein with non-zero usage qualifies."""
    model = _build_toy_ec_model(enzyme_ub=5.0)
    enz, _ = get_conc_control_coeffs(model, limit=0.0)
    assert enz[0]


def test_high_limit_excludes_under_used_enzyme():
    """Set ub very high so usage/ub is small; with limit=0.5, the
    enzyme should be skipped."""
    model = _build_toy_ec_model(enzyme_ub=5.0)
    # Cap upstream so the enzyme isn't used at full capacity:
    # restrict EX_A to deliver only 1 unit/time, so usage flux = 1.
    model.reactions.get_by_id("EX_A").lower_bound = -1.0
    enz, coeffs = get_conc_control_coeffs(model, limit=0.5)
    assert not enz[0]
    assert coeffs[0] == 0.0


# --------------------------------------------------------------------------- #
# fold_change parameter
# --------------------------------------------------------------------------- #

def test_fold_change_affects_probe_size():
    """Coefficient = delta_growth / delta_ub. With enzyme_coeff_in_R=-1,
    delta_growth = delta_ub, so the coefficient is 1.0 regardless of
    fold_change. Verify both default and a different value."""
    model = _build_toy_ec_model(enzyme_ub=5.0)
    _, c2 = get_conc_control_coeffs(model, fold_change=2.0)
    _, c5 = get_conc_control_coeffs(model, fold_change=5.0)
    assert c2[0] == pytest.approx(1.0, rel=1e-6)
    assert c5[0] == pytest.approx(1.0, rel=1e-6)


# --------------------------------------------------------------------------- #
# Reversibility: model state restored after each probe
# --------------------------------------------------------------------------- #

def test_model_bounds_unchanged_after_call():
    model = _build_toy_ec_model(enzyme_ub=5.0)
    original_ub = model.reactions.get_by_id("usage_prot_E").upper_bound
    get_conc_control_coeffs(model)
    assert model.reactions.get_by_id("usage_prot_E").upper_bound == original_ub


# --------------------------------------------------------------------------- #
# Subset of proteins
# --------------------------------------------------------------------------- #

def test_subset_proteins_only_those_evaluated():
    model = _build_toy_ec_model(enzyme_ub=5.0)
    enz, coeffs = get_conc_control_coeffs(model, proteins=["E"])
    assert len(enz) == 1
    assert enz[0]
    assert coeffs[0] > 0


# --------------------------------------------------------------------------- #
# Single-solve alternative (validating the review §3.4 suggestion)
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize(
    "ub, coeff_in_R, expected",
    [(5.0, -1.0, 1.0), (3.0, -1.0, 1.0), (6.0, -2.0, 0.5)],
)
def test_shadow_price_matches_finite_difference_coeff(ub, coeff_in_R, expected):
    """The per-protein finite-difference loop can be replaced by a single
    solve: the prot_<id> metabolite shadow price (|dual|) equals the
    finite-difference control coefficient.

    Note: the *reaction* ``reduced_cost`` of the usage reaction does NOT —
    it is scaled by the number of coupled pseudometabolites (prot_pool +
    prot_<id>), so it reads 2x here. Use the metabolite shadow price.
    """
    model = _build_toy_ec_model(enzyme_ub=ub, enzyme_coeff_in_R=coeff_in_R)
    _, coeffs = get_conc_control_coeffs(model)  # finite-difference
    model.optimize()  # single solve
    shadow = abs(model.enzymes.get_by_id("E").shadow_price)
    assert coeffs[0] == pytest.approx(expected)
    assert shadow == pytest.approx(coeffs[0], rel=1e-6)
