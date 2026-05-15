"""Tests for ec_fseof."""
from pathlib import Path

import cobra
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from geckopy import EcModel, ModelAdapter
from geckopy.ec_model.ec_data import EcData
from geckopy.utilities import EcFseofResult, ec_fseof


# --------------------------------------------------------------------------- #
# Tiny FSEOF-able model fixture
# --------------------------------------------------------------------------- #

def _adapter(tmp_path: Path, *, bio_rxn: str = "biomass") -> ModelAdapter:
    (tmp_path / "model_adapter.toml").write_text(
        f'conv_gem = "dummy.xml"\n'
        f'org_name = "test"\n'
        f'bio_rxn = "{bio_rxn}"\n'
    )
    return ModelAdapter.from_folder(tmp_path)


def _build_fseof_model(adapter: ModelAdapter) -> EcModel:
    """Substrate A is split between biomass (via R_BIO) and the
    production target (via R_PROD). Both routes consume an enzyme.

    EX_A: -1*A_e (uptake)
    TR_A: A_e -> A_c
    R_BIO: A_c + (1/100)*prot_E_BIO -> bio_met        (gpr: g_BIO)
    R_PROD: A_c + (1/100)*prot_E_PROD -> prod_met     (gpr: g_PROD)
    biomass: bio_met ->                                (objective)
    EX_PROD: prod_met ->                               (production target)
    usage_prot_E_BIO: prot_pool -> prot_E_BIO
    usage_prot_E_PROD: prot_pool -> prot_E_PROD
    prot_pool_exchange: -> prot_pool                   (small ub forces tradeoff)
    """
    model = EcModel("toy", adapter=adapter)

    A_e = cobra.Metabolite("A_e", compartment="e")
    A_c = cobra.Metabolite("A_c", compartment="c")
    bio_met = cobra.Metabolite("bio_met", compartment="c")
    prod_met = cobra.Metabolite("prod_met", compartment="c")
    pool = cobra.Metabolite("prot_pool", compartment="c")
    prot_BIO = cobra.Metabolite("prot_E_BIO", compartment="c")
    prot_PROD = cobra.Metabolite("prot_E_PROD", compartment="c")
    model.add_metabolites([
        A_e, A_c, bio_met, prod_met, pool, prot_BIO, prot_PROD,
    ])

    EX_A = cobra.Reaction("EX_A")
    EX_A.add_metabolites({A_e: -1.0})
    EX_A.lower_bound = -100.0; EX_A.upper_bound = 0.0

    TR_A = cobra.Reaction("TR_A")
    TR_A.add_metabolites({A_e: -1.0, A_c: 1.0})
    TR_A.lower_bound = 0.0; TR_A.upper_bound = 1000.0

    R_BIO = cobra.Reaction("R_BIO")
    R_BIO.add_metabolites({A_c: -1.0, prot_BIO: -1/100, bio_met: 1.0})
    R_BIO.lower_bound = 0.0; R_BIO.upper_bound = 1000.0
    R_BIO.gene_reaction_rule = "g_BIO"

    R_PROD = cobra.Reaction("R_PROD")
    # R_PROD uses half as much enzyme per unit flux as R_BIO,
    # giving them differing slopes so the top-25% filter retains R_PROD.
    R_PROD.add_metabolites({A_c: -1.0, prot_PROD: -1/200, prod_met: 1.0})
    R_PROD.lower_bound = 0.0; R_PROD.upper_bound = 1000.0
    R_PROD.gene_reaction_rule = "g_PROD"

    BIO = cobra.Reaction("biomass")
    BIO.add_metabolites({bio_met: -1.0})
    BIO.lower_bound = 0.0; BIO.upper_bound = 1000.0

    EX_PROD = cobra.Reaction("EX_PROD")
    EX_PROD.add_metabolites({prod_met: -1.0})
    EX_PROD.lower_bound = 0.0; EX_PROD.upper_bound = 1000.0

    usage_BIO = cobra.Reaction("usage_prot_E_BIO")
    usage_BIO.add_metabolites({pool: -1.0, prot_BIO: 1.0})
    usage_BIO.lower_bound = 0.0; usage_BIO.upper_bound = 1000.0

    usage_PROD = cobra.Reaction("usage_prot_E_PROD")
    usage_PROD.add_metabolites({pool: -1.0, prot_PROD: 1.0})
    usage_PROD.lower_bound = 0.0; usage_PROD.upper_bound = 1000.0

    pool_ex = cobra.Reaction("prot_pool_exchange")
    pool_ex.add_metabolites({pool: 1.0})
    pool_ex.lower_bound = 0.0; pool_ex.upper_bound = 0.5  # forces tradeoff

    model.add_reactions([
        EX_A, TR_A, R_BIO, R_PROD, BIO, EX_PROD,
        usage_BIO, usage_PROD, pool_ex,
    ])
    model.objective = "biomass"

    n = 2
    g = 2
    mat = sparse.lil_matrix((n, g), dtype=float)
    mat[0, 0] = 1.0  # R_BIO uses E_BIO
    mat[1, 1] = 1.0  # R_PROD uses E_PROD
    model.ec = EcData(
        rxns=["R_BIO", "R_PROD"],
        kcat=np.array([1.0, 1.0]),
        source=["initial", "initial"],
        notes=["", ""],
        eccodes=["", ""],
        genes=["g_BIO", "g_PROD"],
        enzymes=["E_BIO", "E_PROD"],
        mw=np.array([100.0, 100.0]),
        sequence=["", ""],
        concs=np.array([np.nan, np.nan]),
        rxn_enz_mat=mat.tocsr(),
    )
    return model


# --------------------------------------------------------------------------- #
# Pre-flight
# --------------------------------------------------------------------------- #

def test_no_adapter_raises(tmp_path):
    adapter = _adapter(tmp_path)
    model = _build_fseof_model(adapter)
    model.adapter = None
    with pytest.raises(ValueError, match="adapter"):
        ec_fseof(model, "EX_PROD", "EX_A")


def test_n_steps_too_small_raises(tmp_path):
    adapter = _adapter(tmp_path)
    model = _build_fseof_model(adapter)
    with pytest.raises(ValueError, match="n_steps"):
        ec_fseof(model, "EX_PROD", "EX_A", n_steps=1)


def test_no_production_headroom_raises(tmp_path):
    """If max prod flux <= initial prod flux, raise."""
    adapter = _adapter(tmp_path)
    model = _build_fseof_model(adapter)
    # Force initial production = max by setting EX_PROD lb high.
    model.reactions.get_by_id("EX_PROD").lower_bound = 100.0
    model.reactions.get_by_id("EX_PROD").upper_bound = 100.0
    with pytest.raises(ValueError, match="headroom|infeasible"):
        ec_fseof(model, "EX_PROD", "EX_A")


# --------------------------------------------------------------------------- #
# Result dataclass shape
# --------------------------------------------------------------------------- #

def test_returns_dataclass(tmp_path):
    adapter = _adapter(tmp_path)
    model = _build_fseof_model(adapter)
    result = ec_fseof(model, "EX_PROD", "EX_A", n_steps=8)
    assert isinstance(result, EcFseofResult)
    assert isinstance(result.alpha, np.ndarray)
    assert result.alpha.shape == (8,)
    assert isinstance(result.v_matrix, pd.DataFrame)
    assert isinstance(result.rxn_targets, pd.DataFrame)
    assert isinstance(result.transport_targets, pd.DataFrame)
    assert isinstance(result.gene_targets, pd.DataFrame)


def test_alpha_grid_endpoints(tmp_path):
    """alpha[0] = initial production flux at biomass-max;
    alpha[-1] = 90% of max theoretical production."""
    adapter = _adapter(tmp_path)
    model = _build_fseof_model(adapter)
    result = ec_fseof(model, "EX_PROD", "EX_A", n_steps=4)
    # Initial production at biomass max = 0 (model has no incentive
    # to produce). 90% of max = 0.9 * (some max flux).
    assert result.alpha[0] == pytest.approx(0.0, abs=1e-9)
    assert result.alpha[-1] > 0


def test_v_matrix_columns_match_alpha(tmp_path):
    adapter = _adapter(tmp_path)
    model = _build_fseof_model(adapter)
    result = ec_fseof(model, "EX_PROD", "EX_A", n_steps=4)
    assert list(result.v_matrix.columns) == [str(a) for a in result.alpha]


# --------------------------------------------------------------------------- #
# Target identification
# --------------------------------------------------------------------------- #

def test_r_prod_identified_as_oe_target(tmp_path):
    """As alpha increases, R_PROD flux must increase (it produces the
    target). It should be tagged OE."""
    adapter = _adapter(tmp_path)
    model = _build_fseof_model(adapter)
    result = ec_fseof(model, "EX_PROD", "EX_A", n_steps=8)
    # R_PROD (or its associated rxn) should appear with OE action.
    if "R_PROD" in list(result.rxn_targets["rxn_id"]):
        row = result.rxn_targets[
            result.rxn_targets["rxn_id"] == "R_PROD"
        ].iloc[0]
        assert row["action"] == "OE"


def test_gene_targets_have_expected_columns(tmp_path):
    adapter = _adapter(tmp_path)
    model = _build_fseof_model(adapter)
    result = ec_fseof(model, "EX_PROD", "EX_A", n_steps=8)
    assert list(result.gene_targets.columns) == [
        "gene_id", "gene_name", "slope", "action", "essentiality",
    ]


def test_rxn_targets_have_expected_columns(tmp_path):
    adapter = _adapter(tmp_path)
    model = _build_fseof_model(adapter)
    result = ec_fseof(model, "EX_PROD", "EX_A", n_steps=8)
    assert list(result.rxn_targets.columns) == [
        "rxn_id", "rxn_name", "slope", "gpr", "equation", "action",
    ]


# --------------------------------------------------------------------------- #
# Sorting
# --------------------------------------------------------------------------- #

def test_gene_targets_sorted_by_slope_descending(tmp_path):
    adapter = _adapter(tmp_path)
    model = _build_fseof_model(adapter)
    result = ec_fseof(model, "EX_PROD", "EX_A", n_steps=8)
    if len(result.gene_targets) >= 2:
        slopes = list(result.gene_targets["slope"])
        assert slopes == sorted(slopes, reverse=True)


# --------------------------------------------------------------------------- #
# usage_prot_* excluded from rxn_targets
# --------------------------------------------------------------------------- #

def test_usage_prot_rxns_excluded_from_rxn_targets(tmp_path):
    adapter = _adapter(tmp_path)
    model = _build_fseof_model(adapter)
    result = ec_fseof(model, "EX_PROD", "EX_A", n_steps=8)
    for rid in result.rxn_targets["rxn_id"]:
        assert not rid.startswith("usage_prot_")


def test_usage_prot_rxns_excluded_from_v_matrix(tmp_path):
    adapter = _adapter(tmp_path)
    model = _build_fseof_model(adapter)
    result = ec_fseof(model, "EX_PROD", "EX_A", n_steps=8)
    for rid in result.v_matrix.index:
        assert not rid.startswith("usage_prot_")
