"""Tests for constrain_flux_data."""
from pathlib import Path

import cobra
import numpy as np
import pytest
from scipy import sparse

from geckopy import EcModel, ModelAdapter
from geckopy.databases import FluxData
from geckopy.ec_model.ec_data import EcData
from geckopy.limit_proteins import constrain_flux_data


# --------------------------------------------------------------------------- #
# Fixture builders
# --------------------------------------------------------------------------- #

def _adapter(
    tmp_path: Path,
    *,
    bio_rxn: str = "biomass",
    c_source: str = "EX_glc",
) -> ModelAdapter:
    (tmp_path / "model_adapter.toml").write_text(
        f'conv_gem = "dummy.xml"\n'
        f'org_name = "test"\n'
        f'bio_rxn = "{bio_rxn}"\n'
        f'c_source = "{c_source}"\n'
    )
    return ModelAdapter.from_folder(tmp_path)


def _build_simple_model(adapter: ModelAdapter) -> EcModel:
    """Tiny model with biomass + a few exchange rxns."""
    model = EcModel("toy", adapter=adapter)

    glc = cobra.Metabolite("glc_e", compartment="e")
    O2 = cobra.Metabolite("O2_e", compartment="e")
    eth = cobra.Metabolite("eth_e", compartment="e")
    bio_met = cobra.Metabolite("bio", compartment="c")
    model.add_metabolites([glc, O2, eth, bio_met])

    EX_glc = cobra.Reaction("EX_glc")
    EX_glc.add_metabolites({glc: -1.0})
    EX_glc.lower_bound = -1000.0
    EX_glc.upper_bound = 1000.0

    EX_O2 = cobra.Reaction("EX_O2")
    EX_O2.add_metabolites({O2: -1.0})
    EX_O2.lower_bound = -1000.0
    EX_O2.upper_bound = 1000.0

    EX_eth = cobra.Reaction("EX_eth")
    EX_eth.add_metabolites({eth: -1.0})
    EX_eth.lower_bound = -1000.0
    EX_eth.upper_bound = 1000.0

    BIO = cobra.Reaction("biomass")
    BIO.add_metabolites({bio_met: 1.0})
    BIO.lower_bound = 0.0
    BIO.upper_bound = 1000.0

    model.add_reactions([EX_glc, EX_O2, EX_eth, BIO])

    model.ec = EcData(
        rxns=[],
        kcat=np.empty(0, dtype=float),
        source=[],
        notes=[],
        eccodes=[],
        rxn_enz_mat=sparse.csr_matrix((0, 0), dtype=float),
    )
    return model


def _flux_data(
    *,
    conds: list[str],
    p_tot: list[float],
    gr_rate: list[float],
    exch_fluxes: list[list[float]],
    exch_rxn_ids: list[str],
) -> FluxData:
    return FluxData(
        conds=list(conds),
        p_tot=np.array(p_tot, dtype=float),
        gr_rate=np.array(gr_rate, dtype=float),
        exch_fluxes=np.array(exch_fluxes, dtype=float),
        exch_mets=list(exch_rxn_ids),
        exch_rxn_ids=list(exch_rxn_ids),
    )


# --------------------------------------------------------------------------- #
# Pre-flight
# --------------------------------------------------------------------------- #

def test_no_adapter_raises(tmp_path):
    adapter = _adapter(tmp_path)
    model = _build_simple_model(adapter)
    model.adapter = None
    fd = _flux_data(
        conds=["c1"], p_tot=[0.5], gr_rate=[0.4],
        exch_fluxes=[[-5.0]], exch_rxn_ids=["EX_glc"],
    )
    with pytest.raises(ValueError, match="adapter"):
        constrain_flux_data(model, fd)


def test_invalid_max_min_growth_raises(tmp_path):
    adapter = _adapter(tmp_path)
    model = _build_simple_model(adapter)
    fd = _flux_data(
        conds=["c1"], p_tot=[0.5], gr_rate=[0.4],
        exch_fluxes=[[-5.0]], exch_rxn_ids=["EX_glc"],
    )
    with pytest.raises(ValueError, match="max_min_growth"):
        constrain_flux_data(model, fd, max_min_growth="invalid")


def test_invalid_loose_strict_raises(tmp_path):
    adapter = _adapter(tmp_path)
    model = _build_simple_model(adapter)
    fd = _flux_data(
        conds=["c1"], p_tot=[0.5], gr_rate=[0.4],
        exch_fluxes=[[-5.0]], exch_rxn_ids=["EX_glc"],
    )
    with pytest.raises(ValueError, match="loose_strict_flux"):
        constrain_flux_data(model, fd, loose_strict_flux=200)


def test_unknown_condition_name_raises(tmp_path):
    adapter = _adapter(tmp_path)
    model = _build_simple_model(adapter)
    fd = _flux_data(
        conds=["c1"], p_tot=[0.5], gr_rate=[0.4],
        exch_fluxes=[[-5.0]], exch_rxn_ids=["EX_glc"],
    )
    with pytest.raises(ValueError, match="Condition"):
        constrain_flux_data(model, fd, condition="nonexistent")


def test_out_of_range_condition_index_raises(tmp_path):
    adapter = _adapter(tmp_path)
    model = _build_simple_model(adapter)
    fd = _flux_data(
        conds=["c1"], p_tot=[0.5], gr_rate=[0.4],
        exch_fluxes=[[-5.0]], exch_rxn_ids=["EX_glc"],
    )
    with pytest.raises(IndexError, match="condition index"):
        constrain_flux_data(model, fd, condition=5)


def test_unknown_rxn_id_raises(tmp_path):
    adapter = _adapter(tmp_path)
    model = _build_simple_model(adapter)
    fd = _flux_data(
        conds=["c1"], p_tot=[0.5], gr_rate=[0.4],
        exch_fluxes=[[-5.0]], exch_rxn_ids=["EX_nonexistent"],
    )
    with pytest.raises(ValueError, match="not present in the model"):
        constrain_flux_data(model, fd)


# --------------------------------------------------------------------------- #
# Growth bounds
# --------------------------------------------------------------------------- #

def test_max_growth_sets_upper_bound(tmp_path):
    adapter = _adapter(tmp_path)
    model = _build_simple_model(adapter)
    fd = _flux_data(
        conds=["c1"], p_tot=[0.5], gr_rate=[0.4],
        exch_fluxes=[[-5.0]], exch_rxn_ids=["EX_glc"],
    )
    constrain_flux_data(model, fd, max_min_growth="max")
    bio = model.reactions.get_by_id("biomass")
    assert bio.lower_bound == 0.0
    assert bio.upper_bound == pytest.approx(0.4)


def test_min_growth_sets_lower_bound(tmp_path):
    adapter = _adapter(tmp_path)
    model = _build_simple_model(adapter)
    fd = _flux_data(
        conds=["c1"], p_tot=[0.5], gr_rate=[0.4],
        exch_fluxes=[[-5.0]], exch_rxn_ids=["EX_glc"],
    )
    constrain_flux_data(model, fd, max_min_growth="min")
    bio = model.reactions.get_by_id("biomass")
    assert bio.lower_bound == pytest.approx(0.4)
    assert bio.upper_bound == 1000.0


# --------------------------------------------------------------------------- #
# c_source reset
# --------------------------------------------------------------------------- #

def test_c_source_zeroed(tmp_path):
    adapter = _adapter(tmp_path, c_source="EX_eth")
    model = _build_simple_model(adapter)
    fd = _flux_data(
        conds=["c1"], p_tot=[0.5], gr_rate=[0.4],
        exch_fluxes=[[-5.0]], exch_rxn_ids=["EX_glc"],
    )
    constrain_flux_data(model, fd)
    eth_rxn = model.reactions.get_by_id("EX_eth")
    assert eth_rxn.lower_bound == 0.0
    assert eth_rxn.upper_bound == 0.0


def test_empty_c_source_skipped(tmp_path):
    adapter = _adapter(tmp_path, c_source="")
    model = _build_simple_model(adapter)
    fd = _flux_data(
        conds=["c1"], p_tot=[0.5], gr_rate=[0.4],
        exch_fluxes=[[-5.0]], exch_rxn_ids=["EX_glc"],
    )
    # Just shouldn't raise; bounds for non-data rxns left untouched.
    constrain_flux_data(model, fd)


# --------------------------------------------------------------------------- #
# Loose mode
# --------------------------------------------------------------------------- #

def test_loose_negative_flux_caps_lb(tmp_path):
    adapter = _adapter(tmp_path)
    model = _build_simple_model(adapter)
    fd = _flux_data(
        conds=["c1"], p_tot=[0.5], gr_rate=[0.4],
        exch_fluxes=[[-5.0]], exch_rxn_ids=["EX_glc"],
    )
    constrain_flux_data(model, fd, loose_strict_flux="loose")
    glc_rxn = model.reactions.get_by_id("EX_glc")
    assert glc_rxn.lower_bound == -5.0
    assert glc_rxn.upper_bound == 0.0


def test_loose_positive_flux_caps_ub(tmp_path):
    adapter = _adapter(tmp_path)
    model = _build_simple_model(adapter)
    fd = _flux_data(
        conds=["c1"], p_tot=[0.5], gr_rate=[0.4],
        exch_fluxes=[[3.0]], exch_rxn_ids=["EX_eth"],
    )
    constrain_flux_data(model, fd, loose_strict_flux="loose")
    eth_rxn = model.reactions.get_by_id("EX_eth")
    assert eth_rxn.lower_bound == 0.0
    assert eth_rxn.upper_bound == 3.0


def test_loose_nan_flux_skipped(tmp_path):
    adapter = _adapter(tmp_path)
    model = _build_simple_model(adapter)
    original_lb = model.reactions.get_by_id("EX_O2").lower_bound
    fd = _flux_data(
        conds=["c1"], p_tot=[0.5], gr_rate=[0.4],
        exch_fluxes=[[np.nan]], exch_rxn_ids=["EX_O2"],
    )
    constrain_flux_data(model, fd)
    # NaN means "no measurement"; bounds untouched.
    assert model.reactions.get_by_id("EX_O2").lower_bound == original_lb


# --------------------------------------------------------------------------- #
# Percentage variance mode
# --------------------------------------------------------------------------- #

def test_percentage_brackets_positive_value(tmp_path):
    """pct=10 -> lb = val*0.95, ub = val*1.05."""
    adapter = _adapter(tmp_path)
    model = _build_simple_model(adapter)
    fd = _flux_data(
        conds=["c1"], p_tot=[0.5], gr_rate=[0.4],
        exch_fluxes=[[3.0]], exch_rxn_ids=["EX_eth"],
    )
    constrain_flux_data(model, fd, loose_strict_flux=10)
    eth_rxn = model.reactions.get_by_id("EX_eth")
    assert eth_rxn.lower_bound == pytest.approx(2.85)
    assert eth_rxn.upper_bound == pytest.approx(3.15)


def test_percentage_handles_negative_value(tmp_path):
    """For negative val, lb and ub are swapped to keep lb <= ub."""
    adapter = _adapter(tmp_path)
    model = _build_simple_model(adapter)
    fd = _flux_data(
        conds=["c1"], p_tot=[0.5], gr_rate=[0.4],
        exch_fluxes=[[-5.0]], exch_rxn_ids=["EX_glc"],
    )
    constrain_flux_data(model, fd, loose_strict_flux=10)
    glc_rxn = model.reactions.get_by_id("EX_glc")
    # val=-5, pct=10 -> raw [-5*0.95, -5*1.05] = [-4.75, -5.25]
    # After swap: lb = -5.25, ub = -4.75.
    assert glc_rxn.lower_bound == pytest.approx(-5.25)
    assert glc_rxn.upper_bound == pytest.approx(-4.75)


# --------------------------------------------------------------------------- #
# Extreme +/-1000 sentinel
# --------------------------------------------------------------------------- #

def test_minus_thousand_is_unconstrained_uptake(tmp_path):
    adapter = _adapter(tmp_path)
    model = _build_simple_model(adapter)
    fd = _flux_data(
        conds=["c1"], p_tot=[0.5], gr_rate=[0.4],
        exch_fluxes=[[-1000.0]], exch_rxn_ids=["EX_O2"],
    )
    constrain_flux_data(model, fd, loose_strict_flux=10)
    o2_rxn = model.reactions.get_by_id("EX_O2")
    assert o2_rxn.lower_bound == -1000.0
    assert o2_rxn.upper_bound == 0.0


def test_plus_thousand_is_unconstrained_excretion(tmp_path):
    adapter = _adapter(tmp_path)
    model = _build_simple_model(adapter)
    fd = _flux_data(
        conds=["c1"], p_tot=[0.5], gr_rate=[0.4],
        exch_fluxes=[[1000.0]], exch_rxn_ids=["EX_eth"],
    )
    constrain_flux_data(model, fd, loose_strict_flux=10)
    eth_rxn = model.reactions.get_by_id("EX_eth")
    assert eth_rxn.lower_bound == 0.0
    assert eth_rxn.upper_bound == 1000.0


# --------------------------------------------------------------------------- #
# Condition selection
# --------------------------------------------------------------------------- #

def test_condition_by_index(tmp_path):
    adapter = _adapter(tmp_path)
    model = _build_simple_model(adapter)
    fd = _flux_data(
        conds=["c1", "c2"],
        p_tot=[0.5, 0.6],
        gr_rate=[0.4, 0.7],
        exch_fluxes=[[-5.0], [-8.0]],
        exch_rxn_ids=["EX_glc"],
    )
    constrain_flux_data(model, fd, condition=1)
    bio = model.reactions.get_by_id("biomass")
    assert bio.upper_bound == pytest.approx(0.7)
    glc_rxn = model.reactions.get_by_id("EX_glc")
    assert glc_rxn.lower_bound == -8.0


def test_condition_by_name(tmp_path):
    adapter = _adapter(tmp_path)
    model = _build_simple_model(adapter)
    fd = _flux_data(
        conds=["c1", "c2"],
        p_tot=[0.5, 0.6],
        gr_rate=[0.4, 0.7],
        exch_fluxes=[[-5.0], [-8.0]],
        exch_rxn_ids=["EX_glc"],
    )
    constrain_flux_data(model, fd, condition="c2")
    glc_rxn = model.reactions.get_by_id("EX_glc")
    assert glc_rxn.lower_bound == -8.0


# --------------------------------------------------------------------------- #
# Multi-rxn flux scenarios
# --------------------------------------------------------------------------- #

def test_multiple_rxns_constrained_correctly(tmp_path):
    adapter = _adapter(tmp_path)
    model = _build_simple_model(adapter)
    fd = _flux_data(
        conds=["c1"], p_tot=[0.5], gr_rate=[0.4],
        exch_fluxes=[[-5.0, -2.0, 3.0]],
        exch_rxn_ids=["EX_glc", "EX_O2", "EX_eth"],
    )
    constrain_flux_data(model, fd, loose_strict_flux="loose")
    glc = model.reactions.get_by_id("EX_glc")
    o2 = model.reactions.get_by_id("EX_O2")
    eth = model.reactions.get_by_id("EX_eth")
    assert (glc.lower_bound, glc.upper_bound) == (-5.0, 0.0)
    assert (o2.lower_bound, o2.upper_bound) == (-2.0, 0.0)
    assert (eth.lower_bound, eth.upper_bound) == (0.0, 3.0)
