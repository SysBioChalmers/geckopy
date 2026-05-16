"""Tests for SBML I/O (MW_KCAT encoding)."""
from pathlib import Path

import cobra
import libsbml
import numpy as np
import pytest

from geckopy import EcModel, ModelAdapter, make_ec_model
from geckopy.ec_model.pipeline import apply_kcat_constraints
from geckopy.io.sbml import read_sbml_ec_model, write_sbml_ec_model
from geckopy.utilities import load_ec_model, save_ec_model

EXAMPLE_DIR = Path(__file__).parents[1] / "examples" / "ecTestGEM"
DATA_DIR = Path(__file__).parent / "data"
LEGACY_EC_COLI_CORE = DATA_DIR / "ec_coli_core.xml"


def _ectestgem_ec_model_with_kcats() -> EcModel:
    adapter = ModelAdapter.from_folder(EXAMPLE_DIR)
    cobra_model = cobra.io.read_sbml_model(str(adapter.params.conv_gem))
    model = make_ec_model(cobra_model, adapter)
    model.ec.kcat[:] = 10.0
    model.ec.concs[:] = np.nan
    apply_kcat_constraints(model)
    model.objective = "R3"
    return model, adapter


# --------------------------------------------------------------------------- #
# Round-trip
# --------------------------------------------------------------------------- #

def test_write_read_roundtrip_enzymes(tmp_path):
    model, adapter = _ectestgem_ec_model_with_kcats()
    out = tmp_path / "rt.xml"
    write_sbml_ec_model(model, out)
    reloaded = read_sbml_ec_model(out, adapter=adapter)
    # Enzyme set preserved.
    assert sorted(reloaded.ec.enzymes) == sorted(model.ec.enzymes)


def test_write_read_roundtrip_mw(tmp_path):
    model, adapter = _ectestgem_ec_model_with_kcats()
    out = tmp_path / "rt.xml"
    write_sbml_ec_model(model, out)
    reloaded = read_sbml_ec_model(out, adapter=adapter)
    for u in model.ec.enzymes:
        i_old = model.ec.enzymes.index(u)
        i_new = reloaded.ec.enzymes.index(u)
        assert reloaded.ec.mw[i_new] == pytest.approx(
            float(model.ec.mw[i_old]), abs=1e-9,
        )


def test_write_read_roundtrip_concs(tmp_path):
    model, adapter = _ectestgem_ec_model_with_kcats()
    # Set a couple of measured concs.
    model.ec.concs[0] = 0.05
    model.ec.concs[1] = 0.01
    out = tmp_path / "rt.xml"
    write_sbml_ec_model(model, out)
    reloaded = read_sbml_ec_model(out, adapter=adapter)
    u0, u1 = model.ec.enzymes[0], model.ec.enzymes[1]
    assert reloaded.ec.concs[reloaded.ec.enzymes.index(u0)] == pytest.approx(0.05)
    assert reloaded.ec.concs[reloaded.ec.enzymes.index(u1)] == pytest.approx(0.01)


def test_write_read_roundtrip_kcats(tmp_path):
    model, adapter = _ectestgem_ec_model_with_kcats()
    # Vary kcats so we can detect any rounding bugs.
    rng = np.random.default_rng(seed=42)
    for i in range(len(model.ec.kcat)):
        model.ec.kcat[i] = float(rng.uniform(1, 1000))
    apply_kcat_constraints(model)

    out = tmp_path / "rt.xml"
    write_sbml_ec_model(model, out)
    reloaded = read_sbml_ec_model(out, adapter=adapter)
    # For each (rxn, enzyme) in original, find the recovered kcat.
    for ri, rxn_id in enumerate(model.ec.rxns):
        kcat_orig = float(model.ec.kcat[ri])
        if not np.isfinite(kcat_orig):
            continue
        if rxn_id not in reloaded.ec.rxns:
            continue
        ri_new = reloaded.ec.rxns.index(rxn_id)
        assert reloaded.ec.kcat[ri_new] == pytest.approx(kcat_orig, rel=1e-5)


# --------------------------------------------------------------------------- #
# Group + notes structure
# --------------------------------------------------------------------------- #

def test_protein_group_present(tmp_path):
    model, _ = _ectestgem_ec_model_with_kcats()
    out = tmp_path / "rt.xml"
    write_sbml_ec_model(model, out)

    doc = libsbml.SBMLReader().readSBML(str(out))
    sbml_model = doc.getModel()
    groups_plugin = sbml_model.getPlugin("groups")
    assert groups_plugin is not None
    group_names = []
    for i in range(groups_plugin.getNumGroups()):
        group_names.append(groups_plugin.getGroup(i).getId())
    assert "Protein" in group_names

    # Verify members reference each prot_<id> species.
    protein_group = None
    for i in range(groups_plugin.getNumGroups()):
        g = groups_plugin.getGroup(i)
        if g.getId() == "Protein":
            protein_group = g
            break
    member_ids = {
        protein_group.getMember(i).getIdRef()
        for i in range(protein_group.getNumMembers())
    }
    for u in model.ec.enzymes:
        assert f"M_prot_{u}" in member_ids


# --------------------------------------------------------------------------- #
# Dispatch via load_ec_model / save_ec_model
# --------------------------------------------------------------------------- #

def test_load_dispatches_sbml(tmp_path):
    model, adapter = _ectestgem_ec_model_with_kcats()
    out = tmp_path / "by_load.xml"
    write_sbml_ec_model(model, out)
    reloaded = load_ec_model(out, adapter=adapter)
    assert isinstance(reloaded, EcModel)
    assert len(reloaded.ec.enzymes) == len(model.ec.enzymes)


def test_save_dispatches_sbml(tmp_path):
    model, adapter = _ectestgem_ec_model_with_kcats()
    out = tmp_path / "by_save.xml"
    written = save_ec_model(model, out, adapter=adapter)
    assert written == out
    assert out.is_file()
    # Sanity: file is a valid SBML and contains the Protein group.
    doc = libsbml.SBMLReader().readSBML(str(out))
    sbml_model = doc.getModel()
    assert sbml_model is not None
    groups_plugin = sbml_model.getPlugin("groups")
    group_names = [
        groups_plugin.getGroup(i).getId()
        for i in range(groups_plugin.getNumGroups())
    ]
    assert "Protein" in group_names


def test_legacy_geckopy_old_file_reads(tmp_path):
    """The fixture from geckopy_old must load without error and
    produce a sensible cobra model. The MW_KCAT encoding may not
    exactly match (legacy file uses the old conventions), but the
    cobra portion should always come through."""
    if not LEGACY_EC_COLI_CORE.is_file():
        pytest.skip("ec_coli_core.xml fixture missing")
    # We don't have a real adapter for E. coli; build a stub.
    (tmp_path / "model_adapter.toml").write_text(
        f'conv_gem = "{LEGACY_EC_COLI_CORE}"\n'
        'org_name = "Escherichia coli"\n'
    )
    adapter = ModelAdapter.from_folder(tmp_path)
    reloaded = read_sbml_ec_model(LEGACY_EC_COLI_CORE, adapter=adapter)
    assert isinstance(reloaded, EcModel)
    # Should at least produce a cobra model with some reactions.
    assert len(reloaded.reactions) > 50
