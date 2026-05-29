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


_ECTESTGEM_CACHE: tuple[EcModel, ModelAdapter] | None = None


def _ectestgem_ec_model_with_kcats() -> tuple[EcModel, ModelAdapter]:
    """Cached build of the ecTestGEM ecModel with kcat=10 + R3 objective;
    deep-copied per call (only the model -- the adapter is shared)."""
    import copy as _copy
    global _ECTESTGEM_CACHE
    if _ECTESTGEM_CACHE is None:
        adapter = ModelAdapter.from_folder(EXAMPLE_DIR)
        cobra_model = cobra.io.read_sbml_model(str(adapter.params.conv_gem))
        model = make_ec_model(cobra_model, adapter)
        model.ec.kcat[:] = 10.0
        model.ec.concs[:] = np.nan
        apply_kcat_constraints(model)
        model.objective = "R3"
        _ECTESTGEM_CACHE = (model, adapter)
    cached_model, adapter = _ECTESTGEM_CACHE
    return _copy.deepcopy(cached_model), adapter


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


def test_read_without_adapter_for_inspection(tmp_path):
    """A model can be read for inspection without supplying an adapter."""
    model, _ = _ectestgem_ec_model_with_kcats()
    out = tmp_path / "rt.xml"
    write_sbml_ec_model(model, out)
    reloaded = read_sbml_ec_model(out, adapter=None)
    assert reloaded.adapter is None
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


def test_write_read_roundtrip_provenance(tmp_path):
    """source / eccodes / notes / sequence and >1 subunit counts survive
    the round-trip (previously silently dropped)."""
    model, adapter = _ectestgem_ec_model_with_kcats()

    # Pick a catalysed reaction with at least one coupled enzyme.
    mat = model.ec.rxn_enz_mat.tolil()
    ri = next(i for i in range(mat.shape[0]) if mat.getrow(i).nnz)
    rxn_id = model.ec.rxns[ri]
    j = int(mat.getrow(ri).indices[0])
    enzyme_id = model.ec.enzymes[j]

    # Stamp provenance and a 2-subunit coupling.
    model.ec.source[ri] = "brenda"
    model.ec.eccodes[ri] = "1.1.1.1;2.7.1.-"
    model.ec.notes[ri] = "preTuneKcat=12.5 | source:dlkcat"
    model.ec.kcat[ri] = 50.0
    model.ec.sequence[j] = "MKVLA"
    mat[ri, j] = 2.0
    model.ec.rxn_enz_mat = mat.tocsr()

    out = tmp_path / "rt.xml"
    write_sbml_ec_model(model, out)
    reloaded = read_sbml_ec_model(out, adapter=adapter)

    ri_new = reloaded.ec.rxns.index(rxn_id)
    j_new = reloaded.ec.enzymes.index(enzyme_id)
    assert reloaded.ec.source[ri_new] == "brenda"
    assert reloaded.ec.eccodes[ri_new] == "1.1.1.1;2.7.1.-"
    assert reloaded.ec.notes[ri_new] == "preTuneKcat=12.5 | source:dlkcat"
    assert reloaded.ec.sequence[j_new] == "MKVLA"
    # kcat recovered exactly, not divided by the subunit count.
    assert reloaded.ec.kcat[ri_new] == pytest.approx(50.0, rel=1e-9)
    assert reloaded.ec.rxn_enz_mat[ri_new, j_new] == pytest.approx(2.0)


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


def test_legacy_carrasco_file_reads(tmp_path):
    """An SBML written by the legacy geckopy (Carrasco et al., 2023,
    https://doi.org/10.1128/spectrum.01705-23) must load without
    error and produce a sensible cobra model. The MW_KCAT encoding
    may not exactly match (legacy file uses the old conventions),
    but the cobra portion should always come through."""
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
