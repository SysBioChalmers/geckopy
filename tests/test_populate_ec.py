"""Tests for stages 6, 7, 8 of make_ec_model.

These three stages populate the ec substructure. Since populate_ec
operates on the result of stages 1-5, the tests that touch real model
fixtures pipe through the earlier stages as well.
"""
from pathlib import Path

import cobra
import numpy as np
import pytest
from scipy import sparse

from geckopy import EcModel, ModelAdapter
from geckopy.databases import UniprotDB, load_uniprot_tsv
from geckopy.ec_model.pipeline import (
    allocate_ec_for_catalyzed_reactions,
    build_rxn_enzyme_coupling,
    convert_to_irreversible,
    expand_model,
    invert_backwards_only_reactions,
    populate_enzyme_data,
    remove_pseudoreaction_gprs,
)

EXAMPLE_DIR = Path(__file__).parents[1] / "examples" / "ecTestGEM"


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #

def _minimal_adapter(tmp_path: Path) -> ModelAdapter:
    (tmp_path / "model_adapter.toml").write_text(
        'conv_gem = "dummy.xml"\norg_name = "test"\n'
    )
    return ModelAdapter.from_folder(tmp_path)


def _build_model(
    reactions: list[tuple[str, dict[str, float], float, float, str]],
) -> cobra.Model:
    model = cobra.Model("test")
    mets: dict[str, cobra.Metabolite] = {}
    for _, stoich, _, _, _ in reactions:
        for met_id in stoich:
            if met_id not in mets:
                mets[met_id] = cobra.Metabolite(met_id, compartment="c")
    for rxn_id, stoich, lb, ub, gpr in reactions:
        rxn = cobra.Reaction(rxn_id)
        rxn.lower_bound = lb
        rxn.upper_bound = ub
        rxn.add_metabolites({mets[m]: c for m, c in stoich.items()})
        if gpr:
            rxn.gene_reaction_rule = gpr
        model.add_reactions([rxn])
    return model


def _synthetic_uniprot(rows: list[tuple[str, str, float, str]]) -> UniprotDB:
    """Build a UniprotDB from (id, gene, mw_kda, seq) tuples, no EC codes."""
    return UniprotDB(
        ids=[r[0] for r in rows],
        genes=[r[1] for r in rows],
        eccodes=[""] * len(rows),
        mw=np.array([r[2] for r in rows], dtype=float),
        sequences=[r[3] for r in rows],
    )


def _run_preprocess(model: cobra.Model) -> None:
    """Run stages 1, 2, 4, 5 in order. Stage 1 needs an adapter;
    we pass one with a path whose data/pseudoRxns.tsv doesn't exist."""
    # Using a throwaway adapter with no pseudoRxns.tsv, so stage 1 is a no-op
    # unless a reaction name literally contains "pseudoreaction".


# --------------------------------------------------------------------------- #
# Stage 6: allocate_ec_for_catalyzed_reactions
# --------------------------------------------------------------------------- #

def test_stage6_selects_only_reactions_with_genes():
    model = _build_model([
        ("r1", {"A": -1.0, "B": 1.0}, 0.0, 1000.0, "g1"),
        ("r2", {"B": -1.0, "C": 1.0}, 0.0, 1000.0, ""),       # no GPR
        ("r3", {"C": -1.0, "D": 1.0}, 0.0, 1000.0, "g2 and g3"),
    ])
    ec_model = EcModel.from_cobra(model, adapter=None)
    allocate_ec_for_catalyzed_reactions(ec_model)

    assert ec_model.ec.rxns == ["r1", "r3"]
    assert ec_model.ec.n_rxns == 2
    assert np.isnan(ec_model.ec.kcat).all()
    assert ec_model.ec.source == ["", ""]
    assert ec_model.ec.notes == ["", ""]
    assert ec_model.ec.eccodes == ["", ""]


def test_stage6_empty_model():
    ec_model = EcModel("empty")
    rxn_ids = allocate_ec_for_catalyzed_reactions(ec_model)
    assert rxn_ids == []
    assert ec_model.ec.n_rxns == 0


def test_stage6_no_gene_associated_reactions():
    model = _build_model([
        ("r1", {"A": -1.0, "B": 1.0}, 0.0, 1000.0, ""),
    ])
    ec_model = EcModel.from_cobra(model, adapter=None)
    rxn_ids = allocate_ec_for_catalyzed_reactions(ec_model)
    assert rxn_ids == []
    assert ec_model.ec.n_rxns == 0


def test_stage6_preserves_model_order():
    model = _build_model([
        ("r3", {"A": -1.0, "B": 1.0}, 0.0, 1000.0, "g1"),
        ("r1", {"B": -1.0, "C": 1.0}, 0.0, 1000.0, "g2"),
        ("r2", {"C": -1.0, "D": 1.0}, 0.0, 1000.0, "g3"),
    ])
    ec_model = EcModel.from_cobra(model, adapter=None)
    allocate_ec_for_catalyzed_reactions(ec_model)
    # Order matches model.reactions, not sorted.
    assert ec_model.ec.rxns == ["r3", "r1", "r2"]


# --------------------------------------------------------------------------- #
# Stage 7: populate_enzyme_data
# --------------------------------------------------------------------------- #

def test_stage7_populates_matched_genes(tmp_path):
    adapter = _minimal_adapter(tmp_path)
    model = _build_model([
        ("r1", {"A": -1.0, "B": 1.0}, 0.0, 1000.0, "g1 and g2"),
    ])
    ec_model = EcModel.from_cobra(model, adapter=adapter)
    allocate_ec_for_catalyzed_reactions(ec_model)

    db = _synthetic_uniprot([
        ("P1", "g1", 10.0, "MKAL"),
        ("P2", "g2", 20.0, "MNTD"),
    ])

    no_uniprot = populate_enzyme_data(ec_model, db)
    assert no_uniprot == []
    assert ec_model.ec.genes == ["g1", "g2"]   # alphabetically sorted
    assert ec_model.ec.enzymes == ["P1", "P2"]
    np.testing.assert_array_equal(ec_model.ec.mw, [10.0, 20.0])
    assert ec_model.ec.sequence == ["MKAL", "MNTD"]
    assert np.isnan(ec_model.ec.concs).all()


def test_stage7_alphabetical_gene_order(tmp_path):
    """Stage 7 sorts genes alphabetically regardless of model order."""
    adapter = _minimal_adapter(tmp_path)
    model = _build_model([
        ("r1", {"A": -1.0, "B": 1.0}, 0.0, 1000.0, "gZ and gA and gM"),
    ])
    ec_model = EcModel.from_cobra(model, adapter=adapter)
    allocate_ec_for_catalyzed_reactions(ec_model)

    db = _synthetic_uniprot([
        ("P1", "gZ", 10.0, "S1"),
        ("P2", "gA", 20.0, "S2"),
        ("P3", "gM", 30.0, "S3"),
    ])
    populate_enzyme_data(ec_model, db)

    assert ec_model.ec.genes == ["gA", "gM", "gZ"]
    assert ec_model.ec.enzymes == ["P2", "P3", "P1"]


def test_stage7_reports_unmatched_genes(tmp_path):
    adapter = _minimal_adapter(tmp_path)
    model = _build_model([
        ("r1", {"A": -1.0, "B": 1.0}, 0.0, 1000.0, "g1 and g2"),
        ("r2", {"B": -1.0, "C": 1.0}, 0.0, 1000.0, "g3"),
    ])
    ec_model = EcModel.from_cobra(model, adapter=adapter)
    allocate_ec_for_catalyzed_reactions(ec_model)

    # Only g1 is in UniProt. g2 and g3 are missing.
    db = _synthetic_uniprot([("P1", "g1", 10.0, "S1")])
    no_uniprot = populate_enzyme_data(ec_model, db)

    assert sorted(no_uniprot) == ["g2", "g3"]
    assert ec_model.ec.genes == ["g1"]
    assert ec_model.ec.enzymes == ["P1"]


def test_stage7_raises_when_nothing_matches(tmp_path):
    adapter = _minimal_adapter(tmp_path)
    model = _build_model([
        ("r1", {"A": -1.0, "B": 1.0}, 0.0, 1000.0, "g1"),
    ])
    ec_model = EcModel.from_cobra(model, adapter=adapter)
    allocate_ec_for_catalyzed_reactions(ec_model)

    db = _synthetic_uniprot([("P_OTHER", "g_other", 10.0, "S")])

    with pytest.raises(ValueError, match="None of the model genes"):
        populate_enzyme_data(ec_model, db)


def test_stage7_adds_warning_note_to_reactions_with_unmatched_genes(tmp_path):
    adapter = _minimal_adapter(tmp_path)
    model = _build_model([
        ("r1", {"A": -1.0, "B": 1.0}, 0.0, 1000.0, "g1 and g_missing"),
        ("r2", {"B": -1.0, "C": 1.0}, 0.0, 1000.0, "g1"),
    ])
    ec_model = EcModel.from_cobra(model, adapter=adapter)
    allocate_ec_for_catalyzed_reactions(ec_model)

    db = _synthetic_uniprot([("P1", "g1", 10.0, "S1")])
    populate_enzyme_data(ec_model, db)

    r1 = ec_model.reactions.get_by_id("r1")
    assert "geckopy_warning" in r1.notes
    assert "g_missing" in r1.notes["geckopy_warning"]

    r2 = ec_model.reactions.get_by_id("r2")
    assert "geckopy_warning" not in r2.notes


def test_stage7_warning_lists_all_missing_genes_for_one_reaction(tmp_path):
    adapter = _minimal_adapter(tmp_path)
    model = _build_model([
        ("r1", {"A": -1.0, "B": 1.0}, 0.0, 1000.0,
         "g1 and g_missingA and g_missingB"),
    ])
    ec_model = EcModel.from_cobra(model, adapter=adapter)
    allocate_ec_for_catalyzed_reactions(ec_model)

    db = _synthetic_uniprot([("P1", "g1", 10.0, "S1")])
    populate_enzyme_data(ec_model, db)

    note = ec_model.reactions.get_by_id("r1").notes["geckopy_warning"]
    assert "g_missingA" in note
    assert "g_missingB" in note


def test_stage7_uses_conversion_table_when_present(tmp_path):
    adapter = _minimal_adapter(tmp_path)
    (tmp_path / "data").mkdir()
    (tmp_path / "data" / "uniprotConversion.tsv").write_text(
        "gene\tuniprot_id\n"
        "gene_foo\tP_FOO\n"
        "gene_bar\tP_BAR\n"
    )
    # Rebuild adapter so it picks up the new file.
    adapter = ModelAdapter.from_folder(tmp_path)

    model = _build_model([
        ("r1", {"A": -1.0, "B": 1.0}, 0.0, 1000.0, "gene_foo"),
    ])
    ec_model = EcModel.from_cobra(model, adapter=adapter)
    allocate_ec_for_catalyzed_reactions(ec_model)

    # UniProt DB is keyed by Entry (P_FOO, P_BAR), genes don't match anything
    # literal; the conversion table routes gene_foo -> P_FOO.
    db = _synthetic_uniprot([
        ("P_FOO", "x", 10.0, "SF"),
        ("P_BAR", "y", 20.0, "SB"),
    ])

    populate_enzyme_data(ec_model, db)
    assert ec_model.ec.genes == ["gene_foo"]
    assert ec_model.ec.enzymes == ["P_FOO"]


def test_stage7_model_with_no_genes(tmp_path):
    adapter = _minimal_adapter(tmp_path)
    model = _build_model([
        ("r1", {"A": -1.0, "B": 1.0}, 0.0, 1000.0, ""),
    ])
    ec_model = EcModel.from_cobra(model, adapter=adapter)
    allocate_ec_for_catalyzed_reactions(ec_model)

    db = _synthetic_uniprot([])
    no_uniprot = populate_enzyme_data(ec_model, db)

    assert no_uniprot == []
    assert ec_model.ec.genes == []
    assert ec_model.ec.enzymes == []
    assert ec_model.ec.mw.shape == (0,)


def test_stage7_raises_without_adapter():
    model = _build_model([
        ("r1", {"A": -1.0, "B": 1.0}, 0.0, 1000.0, "g1"),
    ])
    ec_model = EcModel.from_cobra(model, adapter=None)
    allocate_ec_for_catalyzed_reactions(ec_model)

    db = _synthetic_uniprot([("P1", "g1", 10.0, "S")])
    with pytest.raises(ValueError, match="adapter is None"):
        populate_enzyme_data(ec_model, db)


# --------------------------------------------------------------------------- #
# Stage 8: build_rxn_enzyme_coupling
# --------------------------------------------------------------------------- #

def test_stage8_builds_correct_coupling_matrix(tmp_path):
    adapter = _minimal_adapter(tmp_path)
    model = _build_model([
        ("r1", {"A": -1.0, "B": 1.0}, 0.0, 1000.0, "g1 and g2"),
        ("r2", {"B": -1.0, "C": 1.0}, 0.0, 1000.0, "g3"),
    ])
    ec_model = EcModel.from_cobra(model, adapter=adapter)
    allocate_ec_for_catalyzed_reactions(ec_model)
    db = _synthetic_uniprot([
        ("P1", "g1", 10.0, "S1"),
        ("P2", "g2", 20.0, "S2"),
        ("P3", "g3", 30.0, "S3"),
    ])
    populate_enzyme_data(ec_model, db)

    build_rxn_enzyme_coupling(ec_model)

    mat = ec_model.ec.rxn_enz_mat
    assert isinstance(mat, sparse.csr_matrix)
    assert mat.shape == (2, 3)
    dense = mat.toarray()
    # ec.rxns order: [r1, r2]; ec.genes order: [g1, g2, g3]
    np.testing.assert_array_equal(
        dense,
        [[1.0, 1.0, 0.0],
         [0.0, 0.0, 1.0]],
    )


def test_stage8_all_nonzero_entries_are_1(tmp_path):
    adapter = _minimal_adapter(tmp_path)
    model = _build_model([
        ("r1", {"A": -1.0, "B": 1.0}, 0.0, 1000.0, "g1 and g2 and g3"),
    ])
    ec_model = EcModel.from_cobra(model, adapter=adapter)
    allocate_ec_for_catalyzed_reactions(ec_model)
    db = _synthetic_uniprot([
        ("P1", "g1", 10.0, "S1"),
        ("P2", "g2", 20.0, "S2"),
        ("P3", "g3", 30.0, "S3"),
    ])
    populate_enzyme_data(ec_model, db)
    build_rxn_enzyme_coupling(ec_model)

    # All nonzero entries should be exactly 1.0.
    assert set(ec_model.ec.rxn_enz_mat.data) == {1.0}


def test_stage8_empty_ec_produces_empty_matrix(tmp_path):
    adapter = _minimal_adapter(tmp_path)
    ec_model = EcModel("empty", adapter=adapter)
    build_rxn_enzyme_coupling(ec_model)
    assert ec_model.ec.rxn_enz_mat.shape == (0, 0)


def test_stage8_unmatched_gene_produces_zero_row(tmp_path):
    """A reaction whose only gene is unmatched should have an all-zero row."""
    adapter = _minimal_adapter(tmp_path)
    model = _build_model([
        ("r1", {"A": -1.0, "B": 1.0}, 0.0, 1000.0, "g1"),
        ("r2", {"B": -1.0, "C": 1.0}, 0.0, 1000.0, "g_missing"),
    ])
    ec_model = EcModel.from_cobra(model, adapter=adapter)
    allocate_ec_for_catalyzed_reactions(ec_model)
    db = _synthetic_uniprot([("P1", "g1", 10.0, "S1")])
    populate_enzyme_data(ec_model, db)
    build_rxn_enzyme_coupling(ec_model)

    dense = ec_model.ec.rxn_enz_mat.toarray()
    assert dense.shape == (2, 1)
    # r1 row has a 1 for g1, r2 row is all zero.
    np.testing.assert_array_equal(dense, [[1.0], [0.0]])


def test_stage8_partial_match_preserves_matched_gene(tmp_path):
    """A reaction with one matched and one unmatched gene should still
    have a 1 for the matched one."""
    adapter = _minimal_adapter(tmp_path)
    model = _build_model([
        ("r1", {"A": -1.0, "B": 1.0}, 0.0, 1000.0, "g1 and g_missing"),
    ])
    ec_model = EcModel.from_cobra(model, adapter=adapter)
    allocate_ec_for_catalyzed_reactions(ec_model)
    db = _synthetic_uniprot([("P1", "g1", 10.0, "S1")])
    populate_enzyme_data(ec_model, db)
    build_rxn_enzyme_coupling(ec_model)

    dense = ec_model.ec.rxn_enz_mat.toarray()
    assert dense.shape == (1, 1)
    np.testing.assert_array_equal(dense, [[1.0]])


# --------------------------------------------------------------------------- #
# End-to-end stages 1-8 on the real ecTestGEM fixture
# --------------------------------------------------------------------------- #

def _run_all_stages_on_ectestgem() -> EcModel:
    """Run stages 1-8 in order on the ecTestGEM fixture."""
    adapter = ModelAdapter.from_folder(EXAMPLE_DIR)
    cobra_model = cobra.io.read_sbml_model(str(adapter.params.conv_gem))

    # Stages 1-5 operate on plain cobra.Model.
    remove_pseudoreaction_gprs(cobra_model, adapter)
    invert_backwards_only_reactions(cobra_model)
    convert_to_irreversible(cobra_model)
    expand_model(cobra_model)

    # Promote to EcModel before stages 6-8.
    ec_model = EcModel.from_cobra(cobra_model, adapter=adapter)
    allocate_ec_for_catalyzed_reactions(ec_model)
    uniprot_db = load_uniprot_tsv(EXAMPLE_DIR / "data" / "uniprot.tsv")
    populate_enzyme_data(ec_model, uniprot_db)
    build_rxn_enzyme_coupling(ec_model)
    return ec_model


def test_ectestgem_stage6_expected_reactions():
    """Expected catalyzed reactions after stages 1-5 on the fixture:
    R2_EXP_1, R2_EXP_2, R2_REV_EXP_1, R2_REV_EXP_2, R3, R5."""
    ec_model = _run_all_stages_on_ectestgem()
    assert sorted(ec_model.ec.rxns) == sorted([
        "R2_EXP_1", "R2_EXP_2",
        "R2_REV_EXP_1", "R2_REV_EXP_2",
        "R3", "R5",
    ])


def test_ectestgem_stage7_expected_genes_and_enzymes():
    """All 5 genes (G1-G5) match UniProt (P1-P5 respectively)."""
    ec_model = _run_all_stages_on_ectestgem()
    assert ec_model.ec.genes == ["G1", "G2", "G3", "G4", "G5"]
    assert ec_model.ec.enzymes == ["P1", "P2", "P3", "P4", "P5"]
    np.testing.assert_array_almost_equal(
        ec_model.ec.mw, [10000.0, 20000.0, 30000.0, 40000.0, 50000.0]
    )
    assert ec_model.ec.sequence == ["MRAL", "MNTD", "MSYN", "MDFM", "MLFK"]
    assert np.isnan(ec_model.ec.concs).all()


def test_ectestgem_stage7_no_missing_genes():
    ec_model = _run_all_stages_on_ectestgem()
    assert ec_model.ec.n_enzymes == 5
    # No reaction should have the warning note set.
    for rxn in ec_model.reactions:
        assert "geckopy_warning" not in rxn.notes


def test_ectestgem_stage8_coupling_matrix_shape_and_entries():
    ec_model = _run_all_stages_on_ectestgem()
    mat = ec_model.ec.rxn_enz_mat.toarray()
    assert mat.shape == (6, 5)
    # All entries must be 0 or 1.
    assert set(np.unique(mat)).issubset({0.0, 1.0})
    # Expected 8 nonzero entries across all 6 rows:
    # R2_EXP_1 (G1, G2), R2_EXP_2 (G3),
    # R2_REV_EXP_1 (G1, G2), R2_REV_EXP_2 (G3),
    # R3 (G4), R5 (G5) = 2 + 1 + 2 + 1 + 1 + 1 = 8.
    assert int(mat.sum()) == 8


def test_ectestgem_stage8_row_for_specific_reaction():
    """R3 should have exactly one 1, at the column for G4."""
    ec_model = _run_all_stages_on_ectestgem()
    row_idx = ec_model.ec.rxns.index("R3")
    col_idx = ec_model.ec.genes.index("G4")
    mat = ec_model.ec.rxn_enz_mat.toarray()
    assert mat[row_idx, col_idx] == 1.0
    assert mat[row_idx, :].sum() == 1.0


def test_ectestgem_validate_ec_passes():
    """The constructed ecModel must pass the internal consistency check."""
    ec_model = _run_all_stages_on_ectestgem()
    ec_model.validate_ec()  # should not raise
