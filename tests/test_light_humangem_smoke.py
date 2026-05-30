"""Smoke test for a real-scale gecko-light build on Human-GEM.

Mirrors the MATLAB tutorial ``tutorials/light_ecModel/protocol.m``:
loads the unmodified Human-GEM YAML, builds a light ecModel, and
checks the output has the expected shape (thousands of ec rows with
the ``###_`` prefix, only ``prot_pool`` as the protein metabolite,
sane coupling matrix). The test is skipped when the Human-GEM repo
isn't checked out next to geckopy.

Run explicitly with:

    pytest tests/test_light_humangem_smoke.py -m smoke -q
"""
from __future__ import annotations

import csv
import os
from pathlib import Path

import cobra
import numpy as np
import pytest

from geckopy import ModelAdapter, make_ec_model
from geckopy.databases import UniprotDB

# --------------------------------------------------------------------------- #
# Locate Human-GEM
# --------------------------------------------------------------------------- #
#
# The smoke test needs two files from the Human-GEM repo
# (https://github.com/SysBioChalmers/Human-GEM):
#
#   model/Human-GEM.yml   — the cobra model
#   model/genes.tsv       — the ENSG → UniProt mapping
#
# Override either with an env var if your checkout lives elsewhere.

_DEFAULT_REPO = Path.home() / "github" / "Human-GEM"
_HUMANGEM_YML = Path(
    os.environ.get("HUMANGEM_YML", _DEFAULT_REPO / "model" / "Human-GEM.yml")
)
_HUMANGEM_GENES_TSV = Path(
    os.environ.get("HUMANGEM_GENES_TSV", _DEFAULT_REPO / "model" / "genes.tsv")
)

pytestmark = [
    pytest.mark.smoke,
    pytest.mark.skipif(
        not _HUMANGEM_YML.is_file() or not _HUMANGEM_GENES_TSV.is_file(),
        reason=(
            f"Human-GEM not found at {_HUMANGEM_YML} / {_HUMANGEM_GENES_TSV}. "
            "Set HUMANGEM_YML / HUMANGEM_GENES_TSV to override."
        ),
    ),
]


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #

def _build_humangem_uniprot_db() -> UniprotDB:
    """Synthesise a UniprotDB from the ENSG / UniProt columns of
    Human-GEM's genes.tsv. MW is set to a synthetic 50 kDa so the
    kcat math doesn't NaN-out downstream; this build is for shape
    smoke-testing only, not for quantitative kcat-driven simulation.
    """
    ids: list[str] = []
    genes: list[str] = []
    with open(_HUMANGEM_GENES_TSV, "r", encoding="utf-8") as fh:
        reader = csv.reader(fh, delimiter="\t", quotechar='"')
        header = next(reader)
        try:
            gene_col = header.index("genes")
            uniprot_col = header.index("geneUniProtID")
        except ValueError as e:  # column missing — schema drift
            raise RuntimeError(
                f"Unexpected genes.tsv schema: {header}"
            ) from e
        seen: set[str] = set()
        for row in reader:
            if len(row) <= max(gene_col, uniprot_col):
                continue
            ensg = row[gene_col].strip()
            uniprot = row[uniprot_col].strip()
            if not ensg or not uniprot or uniprot in seen:
                continue
            seen.add(uniprot)
            ids.append(uniprot)
            genes.append(ensg)
    n = len(ids)
    return UniprotDB(
        ids=ids,
        genes=genes,
        eccodes=[""] * n,
        mw=np.full(n, 50_000.0, dtype=float),
        sequences=["M"] * n,
    )


def _humangem_adapter(tmp_path: Path) -> ModelAdapter:
    """Write a minimal model_adapter.toml in tmp_path and load it.

    The toml's ``conv_gem`` points at the absolute Human-GEM.yml so the
    adapter's ``path`` (used to find optional data/ files like
    pseudoRxns.tsv) stays under tmp_path — keeps the test from
    polluting the Human-GEM checkout.
    """
    toml = (
        f'conv_gem = "{_HUMANGEM_YML.as_posix()}"\n'
        'org_name = "homo sapiens"\n'
        'sigma = 0.1\n'
        'p_tot = 0.5057\n'
        'f = 0.412\n'
        'c_source = "MAR09034"\n'
        'bio_rxn = "MAR13082"\n'
        'enzyme_comp = "Cytosol"\n'
        '\n'
        '[uniprot]\n'
        'type = "proteome"\n'
        'id = "UP000005640"\n'
        'gene_id_field = "gene_primary"\n'
        'reviewed = true\n'
        '\n'
        '[kegg]\n'
        'id = "hsa"\n'
        'gene_id = "Ensembl"\n'
        '\n'
        '[complex]\n'
        'taxonomic_id = 9606\n'
    )
    (tmp_path / "model_adapter.toml").write_text(toml)
    return ModelAdapter.from_folder(tmp_path)


# --------------------------------------------------------------------------- #
# The smoke test
# --------------------------------------------------------------------------- #

def test_light_humangem_builds_and_has_expected_shape(tmp_path):
    """Build a light ecModel from the unmodified Human-GEM YAML.

    Asserts:
    - the build completes;
    - ec.rxns has thousands of entries, all with the ``###_`` prefix;
    - ``prot_pool`` is the only protein metabolite;
    - no ``usage_prot_*`` reactions exist;
    - the coupling matrix has the right shape and at least one non-zero
      row per ec entry that survived gene matching.
    """
    cobra_model = cobra.io.load_yaml_model(str(_HUMANGEM_YML))
    assert len(cobra_model.reactions) > 5_000, (
        f"Human-GEM should have >5k reactions; got {len(cobra_model.reactions)}"
    )

    adapter = _humangem_adapter(tmp_path)
    uniprot_db = _build_humangem_uniprot_db()
    assert len(uniprot_db) > 1_000, (
        f"Synthetic UniprotDB should have >1k entries; got {len(uniprot_db)}"
    )

    ec_model = make_ec_model(
        cobra_model, adapter, gecko_light=True, uniprot_db=uniprot_db,
    )

    # Light flag set.
    assert ec_model.ec.gecko_light is True

    # ec.rxns: prefixed, thousands of entries.
    assert ec_model.ec.n_rxns > 1_000, (
        f"Expected >1k ec rows; got {ec_model.ec.n_rxns}"
    )
    bad_prefixes = [
        r for r in ec_model.ec.rxns[:200]
        if not (len(r) > 4 and r[3] == "_" and r[:3].isdigit())
    ]
    assert not bad_prefixes, (
        f"First 200 ec.rxns must all carry ``###_`` prefixes; "
        f"violators: {bad_prefixes[:5]}"
    )

    # Only the shared protein pool — no per-enzyme prot_<id> mets, no
    # usage_prot_<id> reactions.
    met_ids = {m.id for m in ec_model.metabolites}
    assert "prot_pool" in met_ids
    per_enz_mets = [
        m for m in met_ids if m.startswith("prot_") and m != "prot_pool"
    ]
    assert not per_enz_mets, (
        f"Light layout must skip per-enzyme prot mets; found "
        f"{len(per_enz_mets)} (e.g. {per_enz_mets[:3]})"
    )
    rxn_ids = {r.id for r in ec_model.reactions}
    assert "prot_pool_exchange" in rxn_ids
    usage_rxns = [r for r in rxn_ids if r.startswith("usage_prot_")]
    assert not usage_rxns, (
        f"Light layout must skip usage_prot_<id> reactions; found "
        f"{len(usage_rxns)}"
    )

    # Coupling matrix shape + density sanity.
    mat = ec_model.ec.rxn_enz_mat
    assert mat.shape == (ec_model.ec.n_rxns, ec_model.ec.n_enzymes)
    # Some rows may be all-zero (genes the synthetic UniprotDB
    # didn't cover), but the average over rows should be at least 1
    # (each gene-bearing isozyme contributes at least one enzyme).
    assert mat.nnz > ec_model.ec.n_rxns // 4, (
        f"Coupling matrix nnz={mat.nnz} suspiciously low for "
        f"{ec_model.ec.n_rxns} rows × {ec_model.ec.n_enzymes} enzymes"
    )

    # Self-consistency check.
    ec_model.ec.validate()
