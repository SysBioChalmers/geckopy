"""Tests for load_pax_db."""
from pathlib import Path

import numpy as np
import pytest

from geckopy.databases import ProtData, UniprotDB, load_pax_db


# --------------------------------------------------------------------------- #
# Fixture builders
# --------------------------------------------------------------------------- #

def _uniprot(rows: list[tuple[str, str, float]]) -> UniprotDB:
    """Build a UniprotDB from (uniprot_id, gene_name, mw) rows."""
    return UniprotDB(
        ids=[r[0] for r in rows],
        genes=[r[1] for r in rows],
        eccodes=[""] * len(rows),
        mw=np.array([r[2] for r in rows], dtype=float),
        sequences=[""] * len(rows),
    )


def _write_pax_db(
    path: Path,
    rows: list[tuple[str, str, float]],
    *,
    n_header: int = 2,
) -> None:
    """Write a paxDB.tsv with `n_header` `#` lines and (internal_id,
    gene_id, level) data rows."""
    lines = [f"# header line {i}" for i in range(n_header)]
    for internal_id, gene_id, level in rows:
        lines.append(f"{internal_id}\t{gene_id}\t{level}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


# --------------------------------------------------------------------------- #
# Trivial cases
# --------------------------------------------------------------------------- #

def test_missing_file_raises(tmp_path):
    db = _uniprot([("P1", "g1", 100.0)])
    with pytest.raises(FileNotFoundError, match="paxDB file"):
        load_pax_db(tmp_path / "missing.tsv", db)


def test_empty_file_yields_empty_data(tmp_path):
    p = tmp_path / "empty.tsv"
    _write_pax_db(p, [], n_header=0)
    db = _uniprot([("P1", "g1", 100.0)])
    out = load_pax_db(p, db)
    assert isinstance(out, ProtData)
    assert out.uniprot_ids == []
    assert out.abundances.size == 0


def test_basic_single_row(tmp_path):
    p = tmp_path / "pax.tsv"
    _write_pax_db(p, [(1, "g1", 2.0)])
    db = _uniprot([("P1", "g1", 100.0)])
    out = load_pax_db(p, db)
    assert out.uniprot_ids == ["P1"]
    np.testing.assert_array_equal(out.abundances, [2.0 * 100.0])


# --------------------------------------------------------------------------- #
# Header skipping
# --------------------------------------------------------------------------- #

def test_arbitrary_header_count_skipped(tmp_path):
    p = tmp_path / "pax.tsv"
    _write_pax_db(p, [(1, "g1", 1.0)], n_header=10)
    db = _uniprot([("P1", "g1", 100.0)])
    out = load_pax_db(p, db)
    assert len(out.uniprot_ids) == 1


def test_no_header_lines_works(tmp_path):
    p = tmp_path / "pax.tsv"
    _write_pax_db(p, [(1, "g1", 1.0)], n_header=0)
    db = _uniprot([("P1", "g1", 100.0)])
    out = load_pax_db(p, db)
    assert out.uniprot_ids == ["P1"]


# --------------------------------------------------------------------------- #
# Gene ID stripping
# --------------------------------------------------------------------------- #

def test_strips_paxdb_internal_prefix(tmp_path):
    """paxDB gene IDs like '9606.ENSP00012345' should have '9606.' stripped."""
    p = tmp_path / "pax.tsv"
    _write_pax_db(p, [(1, "9606.ENSP00012345", 1.0)])
    db = _uniprot([("P1", "ENSP00012345", 100.0)])
    out = load_pax_db(p, db)
    assert out.uniprot_ids == ["P1"]


def test_no_prefix_unchanged(tmp_path):
    p = tmp_path / "pax.tsv"
    _write_pax_db(p, [(1, "myGene", 1.0)])
    db = _uniprot([("P1", "myGene", 100.0)])
    out = load_pax_db(p, db)
    assert out.uniprot_ids == ["P1"]


def test_only_leading_digits_stripped(tmp_path):
    """Internal periods like 'foo.bar.baz' are NOT stripped (only the
    leading `\\d+\\.` pattern)."""
    p = tmp_path / "pax.tsv"
    _write_pax_db(p, [(1, "foo.bar.baz", 1.0)])
    db = _uniprot([("P1", "foo.bar.baz", 100.0)])
    out = load_pax_db(p, db)
    assert out.uniprot_ids == ["P1"]


# --------------------------------------------------------------------------- #
# Unmatched genes / abundance computation
# --------------------------------------------------------------------------- #

def test_unmatched_gene_dropped_silently(tmp_path):
    p = tmp_path / "pax.tsv"
    _write_pax_db(p, [
        (1, "g1", 1.0),
        (2, "ghost", 5.0),  # no UniProt match
        (3, "g2", 3.0),
    ])
    db = _uniprot([("P1", "g1", 100.0), ("P2", "g2", 200.0)])
    out = load_pax_db(p, db)
    assert out.uniprot_ids == ["P1", "P2"]


def test_abundance_is_level_times_mw(tmp_path):
    p = tmp_path / "pax.tsv"
    _write_pax_db(p, [
        (1, "g1", 2.5),
        (2, "g2", 7.0),
    ])
    db = _uniprot([("P1", "g1", 50.0), ("P2", "g2", 200.0)])
    out = load_pax_db(p, db)
    np.testing.assert_array_equal(
        out.abundances, [2.5 * 50.0, 7.0 * 200.0],
    )


def test_first_match_wins_for_duplicate_gene_names(tmp_path):
    """If a gene name appears twice in UniprotDB, the first index wins."""
    p = tmp_path / "pax.tsv"
    _write_pax_db(p, [(1, "g1", 1.0)])
    db = _uniprot([
        ("P_first", "g1", 100.0),
        ("P_second", "g1", 999.0),
    ])
    out = load_pax_db(p, db)
    assert out.uniprot_ids == ["P_first"]
    np.testing.assert_array_equal(out.abundances, [100.0])


# --------------------------------------------------------------------------- #
# Malformed lines
# --------------------------------------------------------------------------- #

def test_short_line_skipped(tmp_path):
    p = tmp_path / "pax.tsv"
    p.write_text(
        "# hdr\n"
        "1\tg1\t1.0\n"
        "2\tg2\n"   # only 2 cols
        "3\tg3\t3.0\n",
        encoding="utf-8",
    )
    db = _uniprot([
        ("P1", "g1", 100.0),
        ("P2", "g2", 100.0),
        ("P3", "g3", 100.0),
    ])
    out = load_pax_db(p, db)
    assert out.uniprot_ids == ["P1", "P3"]


def test_non_numeric_level_skipped(tmp_path):
    p = tmp_path / "pax.tsv"
    p.write_text(
        "# hdr\n"
        "1\tg1\t1.0\n"
        "2\tg2\tNA\n"
        "3\tg3\t3.0\n",
        encoding="utf-8",
    )
    db = _uniprot([
        ("P1", "g1", 100.0),
        ("P2", "g2", 100.0),
        ("P3", "g3", 100.0),
    ])
    out = load_pax_db(p, db)
    assert out.uniprot_ids == ["P1", "P3"]
