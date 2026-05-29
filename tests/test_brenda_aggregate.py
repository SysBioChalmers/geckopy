"""Unit tests for BRENDA aggregation + TSV writing.

The snapshot ships one wide row per (ec, substrate, organism) triple
for kcat and SA (columns ``kcat_max`` / ``kcat_median`` carry the two
statistics computed from the raw measurements that fell into the
triple) and one row per (ec, organism) for MW.
"""
from __future__ import annotations

import hashlib
from pathlib import Path

import pytest

from geckopy.databases.brenda import Row, parse_brenda_json
from geckopy.databases.brenda.aggregate import aggregate_and_write

FIXTURE = Path(__file__).parent / "data" / "brenda_minimal.json"


def _read_tsv(path: Path) -> tuple[str, list[str], list[list[str]]]:
    """Return (release-comment line, column-header tokens, data rows)."""
    text = path.read_text(encoding="utf-8")
    lines = [ln for ln in text.split("\n") if ln]
    comment = lines[0]
    header = lines[1].split("\t")
    data = [ln.split("\t") for ln in lines[2:]]
    return comment, header, data


@pytest.fixture
def written(tmp_path):
    rows = list(parse_brenda_json(FIXTURE))
    paths = aggregate_and_write(
        rows, tmp_path, release="2026.1", date="2026-05-18",
    )
    return paths


def test_three_files_written(written):
    assert set(written.keys()) == {"kcat", "sa", "mw"}
    for p in written.values():
        assert p.exists()


def test_filenames_are_neutral(written):
    assert written["kcat"].name == "kcat.tsv"
    assert written["sa"].name == "sa.tsv"
    assert written["mw"].name == "mw.tsv"


def test_comment_line_format(written):
    comment, _, _ = _read_tsv(written["kcat"])
    assert comment.startswith("# BRENDA release 2026.1 generated 2026-05-18 - CC BY 4.0")
    assert "kcat in 1/s" in comment


def test_kcat_column_header(written):
    _, header, _ = _read_tsv(written["kcat"])
    assert header == [
        "ec_code", "substrate", "organism",
        "kcat_max", "kcat_median", "n", "references",
    ]


def test_sa_column_header(written):
    _, header, _ = _read_tsv(written["sa"])
    assert header == [
        "ec_code", "substrate", "organism",
        "sa_max", "sa_median", "n", "references",
    ]


def test_mw_column_header(written):
    _, header, _ = _read_tsv(written["mw"])
    assert header == [
        "ec_code", "substrate", "organism", "mw", "n", "references",
    ]


def test_kcat_triple_carries_both_max_and_median(tmp_path):
    """3 raw measurements share a triple → one wide row with both
    statistics and n=3."""
    rows = [
        Row("kcat", "1.1.1.1", "ethanol", "yeast", 10.0, ("PMID:1",)),
        Row("kcat", "1.1.1.1", "ethanol", "yeast", 50.0, ("PMID:2",)),
        Row("kcat", "1.1.1.1", "ethanol", "yeast", 30.0, ("PMID:3",)),
    ]
    paths = aggregate_and_write(rows, tmp_path, release="X", date="2026-01-01")
    _, _, data = _read_tsv(paths["kcat"])
    assert len(data) == 1
    row = data[0]
    assert row[0] == "1.1.1.1"
    assert row[3] == "50.0"        # kcat_max
    assert row[4] == "30.0"        # kcat_median
    assert row[5] == "3"           # n
    assert row[6] == "PMID:1;PMID:2;PMID:3"


def test_kcat_single_measurement_max_equals_median(tmp_path):
    """n=1 → max and median both equal the single value."""
    rows = [Row("kcat", "1.1.1.1", "ethanol", "yeast", 7.0, ("PMID:1",))]
    paths = aggregate_and_write(rows, tmp_path, release="X", date="2026-01-01")
    _, _, data = _read_tsv(paths["kcat"])
    assert len(data) == 1
    assert data[0][3] == "7.0"
    assert data[0][4] == "7.0"
    assert data[0][5] == "1"


def test_mw_single_row_no_aggregation_columns(tmp_path):
    """MW is not aggregation-split: one row per (ec, organism), max
    only."""
    rows = [
        Row("mw", "1.1.1.1", "*", "yeast", 50000.0, ("PMID:1",)),
        Row("mw", "1.1.1.1", "*", "yeast", 50100.0, ("PMID:2",)),
    ]
    paths = aggregate_and_write(rows, tmp_path, release="X", date="2026-01-01")
    _, header, data = _read_tsv(paths["mw"])
    assert "kcat_max" not in header
    assert "kcat_median" not in header
    assert len(data) == 1
    assert data[0][3] == "50100.0"  # max
    assert data[0][4] == "2"


def test_empty_references_emit_star(tmp_path):
    rows = [Row("kcat", "1.1.1.1", "ethanol", "yeast", 10.0, ())]
    paths = aggregate_and_write(rows, tmp_path, release="X", date="2026-01-01")
    _, _, data = _read_tsv(paths["kcat"])
    assert data[0][6] == "*"


def test_deterministic_sort(tmp_path):
    """Triples sort by (ec, substrate, organism); one row per triple."""
    rows = [
        Row("kcat", "2.7.1.1", "glucose",  "ecoli", 1.0, ()),
        Row("kcat", "1.1.1.1", "methanol", "yeast", 2.0, ()),
        Row("kcat", "1.1.1.1", "ethanol",  "yeast", 3.0, ()),
        Row("kcat", "1.1.1.1", "ethanol",  "ecoli", 4.0, ()),
    ]
    paths = aggregate_and_write(rows, tmp_path, release="X", date="2026-01-01")
    _, _, data = _read_tsv(paths["kcat"])
    triple_order = [(r[0], r[1], r[2]) for r in data]
    assert triple_order == sorted(triple_order)


def test_byte_identical_across_runs(tmp_path):
    rows1 = list(parse_brenda_json(FIXTURE))
    d1 = tmp_path / "run1"
    d2 = tmp_path / "run2"
    aggregate_and_write(rows1, d1, release="2026.1", date="2026-05-18")
    aggregate_and_write(rows1, d2, release="2026.1", date="2026-05-18")
    for name in ("kcat.tsv", "sa.tsv", "mw.tsv"):
        h1 = hashlib.sha256((d1 / name).read_bytes()).hexdigest()
        h2 = hashlib.sha256((d2 / name).read_bytes()).hexdigest()
        assert h1 == h2, f"{name} not byte-identical across runs"


def test_lf_line_endings(written):
    raw = written["kcat"].read_bytes()
    assert b"\r\n" not in raw, "CRLF line ending detected"


def test_sa_substrate_column_is_star(written):
    _, _, data = _read_tsv(written["sa"])
    if data:
        assert all(r[1] == "*" for r in data)


def test_mw_substrate_column_is_star(written):
    _, _, data = _read_tsv(written["mw"])
    if data:
        assert all(r[1] == "*" for r in data)


def test_default_date_today(tmp_path):
    rows = [Row("kcat", "1.1.1.1", "x", "y", 1.0, ())]
    paths = aggregate_and_write(rows, tmp_path, release="X")
    comment, _, _ = _read_tsv(paths["kcat"])
    import datetime as dt
    assert dt.date.today().isoformat() in comment
