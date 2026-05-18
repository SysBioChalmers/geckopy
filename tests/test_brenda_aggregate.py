"""Unit tests for BRENDA aggregation + TSV writing."""
from __future__ import annotations

import hashlib
from pathlib import Path

import pytest

from geckopy.databases.brenda import Row, parse_brenda_json
from geckopy.databases.brenda.aggregate import aggregate_and_write

FIXTURE = Path(__file__).parent / "data" / "brenda_minimal.json"


def _read_tsv(path: Path) -> tuple[str, list[list[str]]]:
    text = path.read_text(encoding="utf-8")
    lines = text.split("\n")
    header = lines[0]
    data = [line.split("\t") for line in lines[1:] if line]
    return header, data


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


def test_header_line_format(written):
    header, _ = _read_tsv(written["kcat"])
    assert header.startswith("# BRENDA release 2026.1 generated 2026-05-18 - CC BY 4.0")
    assert "kcat in 1/s" in header


def test_kcat_columns_and_values(written):
    _, data = _read_tsv(written["kcat"])
    assert all(len(row) == 5 for row in data)
    sce_row = [r for r in data if r[2] == "saccharomyces cerevisiae" and r[1] == "ethanol"]
    assert len(sce_row) == 1
    assert sce_row[0][0] == "1.1.1.1"
    assert sce_row[0][3] == "23.5"
    assert sce_row[0][4] == "PMID:11111;PMID:22222"


def test_max_value_wins_on_collision(tmp_path):
    rows = [
        Row("kcat", "1.1.1.1", "ethanol", "yeast", 10.0, ("PMID:1",)),
        Row("kcat", "1.1.1.1", "ethanol", "yeast", 50.0, ("PMID:2",)),
        Row("kcat", "1.1.1.1", "ethanol", "yeast", 30.0, ("PMID:3",)),
    ]
    paths = aggregate_and_write(rows, tmp_path, release="X", date="2026-01-01")
    _, data = _read_tsv(paths["kcat"])
    assert len(data) == 1
    assert data[0][3] == "50.0"
    assert data[0][4] == "PMID:1;PMID:2;PMID:3"


def test_empty_references_emit_star(tmp_path):
    rows = [Row("kcat", "1.1.1.1", "ethanol", "yeast", 10.0, ())]
    paths = aggregate_and_write(rows, tmp_path, release="X", date="2026-01-01")
    _, data = _read_tsv(paths["kcat"])
    assert data[0][4] == "*"


def test_deterministic_sort(tmp_path):
    rows = [
        Row("kcat", "2.7.1.1", "glucose", "ecoli", 1.0, ()),
        Row("kcat", "1.1.1.1", "methanol", "yeast", 2.0, ()),
        Row("kcat", "1.1.1.1", "ethanol",  "yeast", 3.0, ()),
        Row("kcat", "1.1.1.1", "ethanol",  "ecoli", 4.0, ()),
    ]
    paths = aggregate_and_write(rows, tmp_path, release="X", date="2026-01-01")
    _, data = _read_tsv(paths["kcat"])
    ordered = [(r[0], r[1], r[2]) for r in data]
    assert ordered == sorted(ordered)


def test_byte_identical_across_runs(tmp_path):
    rows1 = list(parse_brenda_json(FIXTURE))
    d1 = tmp_path / "run1"
    d2 = tmp_path / "run2"
    aggregate_and_write(rows1, d1, release="2026.1", date="2026-05-18")
    aggregate_and_write(rows1, d2, release="2026.1", date="2026-05-18")
    for name in ("max_kcat.tsv", "max_sa.tsv", "max_mw.tsv"):
        h1 = hashlib.sha256((d1 / name).read_bytes()).hexdigest()
        h2 = hashlib.sha256((d2 / name).read_bytes()).hexdigest()
        assert h1 == h2, f"{name} not byte-identical across runs"


def test_lf_line_endings(written):
    raw = written["kcat"].read_bytes()
    assert b"\r\n" not in raw, "CRLF line ending detected"


def test_sa_and_mw_substrate_column_is_star(written):
    for kind in ("sa", "mw"):
        _, data = _read_tsv(written[kind])
        if data:
            assert all(r[1] == "*" for r in data)


def test_default_date_today(tmp_path):
    rows = [Row("kcat", "1.1.1.1", "x", "y", 1.0, ())]
    paths = aggregate_and_write(rows, tmp_path, release="X")
    header, _ = _read_tsv(paths["kcat"])
    import datetime as dt
    assert dt.date.today().isoformat() in header
