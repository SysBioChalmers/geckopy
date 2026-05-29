"""Tests for download_kegg using mocked HTTP."""
from __future__ import annotations

import csv
from pathlib import Path
from unittest.mock import MagicMock

import pytest

from geckopy.databases.kegg_download import download_kegg


def _entry_text(
    kegg_id: str,
    uniprot: str | None,
    ec: str | None,
    pathway: str | None,
    sequence: str | None,
    *,
    extra_field: str = "",
) -> str:
    """Build a synthetic KEGG flat-file entry."""
    lines = [
        f"ENTRY       {kegg_id}           CDS       T00005",
        f"NAME        gene_{kegg_id}",
    ]
    if extra_field:
        lines.append(extra_field)
    if ec:
        lines.append(
            f"ORTHOLOGY   K00001  alcohol dehydrogenase [EC:{ec}]"
        )
    if pathway:
        lines.append(f"PATHWAY     {pathway}")
        lines.append("BRITE       enzyme classification")
    if uniprot:
        lines.append(f"DBLINKS     UniProt: {uniprot}")
    if sequence:
        lines.append(f"AASEQ       {len(sequence)}")
        lines.append("            " + sequence)
    return "\n".join(lines) + "\n"


def _make_session(
    list_response: str, get_response: str
) -> MagicMock:
    session = MagicMock()

    def fake_get(url, timeout=None):
        resp = MagicMock()
        resp.status_code = 200
        resp.raise_for_status = MagicMock()
        if "/list/" in url:
            resp.text = list_response
        elif "/get/" in url:
            resp.text = get_response
        else:
            resp.text = ""
        return resp

    session.get.side_effect = fake_get
    return session


def test_basic_download_two_genes(tmp_path):
    list_resp = "sce:G1\tgene 1\nsce:G2\tgene 2\n"
    get_resp = "\n///\n".join([
        _entry_text("G1", "P00001", "1.1.1.1", "sce00010 Glycolysis", "MVQR"),
        _entry_text("G2", "P00002", "2.2.2.2", "sce00020 TCA", "MAEK"),
    ]) + "\n///\n"
    sess = _make_session(list_resp, get_resp)

    out = tmp_path / "kegg.tsv"
    download_kegg("sce", out, session=sess)
    rows = list(csv.reader(out.open("r", encoding="utf-8")))
    assert len(rows) == 2
    assert rows[0][:4] == ["P00001", "G1", "G1", "1.1.1.1"]
    assert rows[1][:4] == ["P00002", "G2", "G2", "2.2.2.2"]
    assert rows[0][6] == "MVQR"


def test_entry_without_aaseq_dropped(tmp_path):
    list_resp = "sce:G1\tgene 1\nsce:G2\tgene 2\n"
    get_resp = "\n///\n".join([
        _entry_text("G1", "P00001", "1.1.1.1", "p", "MVQR"),
        _entry_text("G2", "P00002", "2.2.2.2", "p", None),
    ]) + "\n///\n"
    sess = _make_session(list_resp, get_resp)

    out = tmp_path / "kegg.tsv"
    download_kegg("sce", out, session=sess)
    rows = list(csv.reader(out.open("r", encoding="utf-8")))
    assert len(rows) == 1
    assert rows[0][1] == "G1"


def test_entry_without_uniprot_dropped(tmp_path):
    list_resp = "sce:G1\tgene 1\nsce:G2\tgene 2\n"
    get_resp = "\n///\n".join([
        _entry_text("G1", "P00001", "1.1.1.1", "p", "MVQR"),
        _entry_text("G2", None, "2.2.2.2", "p", "MAEK"),
    ]) + "\n///\n"
    sess = _make_session(list_resp, get_resp)

    out = tmp_path / "kegg.tsv"
    download_kegg("sce", out, session=sess)
    rows = list(csv.reader(out.open("r", encoding="utf-8")))
    assert len(rows) == 1
    assert rows[0][1] == "G1"


def test_batching_more_than_ten_genes(tmp_path):
    gene_ids = [f"G{i}" for i in range(12)]
    list_resp = "\n".join(f"sce:{g}\tgene {g}" for g in gene_ids) + "\n"
    sess = MagicMock()
    calls = {"get_url": []}

    def fake_get(url, timeout=None):
        resp = MagicMock()
        resp.status_code = 200
        resp.raise_for_status = MagicMock()
        if "/list/" in url:
            resp.text = list_resp
        else:
            calls["get_url"].append(url)
            requested = [
                seg.split(":")[1] for seg in url.rsplit("/", 1)[1].split("+")
            ]
            entries = [
                _entry_text(g, f"P{g}", "1.1.1.1", "p", "MVQR")
                for g in requested
            ]
            resp.text = "\n///\n".join(entries) + "\n///\n"
        return resp

    sess.get.side_effect = fake_get
    out = tmp_path / "kegg.tsv"
    download_kegg("sce", out, session=sess)
    rows = list(csv.reader(out.open("r", encoding="utf-8")))
    assert len(rows) == 12
    # Two GET-batch calls expected for 12 genes / 10 per batch
    assert len(calls["get_url"]) == 2


def test_mw_computed_from_sequence(tmp_path):
    list_resp = "sce:G1\tgene\n"
    get_resp = (
        _entry_text("G1", "P1", "1.1.1.1", "p", "MAEKR") + "///\n"
    )
    sess = _make_session(list_resp, get_resp)
    out = tmp_path / "kegg.tsv"
    download_kegg("sce", out, session=sess)
    rows = list(csv.reader(out.open("r", encoding="utf-8")))
    mw = float(rows[0][4])
    # MAEKR: M=131.20, A=71.08, E=129.11, K=128.17, R=156.19, +water 18.02
    # ~= 633.77 -> round -> 634
    assert mw == pytest.approx(634, abs=2)


def test_pathway_multiline_joined(tmp_path):
    list_resp = "sce:G1\tgene\n"
    # PATHWAY block followed by additional non-PATHWAY content
    pathway_block = (
        "sce00010  Glycolysis\n"
        "            sce00020  TCA cycle"
    )
    get_resp = (
        _entry_text("G1", "P1", "1.1.1.1", pathway_block, "MVQR") + "///\n"
    )
    sess = _make_session(list_resp, get_resp)
    out = tmp_path / "kegg.tsv"
    download_kegg("sce", out, session=sess)
    rows = list(csv.reader(out.open("r", encoding="utf-8")))
    pathway = rows[0][5]
    assert "Glycolysis" in pathway
    assert "TCA cycle" in pathway


def test_empty_list_response_raises(tmp_path):
    sess = _make_session("", "")
    with pytest.raises(RuntimeError, match="no genes"):
        download_kegg("badcode", tmp_path / "kegg.tsv", session=sess)


def test_invalid_organism_code_raises(tmp_path):
    sess = MagicMock()

    def fake_get(url, timeout=None):
        resp = MagicMock()
        resp.status_code = 400
        resp.text = ""
        resp.raise_for_status = MagicMock(side_effect=Exception("400"))
        return resp

    sess.get.side_effect = fake_get
    with pytest.raises(RuntimeError, match="400"):
        download_kegg("badcode", tmp_path / "kegg.tsv", session=sess)


def test_entry_without_ec_emits_empty_ec(tmp_path):
    list_resp = "sce:G1\tgene\n"
    get_resp = (
        _entry_text("G1", "P1", None, "p", "MVQR") + "///\n"
    )
    sess = _make_session(list_resp, get_resp)
    out = tmp_path / "kegg.tsv"
    download_kegg("sce", out, session=sess)
    rows = list(csv.reader(out.open("r", encoding="utf-8")))
    assert rows[0][3] == ""


def test_round_trip_with_load_kegg_tsv(tmp_path):
    """Download via mocked session, then parse with load_kegg_tsv."""
    from geckopy.databases.kegg_loader import load_kegg_tsv
    list_resp = "sce:G1\tgene\n"
    get_resp = (
        _entry_text("G1", "P00001", "1.1.1.1", "sce00010 Glycolysis", "MAEK")
        + "///\n"
    )
    sess = _make_session(list_resp, get_resp)
    out = tmp_path / "kegg.tsv"
    download_kegg("sce", out, session=sess)
    db = load_kegg_tsv(out, auto_download=False)
    assert len(db) == 1
    assert db.uniprot_ids == ["P00001"]
    assert db.genes == ["G1"]
    assert db.kegg_genes == ["G1"]
    assert db.eccodes == ["1.1.1.1"]
    assert db.sequences == ["MAEK"]
