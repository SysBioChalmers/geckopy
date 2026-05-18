"""Unit tests for BRENDA bulk-JSON cache + extraction."""
from __future__ import annotations

import io
import json
import tarfile
from pathlib import Path

import pytest

from geckopy.databases.brenda.download import (
    BrendaDownloadError,
    ensure_brenda_json,
    extract_brenda_json,
    sha256_of,
)


def _make_tarball(tar_path: Path, json_name: str, payload: dict) -> bytes:
    raw = json.dumps(payload).encode("utf-8")
    with tarfile.open(tar_path, "w:gz") as tar:
        info = tarfile.TarInfo(name=json_name)
        info.size = len(raw)
        tar.addfile(info, io.BytesIO(raw))
    return tar_path.read_bytes()


def test_empty_cache_raises_with_instructions(tmp_path):
    with pytest.raises(BrendaDownloadError, match="download.php"):
        ensure_brenda_json(tmp_path)


def test_creates_cache_dir_if_missing(tmp_path):
    missing = tmp_path / "does-not-exist"
    with pytest.raises(BrendaDownloadError):
        ensure_brenda_json(missing)
    assert missing.is_dir()


def test_returns_existing_unpacked_json(tmp_path):
    target = tmp_path / "brenda_2026_1.json"
    target.write_text('{"release": "2026.1"}', encoding="utf-8")
    result = ensure_brenda_json(tmp_path)
    assert result == target.resolve()


def test_extracts_tarball_when_only_tarball_present(tmp_path):
    tar = tmp_path / "brenda_2026_1.json.tar.gz"
    _make_tarball(tar, "brenda_2026_1.json", {"release": "2026.1"})
    result = ensure_brenda_json(tmp_path)
    assert result.suffix == ".json"
    assert result.exists()
    assert json.loads(result.read_text())["release"] == "2026.1"


def test_sha256_validation_passes(tmp_path):
    tar = tmp_path / "brenda_2026_1.json.tar.gz"
    _make_tarball(tar, "brenda_2026_1.json", {"release": "2026.1"})
    digest = sha256_of(tar)
    result = ensure_brenda_json(tmp_path, expected_sha256=digest)
    assert result.exists()


def test_sha256_validation_fails(tmp_path):
    tar = tmp_path / "brenda_2026_1.json.tar.gz"
    _make_tarball(tar, "brenda_2026_1.json", {"release": "2026.1"})
    with pytest.raises(BrendaDownloadError, match="sha256 mismatch"):
        ensure_brenda_json(tmp_path, expected_sha256="0" * 64)


def test_tarball_with_directory_prefix_flattened(tmp_path):
    tar = tmp_path / "brenda_2026_1.json.tar.gz"
    _make_tarball(tar, "subdir/brenda_2026_1.json", {"release": "2026.1"})
    result = ensure_brenda_json(tmp_path)
    assert result.name == "brenda_2026_1.json"
    assert result.parent == tmp_path.resolve()


def test_tarball_without_json_member_raises(tmp_path):
    tar = tmp_path / "brenda_2026_1.json.tar.gz"
    with tarfile.open(tar, "w:gz") as t:
        info = tarfile.TarInfo(name="readme.txt")
        info.size = 0
        t.addfile(info, io.BytesIO(b""))
    with pytest.raises(BrendaDownloadError, match="no .json member"):
        ensure_brenda_json(tmp_path)


def test_tarball_with_multiple_json_members_raises(tmp_path):
    tar = tmp_path / "brenda_2026_1.json.tar.gz"
    with tarfile.open(tar, "w:gz") as t:
        for name in ("a.json", "b.json"):
            raw = b"{}"
            info = tarfile.TarInfo(name=name)
            info.size = len(raw)
            t.addfile(info, io.BytesIO(raw))
    with pytest.raises(BrendaDownloadError, match="exactly one .json member"):
        ensure_brenda_json(tmp_path)


def test_unpacked_json_takes_precedence_over_tarball(tmp_path):
    tar = tmp_path / "brenda_2026_1.json.tar.gz"
    _make_tarball(tar, "brenda_2026_1.json", {"release": "wrong"})
    direct = tmp_path / "brenda_2026_1.json"
    direct.write_text('{"release": "right"}', encoding="utf-8")
    result = ensure_brenda_json(tmp_path)
    assert json.loads(result.read_text())["release"] == "right"


def test_extract_brenda_json_returns_path(tmp_path):
    tar = tmp_path / "bundle.tar.gz"
    _make_tarball(tar, "data.json", {"k": "v"})
    out = extract_brenda_json(tar, tmp_path)
    assert out.name == "data.json"
    assert json.loads(out.read_text())["k"] == "v"
