"""Tests for the ComplexPortal JSON loader."""
import json
from pathlib import Path

import pytest

from geckopy.databases import (
    ComplexPortalEntry,
    load_complex_portal_json,
)

EXAMPLE_JSON = (
    Path(__file__).parents[1]
    / "examples" / "ecTestGEM" / "data" / "ComplexPortal.json"
)


def _write(path: Path, payload) -> Path:
    path.write_text(json.dumps(payload))
    return path


def test_loads_ectestgem_fixture():
    entries = load_complex_portal_json(EXAMPLE_JSON)
    assert len(entries) == 1
    e = entries[0]
    assert e.complex_id == "R2Compl"
    assert e.protein_ids == ["P1", "P2"]
    assert e.stoichiometry == [1, 2]
    assert e.gene_names == ["G1", "G2"]


def test_loads_array_of_complexes(tmp_path):
    path = _write(tmp_path / "c.json", [
        {"complexID": "C1", "name": "n1", "specie": "s",
         "geneName": ["G1"], "protID": ["P1"], "stochiometry": [1]},
        {"complexID": "C2", "name": "n2", "specie": "s",
         "geneName": ["G2", "G3"], "protID": ["P2", "P3"],
         "stochiometry": [1, 1]},
    ])
    entries = load_complex_portal_json(path)
    assert len(entries) == 2
    assert entries[1].complex_id == "C2"


def test_accepts_single_object_as_list(tmp_path):
    """If the JSON is a single object, load it as a one-element list."""
    path = _write(tmp_path / "c.json", {
        "complexID": "C1", "name": "n", "specie": "s",
        "geneName": ["G1"], "protID": ["P1"], "stochiometry": [1],
    })
    entries = load_complex_portal_json(path)
    assert len(entries) == 1


def test_renames_stochiometry_to_stoichiometry(tmp_path):
    path = _write(tmp_path / "c.json", [{
        "complexID": "C1", "name": "n", "specie": "s",
        "geneName": ["G1"], "protID": ["P1"], "stochiometry": [3],
    }])
    entries = load_complex_portal_json(path)
    assert entries[0].stoichiometry == [3]


def test_missing_file_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        load_complex_portal_json(tmp_path / "nope.json")


def test_invalid_json_raises(tmp_path):
    path = tmp_path / "bad.json"
    path.write_text("not json {")
    with pytest.raises(ValueError, match="not valid JSON"):
        load_complex_portal_json(path)


def test_missing_required_field_raises(tmp_path):
    """complexID is required; missing it raises ValueError."""
    path = _write(tmp_path / "c.json", [{
        "name": "n", "specie": "s",
        "geneName": [], "protID": [], "stochiometry": [],
    }])
    with pytest.raises(ValueError, match="malformed"):
        load_complex_portal_json(path)


def test_optional_fields_default_to_empty(tmp_path):
    path = _write(tmp_path / "c.json", [{
        "complexID": "C1", "name": "n", "specie": "s",
    }])
    entries = load_complex_portal_json(path)
    assert entries[0].gene_names == []
    assert entries[0].protein_ids == []
    assert entries[0].stoichiometry == []


def test_returns_dataclass_instances(tmp_path):
    entries = load_complex_portal_json(EXAMPLE_JSON)
    assert isinstance(entries[0], ComplexPortalEntry)


def test_scalar_stochiometry_wrapped_to_list(tmp_path):
    """ComplexPortal exports collapse single-element arrays to scalars;
    the loader normalises them back to lists."""
    path = _write(tmp_path / "c.json", [{
        "complexID": "C1", "name": "n", "specie": "s",
        "geneName": ["G1"], "protID": ["P1"], "stochiometry": 0,
    }])
    entries = load_complex_portal_json(path)
    assert entries[0].stoichiometry == [0]


def test_scalar_gene_and_protein_id_wrapped_to_list(tmp_path):
    path = _write(tmp_path / "c.json", [{
        "complexID": "C1", "name": "n", "specie": "s",
        "geneName": "G1", "protID": "P1", "stochiometry": [1],
    }])
    entries = load_complex_portal_json(path)
    assert entries[0].gene_names == ["G1"]
    assert entries[0].protein_ids == ["P1"]
