"""Tests for the ModelAdapter and ModelParameters."""
from pathlib import Path

import pytest
from pydantic import ValidationError

from geckopy import ModelAdapter

EXAMPLES_DIR = Path(__file__).resolve().parent.parent / "examples"


MINIMAL_TOML = """\
conv_gem = "models/test.xml"
org_name = "Testus organismus"
"""

FULL_TOML = """\
conv_gem = "models/test.xml"
org_name = "Saccharomyces cerevisiae"
sigma = 0.5
p_tot = 0.5
f = 0.4461
gr_exp = 0.41

[kegg]
id = "sce"

[uniprot]
type = "taxonomy"
id = "559292"
reviewed = true

[complex]
taxonomic_id = 559292
"""


def _write_adapter_folder(tmp_path: Path, toml_content: str) -> Path:
    """Create a temporary adapter folder with the given TOML content."""
    (tmp_path / "model_adapter.toml").write_text(toml_content)
    return tmp_path


def test_from_folder_with_minimal_toml(tmp_path):
    folder = _write_adapter_folder(tmp_path, MINIMAL_TOML)
    adapter = ModelAdapter.from_folder(folder)

    assert adapter.params.org_name == "Testus organismus"
    assert adapter.params.sigma == 0.5  # default
    assert adapter.params.path == folder.resolve()
    assert adapter.params.conv_gem == (folder / "models/test.xml").resolve()


def test_from_folder_with_full_toml(tmp_path):
    folder = _write_adapter_folder(tmp_path, FULL_TOML)
    adapter = ModelAdapter.from_folder(folder)

    assert adapter.params.org_name == "Saccharomyces cerevisiae"
    assert adapter.params.f == 0.4461
    assert adapter.params.kegg.id == "sce"
    assert adapter.params.uniprot.reviewed is True
    assert adapter.params.uniprot.id == "559292"


def test_missing_toml_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        ModelAdapter.from_folder(tmp_path)


def test_unknown_field_rejected(tmp_path):
    bad_toml = MINIMAL_TOML + '\nsgma = 0.3\n'  # typo
    folder = _write_adapter_folder(tmp_path, bad_toml)
    with pytest.raises(ValidationError):
        ModelAdapter.from_folder(folder)


def test_missing_required_field_rejected(tmp_path):
    # org_name is required, leaving it out should fail
    bad_toml = 'conv_gem = "models/test.xml"\n'
    folder = _write_adapter_folder(tmp_path, bad_toml)
    with pytest.raises(ValidationError):
        ModelAdapter.from_folder(folder)


def test_default_overridable_methods(tmp_path):
    folder = _write_adapter_folder(tmp_path, MINIMAL_TOML)
    adapter = ModelAdapter.from_folder(folder)

    # Defaults should behave like the MATLAB base class
    assert adapter.get_spontaneous_reactions(model=None) == []
    assert adapter.get_uniprot_compatible_genes(["YAL001C"]) == ["YAL001C"]


def test_uniprot_table_lookup_without_file(tmp_path):
    folder = _write_adapter_folder(tmp_path, MINIMAL_TOML)
    adapter = ModelAdapter.from_folder(folder)

    result = adapter.get_uniprot_ids_from_table(["YAL001C", "YAL002W"])
    assert result == ["YAL001C", "YAL002W"]


def test_uniprot_table_lookup_with_file(tmp_path):
    folder = _write_adapter_folder(tmp_path, MINIMAL_TOML)
    data_dir = folder / "data"
    data_dir.mkdir()
    (data_dir / "uniprotConversion.tsv").write_text(
        "gene_id\tuniprot_id\n"
        "YAL001C\tP12345\n"
        "YAL002W\tP67890\n"
    )

    adapter = ModelAdapter.from_folder(folder)
    result = adapter.get_uniprot_ids_from_table(["YAL001C", "YAL002W", "YAL999W"])
    assert result == ["P12345", "P67890", "YAL999W"]  # unknown gene passes through


def test_subclass_from_folder_returns_subclass(tmp_path):
    folder = _write_adapter_folder(tmp_path, MINIMAL_TOML)

    class CustomAdapter(ModelAdapter):
        def get_spontaneous_reactions(self, model):
            return ["spont_rxn_1"]

    adapter = CustomAdapter.from_folder(folder)
    assert isinstance(adapter, CustomAdapter)
    assert adapter.get_spontaneous_reactions(model=None) == ["spont_rxn_1"]


@pytest.mark.parametrize(
    "example_dir", [p for p in EXAMPLES_DIR.iterdir() if p.is_dir()],
    ids=lambda p: p.name,
)
def test_shipped_example_adapter_toml_is_valid(example_dir):
    """Every folder under examples/ must have a model_adapter.toml that
    actually validates against the current ModelParameters schema.

    examples/yeast-GEM/model_adapter.toml drifted from the schema
    (`uniprot.tax_id`, rejected by `extra="forbid"`, instead of
    `uniprot.id`) and nothing caught it until a downstream project tried
    to load it directly.
    """
    ModelAdapter.from_folder(example_dir)
