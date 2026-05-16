"""Tests for the template generator and the CLI scaffold."""
import tomllib

from geckopy import ModelAdapter
from geckopy.adapter.template import generate_template_toml
from geckopy.cli import main


def test_generated_template_parses_as_toml():
    text = generate_template_toml()
    parsed = tomllib.loads(text)
    # Required fields are uncommented, so they parse into the dict.
    assert "conv_gem" in parsed
    assert "org_name" in parsed
    # Optional fields are commented out, so they do not appear.
    assert "sigma" not in parsed


def test_template_without_advanced_excludes_bayesian():
    text = generate_template_toml(advanced=False)
    assert "[bayesian]" not in text
    assert "advanced" in text.lower()


def test_template_with_advanced_includes_bayesian():
    text = generate_template_toml(advanced=True)
    assert "[bayesian]" in text


def test_template_includes_nested_sections():
    text = generate_template_toml()
    assert "[kegg]" in text
    assert "[uniprot]" in text
    assert "[complex]" in text


def test_template_includes_field_descriptions_as_comments():
    text = generate_template_toml()
    assert "saturation" in text.lower()


def test_cli_init_creates_folder_structure(tmp_path):
    target = tmp_path / "my-ec"
    exit_code = main(["init", str(target)])
    assert exit_code == 0
    assert (target / "model_adapter.toml").is_file()
    assert (target / "adapter.py").is_file()
    assert (target / "data").is_dir()
    assert (target / "models").is_dir()
    assert (target / "output").is_dir()


def test_cli_init_refuses_nonempty_folder(tmp_path, capsys):
    target = tmp_path / "existing"
    target.mkdir()
    (target / "something.txt").write_text("hello")

    exit_code = main(["init", str(target)])
    assert exit_code == 1
    captured = capsys.readouterr()
    assert "not empty" in captured.err


def test_cli_init_accepts_advanced_flag(tmp_path):
    target = tmp_path / "my-ec"
    main(["init", str(target), "--advanced"])
    toml_text = (target / "model_adapter.toml").read_text()
    assert "[bayesian]" in toml_text


def test_cli_init_produces_loadable_adapter_after_filling_required(tmp_path):
    """End-to-end: scaffold, fill in required fields, load via ModelAdapter."""
    target = tmp_path / "my-ec"
    main(["init", str(target)])

    toml_path = target / "model_adapter.toml"
    toml_text = toml_path.read_text()
    toml_text = toml_text.replace(
        '"models/REPLACE_ME.xml"', '"models/real.xml"'
    ).replace('"REPLACE_ME"', '"Testus testus"')
    toml_path.write_text(toml_text)

    adapter = ModelAdapter.from_folder(target)
    assert adapter.params.org_name == "Testus testus"
    # Defaults kick in for fields the user did not uncomment.
    assert adapter.params.sigma == 0.5
    assert adapter.params.kegg.id == "sce"


def test_cli_init_stub_adapter_py_is_commented_out(tmp_path):
    """Generated adapter.py should define nothing at import time."""
    target = tmp_path / "my-ec"
    main(["init", str(target)])
    namespace: dict = {}
    exec((target / "adapter.py").read_text(), namespace)
    user_classes = [
        v for k, v in namespace.items()
        if not k.startswith("_") and isinstance(v, type)
    ]
    assert user_classes == []
