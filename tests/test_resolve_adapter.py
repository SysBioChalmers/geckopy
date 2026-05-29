"""Tests for the shared ``resolve_adapter`` helper."""
import pytest

from geckopy.adapter import ModelAdapter, resolve_adapter


class _ModelStub:
    """Minimal stand-in: anything with an ``.adapter`` attribute."""

    def __init__(self, adapter=None):
        self.adapter = adapter


def _make_adapter(folder):
    folder.mkdir(parents=True, exist_ok=True)
    (folder / "model_adapter.toml").write_text(
        'conv_gem = "x.xml"\norg_name = "test"\n'
    )
    return ModelAdapter.from_folder(folder)


def test_explicit_adapter_takes_priority(tmp_path):
    """If both are non-None, the explicit ``adapter`` argument wins."""
    a1 = _make_adapter(tmp_path / "a")
    a2 = _make_adapter(tmp_path / "b")
    model = _ModelStub(adapter=a1)
    assert resolve_adapter(model, a2) is a2


def test_falls_back_to_model_adapter(tmp_path):
    """With no explicit argument, ``model.adapter`` is returned."""
    a = _make_adapter(tmp_path)
    model = _ModelStub(adapter=a)
    assert resolve_adapter(model) is a


def test_explicit_none_falls_through_to_model(tmp_path):
    """Passing ``adapter=None`` is the same as omitting it."""
    a = _make_adapter(tmp_path)
    model = _ModelStub(adapter=a)
    assert resolve_adapter(model, None) is a


def test_raises_when_both_are_none():
    """The whole point: surface a clear error when nothing's available."""
    model = _ModelStub(adapter=None)
    with pytest.raises(ValueError, match="No ModelAdapter available"):
        resolve_adapter(model)


def test_purpose_appears_in_error_message():
    """The ``purpose`` keyword is appended to help the caller diagnose."""
    model = _ModelStub(adapter=None)
    with pytest.raises(ValueError, match="reading params.foo from the adapter"):
        resolve_adapter(model, purpose="reading params.foo from the adapter")


def test_model_without_adapter_attribute_handled():
    """Plain cobra models don't have an ``.adapter`` attribute. The
    helper should treat that as None, not crash."""
    class _Plain:
        pass

    with pytest.raises(ValueError, match="No ModelAdapter available"):
        resolve_adapter(_Plain())
