"""Smoke test that every name in ``geckopy.__all__`` resolves."""


def test_top_level_exports():
    import geckopy
    for name in geckopy.__all__:
        assert hasattr(geckopy, name), (
            f"Missing top-level export: {name}"
        )
        # Also assert it's not None (covers the SBML import fallback path).
        assert getattr(geckopy, name) is not None, (
            f"Top-level export {name} is None (import fell back)"
        )
