"""Tests for canonicalize_rxn_id (shared _REV/_EXP suffix stripping)."""
import pytest

from geckopy.ec_model.constants import canonicalize_rxn_id


@pytest.mark.parametrize(
    "rxn_id, expected",
    [
        ("R2", ("R2", False)),
        ("R2_EXP_1", ("R2", False)),
        ("R2_EXP_12", ("R2", False)),
        ("R2_REV", ("R2", True)),
        ("R2_REV_EXP_3", ("R2", True)),
        # ids that merely contain "_REV" must not be mangled
        ("R_REVERTASE", ("R_REVERTASE", False)),
        ("FOO_REV_BAR", ("FOO_REV_BAR", False)),
        # only a trailing _REV is stripped, not an internal one
        ("R_REV_X_REV", ("R_REV_X", True)),
    ],
)
def test_canonicalize_rxn_id(rxn_id, expected):
    assert canonicalize_rxn_id(rxn_id) == expected
