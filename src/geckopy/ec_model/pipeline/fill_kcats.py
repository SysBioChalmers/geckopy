"""Fill missing kcats by averaging across isozymes.

Ported from GECKO MATLAB:
src/geckomat/change_model/getKcatAcrossIsozymes.m.
"""
from __future__ import annotations

import logging
import re
from typing import TYPE_CHECKING

import numpy as np

from .apply_kcat import apply_kcat_constraints

if TYPE_CHECKING:
    from ..ec_model import EcModel


logger = logging.getLogger(__name__)

_EXP_SUFFIX_REGEX = re.compile(r"_EXP_\d+$")
_SOURCE_TAG = "isozymes"


def get_kcat_across_isozymes(
    model: "EcModel",
    *,
    apply: bool = True,
) -> None:
    """Fill in missing kcats by averaging known kcats across sibling isozymes.

    Ported from GECKO MATLAB:
    src/geckomat/change_model/getKcatAcrossIsozymes.m.

    Two reactions are siblings if their ``ec.rxns`` IDs are identical
    after stripping any ``_EXP_<n>`` suffix. The ``_REV`` suffix is
    NOT stripped: forward and reverse reactions share chemistry but
    can have different kcats, so they are kept distinct.

    For each reaction with NaN ``ec.kcat`` (missing), find its siblings
    with non-NaN kcats. If any exist, take their mean and assign it.
    Set ``ec.source`` to ``"isozymes"`` for the filled entries.

    Reactions whose siblings all also have missing kcats remain NaN.
    Reactions with a single isozyme (no siblings beyond themselves)
    also remain NaN.

    MATLAB-COMPAT: MATLAB uses ``kcat == 0`` to mark "missing"; geckopy
    uses ``NaN``. See docs/future_improvements.md for translation
    discussion.

    MATLAB-COMPAT: The function name "get_kcat_across_isozymes" is a
    literal port; "fill_kcats_from_isozymes" would be more descriptive
    since the function modifies ec.kcat in place. A coordinated rename
    is tracked in docs/future_improvements.md.

    Parameters
    ----------
    model
        An EcModel produced by ``make_ec_model``. Mutated in place.
        Must not be a gecko-light model.
    apply
        If True (default), call ``apply_kcat_constraints`` after
        updating ec.kcat so the new values reflect in the S matrix.

    Raises
    ------
    NotImplementedError
        If ``model.ec.gecko_light`` is True. Isozyme-averaging does not
        apply to the light formulation since it does not split
        reactions per isozyme.
    """
    if model.ec.gecko_light:
        raise NotImplementedError(
            "get_kcat_across_isozymes: not applicable to gecko-light models."
        )

    kcat = model.ec.kcat
    if kcat.size == 0 or np.isnan(kcat).all():
        logger.warning(
            "get_kcat_across_isozymes: ec.kcat has no known values to "
            "average from; model unchanged."
        )
        return

    # Strip _EXP_<n> suffix to identify isozyme groups. _REV stays.
    base_ids = [_EXP_SUFFIX_REGEX.sub("", r) for r in model.ec.rxns]

    # Group known kcats by base_id.
    known_by_base: dict[str, list[float]] = {}
    for bid, k in zip(base_ids, kcat):
        if not np.isnan(k):
            known_by_base.setdefault(bid, []).append(float(k))

    # For each missing entry, look up its base group and average if available.
    filled_indices: list[int] = []
    for i, (bid, k) in enumerate(zip(base_ids, kcat)):
        if not np.isnan(k):
            continue
        siblings = known_by_base.get(bid)
        if not siblings:
            continue
        model.ec.kcat[i] = float(np.mean(siblings))
        model.ec.source[i] = _SOURCE_TAG
        filled_indices.append(i)

    if not filled_indices:
        logger.info(
            "get_kcat_across_isozymes: no missing kcats had isozymes with "
            "known kcats; nothing filled."
        )
        return

    logger.info(
        "get_kcat_across_isozymes: filled %d kcat(s) by averaging across "
        "isozymes.", len(filled_indices),
    )

    if apply:
        filled_rxn_ids = [model.ec.rxns[i] for i in filled_indices]
        apply_kcat_constraints(model, update_rxns=filled_rxn_ids)
