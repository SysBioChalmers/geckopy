"""Set kcat values for individual reactions, with optional auto-apply.

Ported from GECKO MATLAB: src/geckomat/change_model/setKcatForReactions.m.
"""
from __future__ import annotations

import re
from typing import TYPE_CHECKING, Sequence, Union

import numpy as np

from .apply_kcat import apply_kcat_constraints

if TYPE_CHECKING:
    from ..ec_model import EcModel


_EXP_SUFFIX_PATTERN = re.compile(r"_EXP_\d+$")
_SOURCE_TAG = "manual"


def set_kcat_for_reactions(
    model: "EcModel",
    rxn_ids: Sequence[str],
    kcat: Union[float, Sequence[float]],
    *,
    apply: bool = True,
) -> list[str]:
    """Set kcat values for one or more reactions in ec.kcat.

    Ported from GECKO MATLAB: src/geckomat/change_model/setKcatForReactions.m.

    Each ID in ``rxn_ids`` is interpreted as follows:

    - If the ID ends in ``_EXP_<n>``, only that exact reaction is matched.
    - Otherwise, the ID is treated as a base name and matches every
      ec.rxns entry whose ID, after stripping any ``_EXP_<n>`` suffix,
      equals the base name. So ``"R2"`` matches ``R2``, ``R2_EXP_1``,
      ``R2_EXP_2``, etc.

    The ``kcat`` argument follows numpy-style broadcasting: a single
    float applies to every matched reaction; a sequence must match the
    total number of matched reactions across all ``rxn_ids``.

    Strict rule: when one un-suffixed ``rxn_id`` expands to multiple
    matches (isozymes), the kcat value for that ID must be a scalar.
    Passing a length-N kcat list to cover N expansions implicitly
    relies on the order of ``ec.rxns``, which is fragile. To set
    different values for different isozymes, pass the suffixed IDs
    explicitly.

    MATLAB-COMPAT: MATLAB allows passing a length-N kcat for an
    un-suffixed ID that expands to N matches; geckopy forbids this.
    MATLAB GECKO should adopt the strict rule.

    MATLAB-COMPAT: MATLAB writes the source string ``'setKcatForReactions'``
    to ec.source for changed reactions. geckopy writes ``'manual'``.
    MATLAB GECKO should adopt ``'manual'`` for round-trippable source
    strings.

    Parameters
    ----------
    model
        An EcModel with ec.rxns populated. Mutated in place.
    rxn_ids
        Reaction identifiers to update.
    kcat
        New kcat value(s) in 1/s. Scalar or sequence (see above).
    apply
        If True (default), call ``apply_kcat_constraints`` after
        updating ec.kcat / ec.source so the change is reflected in the
        S matrix immediately. If False, the caller must invoke
        ``apply_kcat_constraints`` themselves.

    Returns
    -------
    list of str
        IDs of the reactions whose kcat was changed (the post-expansion
        set), in the order in which they were updated.

    Raises
    ------
    ValueError
        If any rxn_id matches zero reactions, or if kcat lengths do not
        match, or if an un-suffixed rxn_id expands to multiple matches
        but its kcat is given as a list.
    """
    rxn_ids = list(rxn_ids)
    if not rxn_ids:
        return []

    ec_rxns = model.ec.rxns
    nonexp_ec_rxns = [_EXP_SUFFIX_PATTERN.sub("", r) for r in ec_rxns]

    # Resolve every input ID to a list of ec-row indices.
    matches_per_input: list[list[int]] = []
    for rxn_id in rxn_ids:
        if _EXP_SUFFIX_PATTERN.search(rxn_id):
            indices = [i for i, r in enumerate(ec_rxns) if r == rxn_id]
        else:
            indices = [
                i for i, r in enumerate(nonexp_ec_rxns) if r == rxn_id
            ]
        if not indices:
            raise ValueError(
                f"rxn_id '{rxn_id}' matched no entries in ec.rxns."
            )
        matches_per_input.append(indices)

    total_matches = sum(len(m) for m in matches_per_input)

    # Resolve kcat to one scalar per ec-row index in resolution order.
    if isinstance(kcat, (int, float)) and not isinstance(kcat, bool):
        kcat_per_index = [float(kcat)] * total_matches
    else:
        kcat_seq = list(kcat)
        # kcat_seq has length len(rxn_ids) (one value per input ID, each
        # of which may broadcast to multiple matches), or length
        # total_matches (one value per resolved index, only allowed when
        # every input has a single match).
        if len(kcat_seq) == len(rxn_ids):
            kcat_per_index = []
            for value, indices in zip(kcat_seq, matches_per_input):
                kcat_per_index.extend([float(value)] * len(indices))
        elif len(kcat_seq) == total_matches:
            for rxn_id, indices in zip(rxn_ids, matches_per_input):
                if len(indices) > 1:
                    raise ValueError(
                        f"rxn_id '{rxn_id}' expands to {len(indices)} "
                        f"isozymes; kcat for it must be a scalar. To set "
                        f"different values per isozyme, pass the suffixed "
                        f"reaction IDs explicitly."
                    )
            kcat_per_index = [float(v) for v in kcat_seq]
        else:
            raise ValueError(
                f"kcat has length {len(kcat_seq)}; expected a scalar, "
                f"length {len(rxn_ids)} (one per input rxn_id), or length "
                f"{total_matches} (one per resolved match)."
            )

    # Apply to ec.kcat and ec.source.
    flat_indices: list[int] = []
    for indices in matches_per_input:
        flat_indices.extend(indices)

    for idx, value in zip(flat_indices, kcat_per_index):
        model.ec.kcat[idx] = value
        model.ec.source[idx] = _SOURCE_TAG

    updated_ids = [ec_rxns[i] for i in flat_indices]

    if apply:
        apply_kcat_constraints(model, update_rxns=updated_ids)

    return updated_ids
