"""Populate model.ec.eccodes from per-reaction `ec-code` annotations.

Ported from GECKO MATLAB:
src/geckomat/get_enzyme_data/getECfromGEM.m.
"""
from __future__ import annotations

import logging
import re
from typing import TYPE_CHECKING, Iterable, Optional

if TYPE_CHECKING:
    from ..ec_model.ec_model import EcModel

logger = logging.getLogger(__name__)


# Pattern for a `;`-joined EC string, e.g. `1.2.3.4;5.6.7.-`. Each token
# is four dot-separated levels of either a non-negative integer or `-`
# (the IUBMB placeholder). The constant is duplicated from
# `ec_string.py` (token-level there, semicolon-joined here) by design;
# the two patterns differ in their separators and a shared helper would
# obscure that.
_EC_TOKEN = r"(?:\d+|-)\.(?:\d+|-)\.(?:\d+|-)\.(?:\d+|-)"
_EC_MULTI_RE = re.compile(rf"^{_EC_TOKEN}(?:;{_EC_TOKEN})*$")


def _normalize_annotation(value) -> str:
    """Convert a cobrapy `ec-code` annotation value to a `;`-joined string.

    cobrapy stores annotations as ``str`` or ``list[str]`` depending on
    SBML loader and source; normalize to the canonical `;`-joined form
    that the MATLAB pipeline expects before validation.
    """
    if value is None or value == "":
        return ""
    if isinstance(value, str):
        return value
    if isinstance(value, (list, tuple)):
        return ";".join(str(v) for v in value if v)
    return str(value)


def get_ec_from_gem(
    model: "EcModel",
    ec_rxns: Optional[Iterable[str]] = None,
) -> None:
    """Populate ``model.ec.eccodes`` from per-reaction ``ec-code`` annotations.

    Ported from GECKO MATLAB:
    src/geckomat/get_enzyme_data/getECfromGEM.m.

    For each entry of ``model.ec.rxns``, looks up the corresponding
    cobra.Reaction (stripping a 4-char prefix first when the model is
    in gecko-light mode), reads ``reaction.annotation['ec-code']``,
    validates the resulting string against the canonical
    ``<ec>(;<ec>)*`` pattern (with ``<ec>`` being four dot-separated
    levels of either a non-negative integer or ``-``), and stores the
    result in ``model.ec.eccodes[i]``. Strings that fail validation
    are replaced with ``""`` and listed in a single
    ``logger.warning``.

    MATLAB-COMPAT: GECKO MATLAB stores EC codes in a top-level
    ``model.eccodes`` cell array. cobrapy stores them per-reaction in
    ``reaction.annotation['ec-code']``; geckopy reads from there.

    MATLAB-COMPAT: GECKO MATLAB returns ``(model, invalidEC,
    invalidECpos)``. geckopy drops both secondary outputs and emits a
    single warning instead.

    MATLAB-COMPAT: The MATLAB validation regex is broken: substituting
    ``$3`` on a valid ``"1.2.3.4"`` yields ``"3"``, which is non-empty
    and therefore flagged as invalid. As written, every non-empty EC
    string is discarded. geckopy implements the docstring intent.
    Tracked in ``docs/future_improvements.md``.

    MATLAB-COMPAT: ``ecRxns`` in MATLAB is a logical mask of length
    ``model.ec.rxns``. geckopy follows the cobrapy idiom of an
    iterable of reaction IDs (see ``model.remove_reactions``,
    ``flux_variability_analysis(reaction_list=...)``, etc.).

    Parameters
    ----------
    model
        EcModel with ``model.ec.rxns`` already allocated (typically by
        ``allocate_ec_for_catalyzed_reactions``). Mutated in place.
    ec_rxns
        Optional iterable of reaction IDs (each must appear in
        ``model.ec.rxns``) selecting which entries of ``ec.eccodes``
        to update. If ``None`` (default), every entry is rewritten.
        Unknown IDs raise ``ValueError``.
    """
    n = len(model.ec.rxns)
    if n == 0:
        return

    if ec_rxns is None:
        positions: list[int] = list(range(n))
    else:
        ec_rxns = list(ec_rxns)
        index_by_id = {rid: i for i, rid in enumerate(model.ec.rxns)}
        unknown = [rid for rid in ec_rxns if rid not in index_by_id]
        if unknown:
            preview = unknown[:5]
            raise ValueError(
                f"{len(unknown)} reaction ID(s) in ec_rxns are not present "
                f"in model.ec.rxns (examples: {preview})"
            )
        positions = [index_by_id[rid] for rid in ec_rxns]

    invalid: list[str] = []

    for i in positions:
        rxn_id = model.ec.rxns[i]
        if model.ec.gecko_light:
            rxn_id = rxn_id[4:]
        rxn = model.reactions.get_by_id(rxn_id)
        raw = _normalize_annotation(rxn.annotation.get("ec-code"))
        if not raw:
            model.ec.eccodes[i] = ""
            continue
        if _EC_MULTI_RE.match(raw):
            model.ec.eccodes[i] = raw
        else:
            invalid.append(raw)
            model.ec.eccodes[i] = ""

    if invalid:
        logger.warning(
            "get_ec_from_gem: skipped %d invalid EC string(s): %s",
            len(invalid),
            ", ".join(repr(s) for s in invalid),
        )
