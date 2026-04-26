"""Format raw EC numbers from a database into a canonical string.

Ported from GECKO MATLAB:
src/geckomat/get_enzyme_data/getECstring.m.
"""
from __future__ import annotations

import logging
import re

logger = logging.getLogger(__name__)


# Valid EC token: four dot-separated levels, each a non-negative integer
# or `-` (the IUBMB placeholder for an unspecified level). Examples:
# `1.1.1.1`, `1.1.1.-`, `1.-.-.-`. The hierarchical IUBMB constraint
# that a `-` at level N implies `-` at every deeper level is NOT
# enforced here; well-shaped but semantically odd tokens like `1.-.1.1`
# will pass.
_EC_TOKEN_RE = re.compile(r"^(?:\d+|-)\.(?:\d+|-)\.(?:\d+|-)\.(?:\d+|-)$")


def get_ec_string(ec_numbers: str) -> str:
    """Format raw EC numbers into a single space-separated, EC-prefixed string.

    Ported from GECKO MATLAB:
    src/geckomat/get_enzyme_data/getECstring.m.

    Splits ``ec_numbers`` on whitespace (consecutive whitespace and
    leading/trailing whitespace are dropped). Each token is stripped of
    any ``;`` characters and any leading ``EC`` / ``ec`` prefix, then
    validated against the canonical EC pattern
    ``\\d+\\.\\d+\\.\\d+\\.\\d+`` (with ``-`` allowed in any level).
    Valid tokens are joined with single spaces, each prefixed with
    ``EC``. Tokens that fail validation are logged via
    ``logger.warning`` and skipped.

    MATLAB-COMPAT: GECKO MATLAB is accumulator-based (the first
    argument is a string the function appends to). geckopy is a pure
    function; callers compose with ``" ".join(...)`` if accumulator
    behaviour is needed.

    MATLAB-COMPAT: GECKO MATLAB returns the bare string ``"EC"`` for
    empty input (a quirk of ``strsplit("", " ")`` returning ``{''}``).
    geckopy returns ``""``.

    MATLAB-COMPAT: GECKO MATLAB does no validation. It will produce
    ``ECEC1.1.1.1`` for already-prefixed input and ``ECnotanec`` for
    junk. geckopy strips a leading ``EC`` (case-insensitive) before
    re-prefixing, and skips invalid tokens with a warning.

    Parameters
    ----------
    ec_numbers
        Whitespace-separated EC tokens as written in a UniProt or
        BRENDA-style database export. May contain ``;`` characters
        (stripped from each token). May already include ``EC`` prefixes
        (stripped before re-prefixing). Empty or all-skipped inputs
        return ``""``.

    Returns
    -------
    str
        Space-separated, ``EC``-prefixed canonical string, e.g.
        ``"EC1.1.1.1 EC2.7.7.7"``. Empty if no tokens passed
        validation.
    """
    if not ec_numbers:
        return ""

    valid_tokens: list[str] = []
    invalid_tokens: list[str] = []

    for raw in ec_numbers.split():
        token = raw.replace(";", "")
        if len(token) >= 2 and token[:2].upper() == "EC":
            token = token[2:]
        if not token or not _EC_TOKEN_RE.match(token):
            invalid_tokens.append(raw)
            continue
        valid_tokens.append(token)

    if invalid_tokens:
        logger.warning(
            "get_ec_string: skipped %d invalid EC token(s): %s",
            len(invalid_tokens),
            ", ".join(repr(t) for t in invalid_tokens),
        )

    return " ".join(f"EC{t}" for t in valid_tokens)
