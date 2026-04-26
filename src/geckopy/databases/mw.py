"""Compute molecular weight of a protein from its amino acid sequence.

Ported from GECKO MATLAB:
src/geckomat/get_enzyme_data/calculateMW.m.
"""
from __future__ import annotations

import logging

logger = logging.getLogger(__name__)


# Average mass of water (g/mol). Added once to account for the two
# H atoms and one OH that remain after peptide bond condensation.
# MATLAB-COMPAT: GECKO MATLAB uses the rounded value 18. geckopy uses
# the precise average mass for accuracy. The difference is ~0.015 Da
# per protein, negligible for ec.mw uses.
_WATER_MASS = 18.01528


# Average residue masses (g/mol). Taken from GECKO MATLAB
# `calculateMW.m`. Standard amino acids agree with widely-cited
# references (e.g., ExPASy ProtParam). Non-standard codes:
# - B = average(D, N).
# - J = same as L and I, which both happen to equal 113.16.
# - O (pyrrolysine) = 255.31.
# - U (selenocysteine) = 150.04.
# - Z = average(E, Q).
#
# X (any amino acid): geckopy uses the unweighted mean of the 20
# standard residues for self-consistency with this exact table.
# MATLAB-COMPAT: GECKO MATLAB uses 126.50 for X, an unsourced
# historical value. MW values for sequences containing X will differ
# from MATLAB by 7.11 Da per X.
_STANDARD_RESIDUE_MASSES: dict[str, float] = {
    "A": 71.08,
    "C": 103.14,
    "D": 115.09,
    "E": 129.11,
    "F": 147.17,
    "G": 57.05,
    "H": 137.14,
    "I": 113.16,
    "K": 128.17,
    "L": 113.16,
    "M": 131.20,
    "N": 114.10,
    "P": 97.12,
    "Q": 128.13,
    "R": 156.19,
    "S": 87.08,
    "T": 101.10,
    "V": 99.13,
    "W": 186.21,
    "Y": 163.17,
}

_X_MASS = sum(_STANDARD_RESIDUE_MASSES.values()) / len(_STANDARD_RESIDUE_MASSES)

_RESIDUE_MASSES: dict[str, float] = {
    **_STANDARD_RESIDUE_MASSES,
    "B": (115.09 + 114.10) / 2,  # D-or-N
    "J": 113.16,                  # L-or-I (both 113.16)
    "O": 255.31,                  # pyrrolysine
    "U": 150.04,                  # selenocysteine
    "X": _X_MASS,                 # any standard residue
    "Z": (129.11 + 128.13) / 2,  # E-or-Q
}


def calculate_mw(sequence: str) -> float:
    """Compute the molecular weight of a protein from its sequence.

    Ported from GECKO MATLAB:
    src/geckomat/get_enzyme_data/calculateMW.m.

    Each residue contributes its average mass; one water mass is added
    once (representing the H + OH at the peptide termini, equivalent to
    one un-eliminated water relative to N residues + N water masses).

    Recognized amino acid codes (case-INsensitive in geckopy):

        ACDEFGHIKLMNPQRSTVWY  - the 20 standard residues
        BJOZUX                - non-standard / ambiguous codes

    Characters not in the table (whitespace, digits, punctuation,
    lowercase letters that are not also recognized) are skipped with
    a warning.

    MATLAB-COMPAT: MATLAB is case-sensitive (lowercase letters silently
    skipped). geckopy uppercases first.

    MATLAB-COMPAT: MATLAB silently skips unknown characters. geckopy
    warns once per call listing the offending character set.

    Parameters
    ----------
    sequence
        Protein sequence as a string. May contain whitespace, digits,
        and other non-letter characters; those are ignored.

    Returns
    -------
    float
        Molecular weight in g/mol (Da). Returns the water mass
        (18.01528) for an empty/all-skipped sequence.
    """
    upper = sequence.upper() if sequence else ""
    mw = _WATER_MASS

    unknown_chars: set[str] = set()
    for ch in upper:
        if ch in _RESIDUE_MASSES:
            mw += _RESIDUE_MASSES[ch]
        elif not ch.isspace() and not ch.isdigit():
            unknown_chars.add(ch)

    if unknown_chars:
        logger.warning(
            "calculate_mw: skipped %d unknown character(s): %s",
            len(unknown_chars),
            ", ".join(repr(c) for c in sorted(unknown_chars)),
        )

    return mw
