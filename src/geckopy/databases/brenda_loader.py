"""Parse BRENDA TSV dumps into kcat and SA-derived-kcat tables.

Ported from GECKO MATLAB:
src/geckomat/get_enzyme_data/loadBRENDAdata.m.

Three tab-delimited files are read from the BRENDA folder:

    max_kcat.tsv   kcat values
    max_sa.tsv     specific activities
    max_mw.tsv     molecular weights

The files are produced by the ``geckopy brenda-refresh`` CLI. Each
has 5 tab-delimited columns: EC number, substrate (``*`` for SA/MW),
organism, numeric value, references (semicolon-joined PMIDs or ``*``).
A ``#`` header line carries the BRENDA release version.
"""
from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)

_KCAT_COLUMNS = ["ec_code", "substrate", "organism", "kcat"]
_SA_COLUMNS = ["ec_code", "organism", "kcat", "mw"]


@dataclass
class BrendaData:
    """Loaded BRENDA tables.

    Attributes
    ----------
    kcat
        DataFrame with columns ``ec_code`` (str), ``substrate`` (str),
        ``organism`` (str), ``kcat`` (float, 1/s).
    sa
        DataFrame with columns ``ec_code`` (str), ``organism`` (str),
        ``kcat`` (float, 1/s, computed as SA * MW for matching
        EC + organism), ``mw`` (float, g/mmol).
    """
    kcat: pd.DataFrame = field(
        default_factory=lambda: pd.DataFrame(columns=_KCAT_COLUMNS)
    )
    sa: pd.DataFrame = field(
        default_factory=lambda: pd.DataFrame(columns=_SA_COLUMNS)
    )


def load_brenda_data(folder: str | Path) -> BrendaData:
    """Load BRENDA kcat, SA, and MW dumps from a folder.

    Ported from GECKO MATLAB:
    src/geckomat/get_enzyme_data/loadBRENDAdata.m.

    Reads ``max_kcat.tsv``, ``max_sa.tsv``, ``max_mw.tsv`` from
    ``folder``, applies the GECKO unit conversions, and joins
    SA + MW on EC + organism (case-insensitive) to produce a
    derived kcat table.

    Unit conventions:

    * KCAT: file value is in 1/s already; no scaling.
    * SA:   file value in ``[umol/min/mg]`` is multiplied by ``1/60``
            to give ``[mmol/s/g]``.
    * MW:   file value in ``[g/mol]`` is multiplied by ``1/1000`` to
            give ``[g/mmol]``.
    * Derived kcat in the SA table is ``SA * MW = [1/s]``.

    MATLAB-COMPAT: GECKO MATLAB takes a ``modelAdapter`` and resolves
    the path via ``modelAdapter.getBrendaDBFolder()``. geckopy takes
    a folder path directly; the caller resolves the path.

    Malformed lines (wrong number of fields, non-numeric column 4)
    are skipped with a ``logger.warning``; the rest of the file is
    parsed normally.

    Parameters
    ----------
    folder
        Directory containing all three BRENDA dump files.

    Returns
    -------
    BrendaData
        Both DataFrames may be empty if the source files were empty
        or contained only malformed lines, or (for ``sa``) if no SA
        rows had a matching MW row.

    Raises
    ------
    FileNotFoundError
        If any of the three expected files is missing.
    """
    folder = Path(folder)
    kcat_path = folder / "max_kcat.tsv"
    sa_path = folder / "max_sa.tsv"
    mw_path = folder / "max_mw.tsv"

    for p in (kcat_path, sa_path, mw_path):
        if not p.is_file():
            raise FileNotFoundError(f"Expected BRENDA file not found: {p}")

    kcat_raw = _load_table(kcat_path, 1.0)
    sa_raw = _load_table(sa_path, 1.0 / 60.0)
    mw_raw = _load_table(mw_path, 1.0 / 1000.0)

    kcat_table = kcat_raw.rename(columns={"value": "kcat"})[_KCAT_COLUMNS]
    sa_table = _join_sa_with_mw(sa_raw, mw_raw)

    return BrendaData(kcat=kcat_table, sa=sa_table)


def _load_table(path: Path, value_scale: float) -> pd.DataFrame:
    """Parse a 5-column tab-delimited BRENDA dump.

    Returns a DataFrame with columns ``ec_code``, ``substrate``,
    ``organism``, ``value``. The 5th column (references) is dropped.
    Lines beginning with ``#`` are treated as comments and skipped.
    """
    rows: list[tuple[str, str, str, float]] = []
    invalid_lines = 0

    with open(path, "r", encoding="utf-8") as f:
        for line_no, line in enumerate(f, start=1):
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 5:
                invalid_lines += 1
                logger.warning(
                    "%s:%d: expected 5 tab-delimited fields, got %d; skipping.",
                    path.name, line_no, len(parts),
                )
                continue
            ec, substrate, organism, value_str, _ = parts[:5]
            try:
                value = float(value_str) * value_scale
            except ValueError:
                invalid_lines += 1
                logger.warning(
                    "%s:%d: column 4 %r is not numeric; skipping.",
                    path.name, line_no, value_str,
                )
                continue
            rows.append((ec, substrate, organism, value))

    if invalid_lines:
        logger.warning(
            "%s: skipped %d malformed line(s).", path.name, invalid_lines,
        )

    return pd.DataFrame(
        rows,
        columns=["ec_code", "substrate", "organism", "value"],
    )


def _join_sa_with_mw(sa: pd.DataFrame, mw: pd.DataFrame) -> pd.DataFrame:
    """Join SA + MW on (EC, organism) case-insensitive, first MW match wins.

    For each SA row, find the first MW row with the same EC and
    organism (case-insensitive). If found, emit a row with the
    SA-derived kcat. SA rows without a matching MW row are dropped,
    matching MATLAB.
    """
    if sa.empty or mw.empty:
        return pd.DataFrame(columns=_SA_COLUMNS)

    mw_lookup: dict[tuple[str, str], float] = {}
    for ec, organism, mw_value in zip(
        mw["ec_code"].str.upper(),
        mw["organism"].str.upper(),
        mw["value"],
    ):
        mw_lookup.setdefault((ec, organism), mw_value)

    out_rows: list[dict[str, object]] = []
    for ec, organism, sa_value in zip(
        sa["ec_code"], sa["organism"], sa["value"]
    ):
        mw_value = mw_lookup.get((ec.upper(), organism.upper()))
        if mw_value is None:
            continue
        out_rows.append({
            "ec_code": ec,
            "organism": organism,
            "kcat": sa_value * mw_value,
            "mw": mw_value,
        })

    if not out_rows:
        return pd.DataFrame(columns=_SA_COLUMNS)
    return pd.DataFrame(out_rows, columns=_SA_COLUMNS)
