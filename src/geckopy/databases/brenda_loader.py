"""Parse BRENDA TSV dumps into kcat and SA-derived-kcat tables.

Ported from GECKO MATLAB:
src/geckomat/get_enzyme_data/loadBRENDAdata.m.

Three tab-delimited files are read from the BRENDA folder:

    kcat.tsv   kcat values (wide: one row per triple, both max and
               median across the raw measurements that fed into it)
    sa.tsv     specific activities (same wide shape)
    mw.tsv     molecular weights (single value per (ec, organism))

The files are produced by the ``geckopy brenda-refresh`` CLI. kcat and
SA have seven tab-delimited columns: EC number, substrate (``*`` for
SA), organism, value-max, value-median, n (number of raw measurements
aggregated), references (semicolon-joined PMIDs or ``*``). MW has six
columns (single value, no aggregation choice). A ``#`` header line
carries the BRENDA release version; the line after it is the TSV
column header (also skipped).
"""
from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal

import pandas as pd

logger = logging.getLogger(__name__)

_KCAT_COLUMNS = ["ec_code", "substrate", "organism", "kcat", "n"]
_SA_COLUMNS = ["ec_code", "organism", "kcat", "mw", "n"]


@dataclass
class BrendaData:
    """Loaded BRENDA tables.

    Each kcat / SA table is exposed as separate ``_max`` and ``_median``
    DataFrames built once at load time. Both views share the same row
    set (one per (ec, substrate, organism) triple); they differ only in
    which aggregation of the raw measurements fills the ``kcat``
    column. Choose between views via :meth:`kcat_for` / :meth:`sa_for`.

    Attributes
    ----------
    kcat_max, kcat_median
        DataFrames with columns ``ec_code`` (str), ``substrate`` (str),
        ``organism`` (str), ``kcat`` (float, 1/s), ``n`` (Int64).
    sa_max, sa_median
        Same shape, joined with MW to give an SA-derived kcat. Columns
        ``ec_code``, ``organism``, ``kcat`` (float, 1/s), ``mw`` (float,
        g/mmol), ``n`` (Int64).
    """
    kcat_max: pd.DataFrame = field(
        default_factory=lambda: pd.DataFrame(columns=_KCAT_COLUMNS)
    )
    kcat_median: pd.DataFrame = field(
        default_factory=lambda: pd.DataFrame(columns=_KCAT_COLUMNS)
    )
    sa_max: pd.DataFrame = field(
        default_factory=lambda: pd.DataFrame(columns=_SA_COLUMNS)
    )
    sa_median: pd.DataFrame = field(
        default_factory=lambda: pd.DataFrame(columns=_SA_COLUMNS)
    )

    def kcat_for(self, aggregation: Literal["max", "median"]) -> pd.DataFrame:
        """Return the kcat view for the requested aggregation."""
        if aggregation == "median":
            return self.kcat_median
        if aggregation == "max":
            return self.kcat_max
        raise ValueError(
            f"aggregation must be 'max' or 'median', got {aggregation!r}"
        )

    def sa_for(self, aggregation: Literal["max", "median"]) -> pd.DataFrame:
        """Return the SA-derived-kcat view for the requested aggregation."""
        if aggregation == "median":
            return self.sa_median
        if aggregation == "max":
            return self.sa_max
        raise ValueError(
            f"aggregation must be 'max' or 'median', got {aggregation!r}"
        )


def load_brenda_data(folder: str | Path) -> BrendaData:
    """Load BRENDA kcat, SA, and MW dumps from a folder.

    Ported from GECKO MATLAB:
    src/geckomat/get_enzyme_data/loadBRENDAdata.m.

    Reads ``kcat.tsv``, ``sa.tsv``, ``mw.tsv`` from ``folder``, splits
    the kcat and SA wide tables into per-aggregation views, applies the
    GECKO unit conversions, and joins SA + MW on EC + organism
    (case-insensitive) to produce derived kcat tables.

    Unit conventions:

    * KCAT: file values are in 1/s already; no scaling.
    * SA:   file values in ``[umol/min/mg]`` are multiplied by ``1/60``
            to give ``[mmol/s/g]``.
    * MW:   file value in ``[g/mol]`` is multiplied by ``1/1000`` to
            give ``[g/mmol]``.
    * Derived kcat in the SA table is ``SA * MW = [1/s]``.

    Both kcat and SA are loaded as separate max and median views over
    the raw measurements aggregated per (ec, substrate, organism)
    triple, so a project can opt into median aggregation via
    ``adapter.params.kcat_aggregate_brenda`` without regenerating the
    snapshot. See ``docs/kcat_aggregation.md`` for the rationale.

    Malformed lines (wrong number of fields, non-numeric value column)
    are skipped with a ``logger.warning``; the rest of the file is
    parsed normally.

    Parameters
    ----------
    folder
        Directory containing all three BRENDA dump files.

    Returns
    -------
    BrendaData
        Both kcat views may be empty if the source files were empty or
        contained only malformed lines; both SA views may be empty if
        no SA row matched any MW row.

    Raises
    ------
    FileNotFoundError
        If any of the three expected files is missing.
    """
    folder = Path(folder)
    kcat_path = folder / "kcat.tsv"
    sa_path = folder / "sa.tsv"
    mw_path = folder / "mw.tsv"

    for p in (kcat_path, sa_path, mw_path):
        if not p.is_file():
            raise FileNotFoundError(f"Expected BRENDA file not found: {p}")

    kcat_wide = _load_wide_table(
        kcat_path, value_columns=("kcat_max", "kcat_median"),
        value_scale=1.0,
    )
    sa_wide = _load_wide_table(
        sa_path, value_columns=("sa_max", "sa_median"),
        value_scale=1.0 / 60.0,
    )
    mw_raw = _load_mw_table(mw_path, value_scale=1.0 / 1000.0)

    kcat_max = _view_from_wide(kcat_wide, "kcat_max", _KCAT_COLUMNS)
    kcat_median = _view_from_wide(kcat_wide, "kcat_median", _KCAT_COLUMNS)
    sa_max = _join_sa_with_mw(
        _sa_view_from_wide(sa_wide, "sa_max"), mw_raw,
    )
    sa_median = _join_sa_with_mw(
        _sa_view_from_wide(sa_wide, "sa_median"), mw_raw,
    )

    return BrendaData(
        kcat_max=kcat_max,
        kcat_median=kcat_median,
        sa_max=sa_max,
        sa_median=sa_median,
    )


# --------------------------------------------------------------------------- #
# File parsing
# --------------------------------------------------------------------------- #

# Wide-format file layout (kcat / SA): 7 columns
# 0=ec_code, 1=substrate, 2=organism, 3=value_max, 4=value_median,
# 5=n, 6=refs
_WIDE_NCOLS = 7
# MW file layout: 6 columns
# 0=ec_code, 1=substrate, 2=organism, 3=value, 4=n, 5=refs
_MW_NCOLS = 6


def _load_wide_table(
    path: Path,
    *,
    value_columns: tuple[str, str],
    value_scale: float,
) -> pd.DataFrame:
    """Parse a 7-column kcat/SA TSV (wide format) into a DataFrame.

    Returns a DataFrame with columns ``ec_code``, ``substrate``,
    ``organism``, ``<value_columns[0]>``, ``<value_columns[1]>``, ``n``.
    The references column is dropped. Comment lines (``#``-prefixed)
    and the column-header row are skipped.
    """
    max_col, med_col = value_columns
    rows: list[tuple[str, str, str, float, float, int]] = []
    invalid = 0
    with open(path, "r", encoding="utf-8") as f:
        for line_no, line in enumerate(f, start=1):
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            # Skip the TSV column-header row (string in a numeric column).
            if line_no <= 3 and parts[0] == "ec_code":
                continue
            if len(parts) < _WIDE_NCOLS:
                invalid += 1
                logger.warning(
                    "%s:%d: expected %d tab-delimited fields, got %d; skipping.",
                    path.name, line_no, _WIDE_NCOLS, len(parts),
                )
                continue
            ec, substrate, organism, vmax_s, vmed_s, n_s, _ = parts[:_WIDE_NCOLS]
            try:
                vmax = float(vmax_s) * value_scale
                vmed = float(vmed_s) * value_scale
                n = int(n_s)
            except ValueError:
                invalid += 1
                logger.warning(
                    "%s:%d: non-numeric value(s) (%r, %r, %r); skipping.",
                    path.name, line_no, vmax_s, vmed_s, n_s,
                )
                continue
            rows.append((ec, substrate, organism, vmax, vmed, n))

    if invalid:
        logger.warning("%s: skipped %d malformed line(s).", path.name, invalid)

    return pd.DataFrame(
        rows,
        columns=["ec_code", "substrate", "organism", max_col, med_col, "n"],
    )


def _load_mw_table(path: Path, *, value_scale: float) -> pd.DataFrame:
    """Parse the 6-column MW TSV (no aggregation choice).

    Returns a DataFrame with columns ``ec_code``, ``organism``,
    ``value``. ``substrate`` is always ``*`` and is dropped on the way
    out; ``n`` and the references column are dropped too.
    """
    rows: list[tuple[str, str, float]] = []
    invalid = 0
    with open(path, "r", encoding="utf-8") as f:
        for line_no, line in enumerate(f, start=1):
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if line_no <= 3 and parts[0] == "ec_code":
                continue
            if len(parts) < _MW_NCOLS:
                invalid += 1
                logger.warning(
                    "%s:%d: expected %d tab-delimited fields, got %d; skipping.",
                    path.name, line_no, _MW_NCOLS, len(parts),
                )
                continue
            ec, _, organism, value_str, _, _ = parts[:_MW_NCOLS]
            try:
                value = float(value_str) * value_scale
            except ValueError:
                invalid += 1
                logger.warning(
                    "%s:%d: column 4 %r is not numeric; skipping.",
                    path.name, line_no, value_str,
                )
                continue
            rows.append((ec, organism, value))

    if invalid:
        logger.warning("%s: skipped %d malformed line(s).", path.name, invalid)

    return pd.DataFrame(rows, columns=["ec_code", "organism", "value"])


def _view_from_wide(
    wide: pd.DataFrame, value_col: str, out_columns: list[str],
) -> pd.DataFrame:
    """Project the wide kcat table down to one of its value columns.

    Returns a DataFrame with the public schema
    ``[ec_code, substrate, organism, kcat, n]``.
    """
    if wide.empty:
        return pd.DataFrame(columns=out_columns)
    sub = wide[["ec_code", "substrate", "organism", value_col, "n"]].copy()
    sub = sub.rename(columns={value_col: "kcat"})
    return sub.reset_index(drop=True)


def _sa_view_from_wide(wide: pd.DataFrame, value_col: str) -> pd.DataFrame:
    """Project the wide SA table to one aggregation, dropping substrate.

    Returns a DataFrame with columns ``ec_code``, ``organism``, ``sa``,
    ``n``. The MW join helper attaches MW + computes ``kcat = sa * mw``.
    """
    if wide.empty:
        return pd.DataFrame(columns=["ec_code", "organism", "sa", "n"])
    sub = wide[["ec_code", "organism", value_col, "n"]].copy()
    sub = sub.rename(columns={value_col: "sa"})
    return sub.reset_index(drop=True)


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
    for ec, organism, sa_value, n in zip(
        sa["ec_code"], sa["organism"], sa["sa"], sa["n"],
    ):
        mw_value = mw_lookup.get((ec.upper(), organism.upper()))
        if mw_value is None:
            continue
        out_rows.append({
            "ec_code": ec,
            "organism": organism,
            "kcat": sa_value * mw_value,
            "mw": mw_value,
            "n": int(n),
        })

    if not out_rows:
        return pd.DataFrame(columns=_SA_COLUMNS)
    return pd.DataFrame(out_rows, columns=_SA_COLUMNS)
