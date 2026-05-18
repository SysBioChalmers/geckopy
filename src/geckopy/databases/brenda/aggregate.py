"""Aggregate parsed BRENDA rows to max-per-key and write the three TSVs.

Rows sharing ``(kind, ec, substrate, organism)`` collapse to one
output row holding the maximum value and the sorted-deduplicated
union of all PMIDs that backed the merged inputs.
"""
from __future__ import annotations

import datetime as dt
import logging
from pathlib import Path
from typing import Iterable

from .parse import Kind, Row

logger = logging.getLogger(__name__)


_OUTPUT_FILES: dict[Kind, tuple[str, str]] = {
    "kcat": ("max_kcat.tsv", "kcat in 1/s"),
    "sa":   ("max_sa.tsv",   "specific activity in umol/min/mg"),
    "mw":   ("max_mw.tsv",   "molecular weight in g/mol"),
}

_HEADER_FMT = (
    "# BRENDA release {release} generated {date} - CC BY 4.0 - {desc}\n"
)


def aggregate_and_write(
    rows: Iterable[Row],
    out_dir: str | Path,
    *,
    release: str,
    date: str | None = None,
) -> dict[Kind, Path]:
    """Group rows to one max-row per (ec, substrate, organism); write three TSVs.

    Parameters
    ----------
    rows
        Iterable of ``Row`` from ``parse_brenda_json``.
    out_dir
        Output directory for ``max_kcat.tsv``, ``max_sa.tsv``,
        ``max_mw.tsv``. Created if missing.
    release
        BRENDA release string, e.g. ``"2026.1"``. Written into the
        header line of each TSV.
    date
        ISO date string for the header. Defaults to today.

    Returns
    -------
    dict
        Mapping of ``Kind`` to written path.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    if date is None:
        date = dt.date.today().isoformat()

    aggregated: dict[tuple[Kind, str, str, str], tuple[float, set[str]]] = {}
    for r in rows:
        key = (r.kind, r.ec, r.substrate, r.organism)
        existing = aggregated.get(key)
        if existing is None:
            aggregated[key] = (r.value, set(r.references))
        else:
            aggregated[key] = (
                max(existing[0], r.value),
                existing[1] | set(r.references),
            )

    paths: dict[Kind, Path] = {}
    for kind, (basename, desc) in _OUTPUT_FILES.items():
        path = out_dir / basename
        paths[kind] = path
        keys = sorted(
            (k for k in aggregated if k[0] == kind),
            key=lambda k: (k[1], k[2], k[3]),
        )
        with path.open("w", encoding="utf-8", newline="\n") as fh:
            fh.write(_HEADER_FMT.format(release=release, date=date, desc=desc))
            for k in keys:
                value, refs = aggregated[k]
                refs_field = ";".join(sorted(refs)) if refs else "*"
                fh.write(
                    f"{k[1]}\t{k[2]}\t{k[3]}\t{_fmt_value(value)}\t{refs_field}\n"
                )
        logger.info("wrote %d rows to %s", len(keys), path)
    return paths


def _fmt_value(v: float) -> str:
    # Python ``str(float)`` always includes a decimal point or exponent,
    # so ``598.0`` and ``2.5`` round-trip identically to the 2018 file.
    return str(v)
