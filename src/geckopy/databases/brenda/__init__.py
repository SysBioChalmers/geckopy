"""Refresh the BRENDA-derived max_kcat / max_sa / max_mw TSVs from the BRENDA bulk JSON.

See ``docs/brenda_refresh_plan.md`` for the design.
"""
from __future__ import annotations

from .parse import Row, parse_brenda_json

__all__ = [
    "Row",
    "parse_brenda_json",
]
