"""Shared pytest configuration.

Selects Gurobi as cobra's default LP solver for the whole test session
when it is importable and licensed, so every FBA a test triggers --
including the ones inside worker processes, which inherit the solver
interface with the pickled model -- goes through the same solver.

Falls back silently to whatever cobra would have picked (typically
GLPK) when gurobipy is absent or its license is unavailable, so the
suite still runs on a machine without Gurobi. Library code stays
solver-agnostic: geckopy reads ``cobra.Configuration()`` for this the
same way it does for ``processes``.
"""
from __future__ import annotations

import cobra


def pytest_configure():
    try:
        import gurobipy  # noqa: F401
    except ImportError:
        return
    try:
        cobra.Configuration().solver = "gurobi"
    except Exception as exc:  # unlicensed, or size-limited license
        print(f"conftest: Gurobi present but unusable ({exc}); leaving cobra's default solver.")
