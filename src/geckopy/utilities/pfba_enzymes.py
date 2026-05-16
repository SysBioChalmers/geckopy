"""Parsimonious FBA minimising enzyme usage.

Ported from the legacy geckopy package described in Carrasco et al.
(2023, https://doi.org/10.1128/spectrum.01705-23), file
geckopy/flux_analysis.py:342-386 (pfba_protein), adapted to the new
substrate (``usage_prot_<id>`` reactions instead of a Protein DictList).
"""
from __future__ import annotations

from typing import TYPE_CHECKING, Union

import cobra
from optlang.symbolics import Zero

from ..ec_model.constants import USAGE_PREFIX

if TYPE_CHECKING:
    from ..ec_model.ec_model import EcModel


_PFBA_OBJ_NAME = "_pfba_enzymes_objective"


def pfba_enzymes(
    model: "EcModel",
    objective: Union[cobra.Reaction, str, None] = None,
    *,
    fraction_of_optimum: float = 1.0,
) -> cobra.Solution:
    """Parsimonious FBA minimising the L1 norm of ``usage_prot_*`` fluxes.

    Sets the objective (if given), fixes it at
    ``fraction_of_optimum * optimal_value``, then minimises
    ``sum(usage_prot_*.forward_variable)``. Usage reactions are
    forward-only in the GECKO 3 layout (``lb=0``); if a user has
    flipped one to allow reverse flux, its ``reverse_variable`` is
    included automatically.

    Returns a standard ``cobra.Solution`` over all reactions.

    Ported from the legacy geckopy package (Carrasco et al., 2023,
    https://doi.org/10.1128/spectrum.01705-23),
    geckopy/flux_analysis.py:342-386.
    """
    if model.ec.gecko_light:
        raise NotImplementedError(
            "pfba_enzymes requires usage reactions (not gecko-light)"
        )

    with model as m:
        if objective is not None:
            m.objective = objective
        if m.objective.name == _PFBA_OBJ_NAME:
            raise ValueError(
                f"Objective is already {_PFBA_OBJ_NAME}; "
                "pfba_enzymes appears to be nested."
            )
        cobra.util.fix_objective_as_constraint(m, fraction=fraction_of_optimum)

        usage_rxns = [r for r in m.reactions if r.id.startswith(USAGE_PREFIX)]
        variables = []
        for r in usage_rxns:
            variables.append(r.forward_variable)
            if r.lower_bound < 0:
                variables.append(r.reverse_variable)

        m.objective = m.problem.Objective(
            Zero, direction="min", sloppy=True, name=_PFBA_OBJ_NAME,
        )
        m.objective.set_linear_coefficients({v: 1.0 for v in variables})
        m.slim_optimize(error_value=None)
        solution = cobra.core.solution.get_solution(m)
    return solution
