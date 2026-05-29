"""Expand reactions with isozymes into one reaction per isozyme.

Ported from GECKO MATLAB: src/geckomat/change_model/makeEcModel.m
(stage 5). Equivalent to RAVEN: expandModel.m, reimplemented using
cobrapy's GPR AST rather than RAVEN's string manipulation.

Only applies to full ecModels (geckoLight skips this stage).
"""
from __future__ import annotations

import ast
from typing import TYPE_CHECKING

import cobra
from cobra.core.gene import GPR

if TYPE_CHECKING:
    pass


def _gpr_to_dnf(gpr: GPR) -> list[list[str]]:
    """Convert a GPR to disjunctive normal form (list of AND-clauses).

    An empty GPR yields an empty list. A single clause (no OR anywhere)
    yields a list of length 1. OR-of-ANDs yields one sublist per
    disjunct, each containing the gene names ANDed together.

    Handles distributivity: ``g1 and (g2 or g3)`` becomes
    ``[[g1, g2], [g1, g3]]``.
    """
    if gpr is None or gpr.body is None:
        return []
    return _node_to_dnf(gpr.body)


def _node_to_dnf(node) -> list[list[str]]:
    """Recursive helper. Returns DNF as list of AND-clauses."""
    if isinstance(node, ast.Name):
        return [[node.id]]
    if isinstance(node, ast.BoolOp):
        if isinstance(node.op, ast.Or):
            result: list[list[str]] = []
            for child in node.values:
                result.extend(_node_to_dnf(child))
            return result
        if isinstance(node.op, ast.And):
            clauses: list[list[str]] = [[]]
            for child in node.values:
                child_dnf = _node_to_dnf(child)
                new_clauses: list[list[str]] = []
                for existing in clauses:
                    for extra in child_dnf:
                        new_clauses.append(existing + extra)
                clauses = new_clauses
            return clauses
    raise ValueError(f"Unexpected GPR node type: {type(node).__name__}")


def expand_model(model: "cobra.Model") -> list[str]:
    """Split reactions with isozymes (OR in GPR) into one reaction per isozyme.

    Ported from GECKO MATLAB: src/geckomat/change_model/makeEcModel.m
    (stage 5). Equivalent to RAVEN: expandModel.m, reimplemented using
    cobrapy.

    For each reaction whose GPR contains at least one OR, the reaction
    is removed and replaced by one copy per disjunctive clause. The new
    reactions get ID suffix ``_EXP_1``, ``_EXP_2``, etc. All other
    fields (stoichiometry, bounds, name, subsystem) are copied verbatim;
    only the GPR is simplified to the single AND-clause for that
    isozyme.

    Reactions with no GPR, or with a GPR that has no OR, are left
    untouched.

    Parameters
    ----------
    model
        A cobra.Model, mutated in place.

    Returns
    -------
    list of str
        Sorted IDs of newly added expanded reactions (those with
        ``_EXP_N`` suffixes). The original reactions that were split
        are no longer in the model.
    """
    expansions: list[tuple[cobra.Reaction, list[list[str]]]] = []

    for rxn in model.reactions:
        if not rxn.gene_reaction_rule:
            continue
        clauses = _gpr_to_dnf(rxn.gpr)
        if len(clauses) <= 1:
            continue
        expansions.append((rxn, clauses))

    added_ids: list[str] = []
    for original_rxn, clauses in expansions:
        new_rxns: list[cobra.Reaction] = []
        for i, clause in enumerate(clauses, start=1):
            new_rxn = cobra.Reaction(
                id=f"{original_rxn.id}_EXP_{i}",
                name=original_rxn.name,
            )
            new_rxn.lower_bound = original_rxn.lower_bound
            new_rxn.upper_bound = original_rxn.upper_bound
            new_rxn.add_metabolites(
                {m: c for m, c in original_rxn.metabolites.items()}
            )
            new_rxn.subsystem = original_rxn.subsystem
            new_rxn.gene_reaction_rule = " and ".join(clause)
            new_rxns.append(new_rxn)

        model.remove_reactions([original_rxn])
        model.add_reactions(new_rxns)
        added_ids.extend(r.id for r in new_rxns)

    return sorted(added_ids)
