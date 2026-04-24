"""Structural preprocessing of the GEM before enzyme extension.

Corresponds to stages 1 to 4 of makeEcModel in GECKO MATLAB. These
stages do not touch the `ec` substructure; they only reshape the
underlying model so subsequent stages operate on a clean irreversible form.
"""
from __future__ import annotations
from typing import TYPE_CHECKING

import cobra

if TYPE_CHECKING:
    from ...adapter import ModelAdapter


def convert_to_irreversible(model: "cobra.Model") -> list[str]:
    """Split non-exchange reversible reactions into forward + reverse pair.

    Corresponds to stage 4 of GECKO MATLAB `makeEcModel`. Implicitly
    covers stage 3 (`rev` field computation) since that field is only
    consumed here and nowhere else in GECKO.

    For each non-exchange reaction with ``lb < 0``:

    - The original reaction is kept as the forward direction. Its
      lower bound is clamped to 0.
    - A new reaction with the same ID plus a ``_REV`` suffix is added,
      representing the reverse direction. Its stoichiometry is the
      negation of the original, its bounds are ``(0, -original_lb)``,
      and it inherits the name (with " (reversible)" appended) and the
      gene-protein rule of the original.

    Exchange reactions (reactions with exactly one participating
    metabolite) are never split, regardless of their bounds, to match
    the MATLAB behavior where exchange reactions are explicitly
    excluded from `convertToIrrev`.

    Parameters
    ----------
    model
        A cobra.Model, mutated in place.

    Returns
    -------
    list of str
        Sorted IDs of newly added reverse reactions (the ones ending in
        ``_REV``). The forward reactions retain their original IDs.
    """
    reverse_rxns_to_add: list[cobra.Reaction] = []
    forward_updates: list[tuple[cobra.Reaction, float]] = []

    for rxn in model.reactions:
        if rxn.boundary:
            continue
        if rxn.lower_bound >= 0:
            continue

        original_lb = rxn.lower_bound

        rev_rxn = cobra.Reaction(
            id=f"{rxn.id}_REV",
            name=(f"{rxn.name} (reversible)" if rxn.name else f"{rxn.id}_REV"),
        )
        rev_rxn.lower_bound = 0.0
        rev_rxn.upper_bound = -original_lb
        rev_rxn.add_metabolites({m: -c for m, c in rxn.metabolites.items()})
        rev_rxn.gene_reaction_rule = rxn.gene_reaction_rule

        reverse_rxns_to_add.append(rev_rxn)
        forward_updates.append((rxn, original_lb))

    for rxn, _ in forward_updates:
        rxn.lower_bound = 0.0

    if reverse_rxns_to_add:
        model.add_reactions(reverse_rxns_to_add)

    return sorted(r.id for r in reverse_rxns_to_add)


def remove_pseudoreaction_gprs(
    model: "cobra.Model", adapter: "ModelAdapter"
) -> list[str]:
    """Clear gene-protein-reaction rules from pseudoreactions.

    Corresponds to stage 1 of GECKO MATLAB `makeEcModel`. A reaction is
    treated as a pseudoreaction if either:

    - its name contains the substring ``pseudoreaction`` (case-sensitive,
      matching MATLAB's default `contains` behavior); or
    - its ID is listed in the first column of ``data/pseudoRxns.tsv``
      inside the adapter's project folder, if that file exists.

    Parameters
    ----------
    model
        A cobra.Model that will be mutated in place: the gene reaction
        rule of each matched reaction is cleared (empty string). This
        automatically updates the reaction's `.genes` set via cobrapy's
        GPR handling.
    adapter
        Provides the project path used to locate the optional TSV.

    Returns
    -------
    list of str
        Sorted IDs of reactions whose GPR was cleared. Intended for
        logging and test assertions.
    """
    pseudo_ids: set[str] = set()

    # By reaction name (MATLAB: contains(model.rxnNames,'pseudoreaction'))
    for rxn in model.reactions:
        if rxn.name and "pseudoreaction" in rxn.name:
            pseudo_ids.add(rxn.id)

    # By optional TSV (first column = reaction ID)
    tsv_path = adapter.params.path / "data" / "pseudoRxns.tsv"
    if tsv_path.is_file():
        with open(tsv_path, "r", encoding="utf-8") as f:
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if parts and parts[0]:
                    pseudo_ids.add(parts[0])

    # Clear GPRs for matched reactions that exist in the model.
    model_rxn_ids = {r.id for r in model.reactions}
    cleared: list[str] = []
    for rxn_id in pseudo_ids & model_rxn_ids:
        rxn = model.reactions.get_by_id(rxn_id)
        rxn.gene_reaction_rule = ""
        cleared.append(rxn_id)

    return sorted(cleared)

def invert_backwards_only_reactions(model: "cobra.Model") -> list[str]:
    """Flip reactions that are constrained to only carry negative flux.

    Corresponds to stage 2 of GECKO MATLAB `makeEcModel`. A reaction
    with ``lb < 0`` and ``ub == 0`` is physically identical to a
    forward-only reaction with inverted stoichiometry and bounds
    ``[0, -lb]``. Rewriting it in that canonical form keeps downstream
    irreversibility handling simple.

    Reactions with ``lb == 0`` (already forward) or with both bounds
    negative or both positive are not touched.

    Parameters
    ----------
    model
        A cobra.Model, mutated in place.

    Returns
    -------
    list of str
        Sorted IDs of reactions that were inverted.
    """
    inverted: list[str] = []
    for rxn in model.reactions:
        if rxn.lower_bound < 0 and rxn.upper_bound == 0:
            # Flip stoichiometry: every coefficient gets negated.
            flipped = {m: -coef for m, coef in rxn.metabolites.items()}
            # add_metabolites with combine=False overwrites the existing values.
            rxn.add_metabolites(flipped, combine=False)
            # Swap bounds. Assign upper_bound first to avoid the
            # lb > ub transient state that cobrapy rejects.
            new_ub = -rxn.lower_bound
            rxn.lower_bound = 0.0
            rxn.upper_bound = new_ub
            inverted.append(rxn.id)
    return sorted(inverted)