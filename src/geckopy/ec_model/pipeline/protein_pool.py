"""Add protein pseudometabolites and the protein pool machinery.

Covers stages 9, 10, 11, and 12 of GECKO MATLAB `makeEcModel.m`.

Reaction directions differ from the MATLAB convention. MATLAB defined
usage and pool reactions with `lb=-1000, ub=0` (negative-flux
reactions); geckopy defines them in the natural forward direction
(`lb=0, ub=1000`). A future I/O layer handles translation between the
two conventions so models remain compatible.
"""
from __future__ import annotations

from typing import TYPE_CHECKING

import cobra

if TYPE_CHECKING:
    from ..ec_model import EcModel


_PROT_PREFIX = "prot_"
_USAGE_PREFIX = "usage_prot_"
_POOL_ID = "prot_pool"
_POOL_EXCHANGE_ID = "prot_pool_exchange"
_PROTEIN_USAGE_SUBSYSTEM = "Protein usage"


# --------------------------------------------------------------------------- #
# Compartment resolution
# --------------------------------------------------------------------------- #

def _resolve_enzyme_compartment_id(model: "EcModel", requested: str) -> str:
    """Resolve ``adapter.params.enzyme_comp`` to a compartment ID.

    Matches ``requested`` first against compartment names (the typical
    case, e.g. "cytoplasm"), then falls back to compartment IDs (e.g.
    "c" when the model happens to use short IDs as names, as in
    ecTestGEM).

    Raises ValueError if neither matches.
    """
    compartments = model.compartments  # dict: id -> name

    # Prefer name match (matches MATLAB's strcmp(compNames, enzyme_comp)).
    for comp_id, comp_name in compartments.items():
        if comp_name == requested:
            return comp_id

    # Fallback: direct ID match.
    if requested in compartments:
        return requested

    raise ValueError(
        f"Enzyme compartment '{requested}' (adapter.params.enzyme_comp) "
        f"not found in model. Available compartments: "
        f"{dict(compartments)}"
    )


# --------------------------------------------------------------------------- #
# Stage 9: add one prot_<enzyme> pseudometabolite per enzyme
# --------------------------------------------------------------------------- #

def add_protein_pseudometabolites(model: "EcModel") -> list[str]:
    """Stage 9: add ``prot_<enzyme>`` pseudometabolites to the model.
    ... (docstring unchanged) ...
    """
    if model.adapter is None:
        raise ValueError("EcModel.adapter is None; cannot resolve enzyme_comp.")

    unique_enzymes = sorted(set(model.ec.enzymes))
    if not unique_enzymes:
        return []

    comp_id = _resolve_enzyme_compartment_id(
        model, model.adapter.params.enzyme_comp
    )

    new_mets: list[cobra.Metabolite] = []
    for enzyme in unique_enzymes:
        met_id = f"{_PROT_PREFIX}{enzyme}"
        if met_id in {m.id for m in model.metabolites}:
            continue
        met = cobra.Metabolite(
            id=met_id,
            name=met_id,
            compartment=comp_id,
        )
        met.annotation["sbo"] = "SBO:0000252"
        met.notes["enzyme_usage"] = "Enzyme-usage pseudometabolite"
        new_mets.append(met)

    if new_mets:
        model.add_metabolites(new_mets)

    return sorted(m.id for m in new_mets)


# --------------------------------------------------------------------------- #
# Stage 10: add prot_pool pseudometabolite
# --------------------------------------------------------------------------- #

def add_protein_pool_pseudometabolite(model: "EcModel") -> None:
    """Stage 10: add the single ``prot_pool`` pseudometabolite.

    Ported from GECKO MATLAB: src/geckomat/change_model/makeEcModel.m
    (stage 10).

    Parameters
    ----------
    model
        An EcModel with adapter set. Mutated in place.
    """
    if model.adapter is None:
        raise ValueError("EcModel.adapter is None; cannot resolve enzyme_comp.")

    if _POOL_ID in {m.id for m in model.metabolites}:
        return  # idempotent

    comp_id = _resolve_enzyme_compartment_id(
        model, model.adapter.params.enzyme_comp
    )

    pool = cobra.Metabolite(
        id=_POOL_ID,
        name=_POOL_ID,
        compartment=comp_id,
    )
    pool.notes["enzyme_usage"] = "Enzyme-usage protein pool"
    model.add_metabolites([pool])


# --------------------------------------------------------------------------- #
# Stage 11: add usage reactions
# --------------------------------------------------------------------------- #

def add_protein_usage_reactions(model: "EcModel") -> list[str]:
    """Stage 11: add ``usage_prot_<enzyme>`` reactions.

    Ported from GECKO MATLAB: src/geckomat/change_model/makeEcModel.m
    (stage 11, full-model branch only).

    Each usage reaction has stoichiometry ``prot_pool -> prot_<enzyme>``
    with bounds ``(0, 1000)``. This differs from the MATLAB convention,
    where the same reaction was written as ``prot_<enzyme> -> prot_pool``
    with bounds ``(-1000, 0)``, so that positive flux represents enzyme
    production in both cases. The two forms are equivalent; the geckopy
    form is the cobrapy-native forward direction.

    The GPR of each usage reaction is set to the single gene whose
    enzyme is the produced pseudometabolite, so that
    ``model.genes`` retains that gene even if no metabolic reaction
    references it anymore after upstream filtering.

    Parameters
    ----------
    model
        An EcModel with stage 9 already run. Mutated in place.

    Returns
    -------
    list of str
        Sorted IDs of added usage reactions.
    """
    pool_met = model.metabolites.get_by_id(_POOL_ID)

    new_rxns: list[cobra.Reaction] = []
    for enzyme, gene in sorted(
        set(zip(model.ec.enzymes, model.ec.genes))
    ):
        rxn_id = f"{_USAGE_PREFIX}{enzyme}"
        if rxn_id in {r.id for r in model.reactions}:
            continue  # idempotent

        met_id = f"{_PROT_PREFIX}{enzyme}"
        prot_met = model.metabolites.get_by_id(met_id)

        rxn = cobra.Reaction(id=rxn_id, name=rxn_id)
        rxn.lower_bound = 0.0
        rxn.upper_bound = 1000.0
        rxn.add_metabolites({pool_met: -1.0, prot_met: 1.0})
        rxn.gene_reaction_rule = gene
        rxn.subsystem = _PROTEIN_USAGE_SUBSYSTEM
        new_rxns.append(rxn)

    if new_rxns:
        model.add_reactions(new_rxns)

    return sorted(r.id for r in new_rxns)


# --------------------------------------------------------------------------- #
# Stage 12: add pool exchange reaction
# --------------------------------------------------------------------------- #

def add_protein_pool_exchange_reaction(model: "EcModel") -> None:
    """Stage 12: add the single ``prot_pool_exchange`` reaction.

    Ported from GECKO MATLAB: src/geckomat/change_model/makeEcModel.m
    (stage 12).

    Stoichiometry: ``(nothing) -> prot_pool`` with bounds ``(0, 1000)``.
    Positive flux imports protein into the pool, mirroring how cobrapy
    writes ordinary exchange reactions for extracellular metabolites.
    This is a direction flip vs the MATLAB convention (which uses
    ``prot_pool -> (nothing)`` with ``(-1000, 0)``); the two are
    equivalent.

    Parameters
    ----------
    model
        An EcModel with stage 10 already run. Mutated in place.
    """
    if _POOL_EXCHANGE_ID in {r.id for r in model.reactions}:
        return  # idempotent

    pool_met = model.metabolites.get_by_id(_POOL_ID)

    rxn = cobra.Reaction(id=_POOL_EXCHANGE_ID, name=_POOL_EXCHANGE_ID)
    rxn.lower_bound = 0.0
    rxn.upper_bound = 1000.0
    rxn.add_metabolites({pool_met: 1.0})
    rxn.subsystem = _PROTEIN_USAGE_SUBSYSTEM
    model.add_reactions([rxn])
