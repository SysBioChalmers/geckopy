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


from ..constants import (
    POOL_EXCHANGE_ID,
    POOL_ID,
    PROT_PREFIX,
    PROTEIN_USAGE_SUBSYSTEM,
    USAGE_PREFIX,
)


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
    from ...adapter import resolve_adapter
    adapter = resolve_adapter(
        model,
        purpose="add_protein_pseudometabolites reads "
        "params.enzyme_comp from the adapter",
    )

    unique_enzymes = sorted(set(model.ec.enzymes))
    if not unique_enzymes:
        return []

    comp_id = _resolve_enzyme_compartment_id(
        model, adapter.params.enzyme_comp
    )

    new_mets: list[cobra.Metabolite] = []
    for enzyme in unique_enzymes:
        met_id = f"{PROT_PREFIX}{enzyme}"
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
    from ...adapter import resolve_adapter
    adapter = resolve_adapter(
        model,
        purpose="add_protein_pool_pseudometabolite reads "
        "params.enzyme_comp from the adapter",
    )

    if POOL_ID in {m.id for m in model.metabolites}:
        return  # idempotent

    comp_id = _resolve_enzyme_compartment_id(
        model, adapter.params.enzyme_comp
    )

    pool = cobra.Metabolite(
        id=POOL_ID,
        name=POOL_ID,
        compartment=comp_id,
    )
    pool.notes["enzyme_usage"] = "Enzyme-usage protein pool"
    model.add_metabolites([pool])


# --------------------------------------------------------------------------- #
# Stage 11: add usage reactions
# --------------------------------------------------------------------------- #

def add_protein_usage_reactions(model: "EcModel") -> list[str]:
    """Stage 11: add ``usage_prot_<enzyme>`` reactions.

    MATLAB-COMPAT: GECKO MATLAB writes these reactions as
    ``prot_<enzyme> -> prot_pool`` with bounds ``(-1000, 0)``, so flux
    is negative when the enzyme is being produced. geckopy uses the
    forward direction ``prot_pool -> prot_<enzyme>`` with bounds
    ``(0, 1000)``, matching cobrapy convention. To make the two
    representations interchangeable, MATLAB GECKO should switch to
    the forward direction in a future release. The model I/O layer
    (when written) detects and translates the old form on read.


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
    pool_met = model.metabolites.get_by_id(POOL_ID)

    new_rxns: list[cobra.Reaction] = []
    for enzyme, gene in sorted(
        set(zip(model.ec.enzymes, model.ec.genes))
    ):
        rxn_id = f"{USAGE_PREFIX}{enzyme}"
        if rxn_id in {r.id for r in model.reactions}:
            continue  # idempotent

        met_id = f"{PROT_PREFIX}{enzyme}"
        prot_met = model.metabolites.get_by_id(met_id)

        rxn = cobra.Reaction(id=rxn_id, name=rxn_id)
        rxn.lower_bound = 0.0
        rxn.upper_bound = 1000.0
        rxn.add_metabolites({pool_met: -1.0, prot_met: 1.0})
        rxn.gene_reaction_rule = gene
        rxn.subsystem = PROTEIN_USAGE_SUBSYSTEM
        new_rxns.append(rxn)

    if new_rxns:
        model.add_reactions(new_rxns)

    return sorted(r.id for r in new_rxns)


# --------------------------------------------------------------------------- #
# Stage 12: add pool exchange reaction
# --------------------------------------------------------------------------- #

def add_protein_pool_exchange_reaction(model: "EcModel") -> None:
    """Stage 12: add the single ``prot_pool_exchange`` reaction.

    MATLAB-COMPAT: GECKO MATLAB writes this as
    ``prot_pool -> (nothing)`` with bounds ``(-1000, 0)``, so negative
    flux imports protein. geckopy uses the forward direction
    ``(nothing) -> prot_pool`` with bounds ``(0, 1000)``. As with the
    usage reactions, MATLAB GECKO should switch to the forward
    direction.


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
    if POOL_EXCHANGE_ID in {r.id for r in model.reactions}:
        return  # idempotent

    pool_met = model.metabolites.get_by_id(POOL_ID)

    rxn = cobra.Reaction(id=POOL_EXCHANGE_ID, name=POOL_EXCHANGE_ID)
    rxn.lower_bound = 0.0
    rxn.upper_bound = 1000.0
    rxn.add_metabolites({pool_met: 1.0})
    rxn.subsystem = PROTEIN_USAGE_SUBSYSTEM
    model.add_reactions([rxn])


# --------------------------------------------------------------------------- #
# setProtPoolSize (post-makeEcModel utility, same subsystem)
# --------------------------------------------------------------------------- #

def set_prot_pool_size(
    model: "EcModel",
    *,
    p_tot: float | None = None,
    f: float | None = None,
    sigma: float | None = None,
) -> float:
    """Set the upper bound of ``prot_pool_exchange`` (the protein budget).

    The protein pool exchange represents the cell's total enzyme
    supply. Its upper bound is what limits how much protein the
    model can deploy at once. The formula is:

        bound = p_tot * f * sigma * 1000

    where:

    - ``p_tot`` is the total cellular protein content (g/gDCW),
    - ``f`` is the fraction of that protein that's made up of
      enzymes in the model (g model-enzyme / g total-protein),
    - ``sigma`` is the average enzyme saturation factor (0-1) —
      how close enzymes run to their Vmax in vivo,
    - the ``1000`` converts units (``ec.mw`` is in Da and fluxes
      are in mmol/gDW/h, so multiplying ``p_tot`` (g/gDW) by 1000
      brings it to mg/gDW, matching the mg/mmol that Da
      represents).

    When called with no arguments, the function reads ``p_tot``,
    ``f``, and ``sigma`` from the adapter. Any of the three can
    be overridden via keyword.

    MATLAB-COMPAT: GECKO MATLAB sets the (negative) lower bound;
    geckopy sets the (positive) upper bound. Numerically the same
    in absolute terms. When MATLAB GECKO switches its pool
    exchange to the forward direction, that MATLAB function will
    also switch from setting lb to setting ub.

    Ported from GECKO MATLAB: src/geckomat/change_model/setProtPoolSize.m.

    Parameters
    ----------
    model
        An EcModel with ``prot_pool_exchange`` already added (stage 12).
    p_tot, f, sigma
        Override values. If None, read from ``model.adapter.params``.

    Returns
    -------
    float
        The bound that was set (for logging or assertions).

    Raises
    ------
    ValueError
        If the model has no adapter and any argument is None, or if
        ``prot_pool_exchange`` is not in the model.
    """
    if p_tot is None or f is None or sigma is None:
        from ...adapter import resolve_adapter
        adapter = resolve_adapter(
            model,
            purpose="set_prot_pool_size reads p_tot / f / sigma "
            "defaults from the adapter (or pass them explicitly)",
        )
        params = adapter.params
        if p_tot is None:
            p_tot = params.p_tot
        if f is None:
            f = params.f
        if sigma is None:
            sigma = params.sigma

    if POOL_EXCHANGE_ID not in {r.id for r in model.reactions}:
        raise ValueError(
            f"Reaction '{POOL_EXCHANGE_ID}' not found. Call make_ec_model "
            f"or add_protein_pool_exchange_reaction first."
        )

    bound = p_tot * f * sigma * 1000.0
    model.reactions.get_by_id(POOL_EXCHANGE_ID).upper_bound = bound
    return bound
