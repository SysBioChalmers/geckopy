"""Top-level orchestrator: build an EcModel from a conventional GEM.

Ported from GECKO MATLAB: src/geckomat/change_model/makeEcModel.m.

This function runs all 12 stages of the original makeEcModel in order,
returning the populated EcModel plus the list of model genes not found
in the UniProt database.
"""
from __future__ import annotations

from typing import TYPE_CHECKING, Optional

import cobra

from ..databases import UniprotDB, load_uniprot_tsv
from .ec_model import EcModel
from .pipeline import (
    add_protein_pool_exchange_reaction,
    add_protein_pool_pseudometabolite,
    add_protein_pseudometabolites,
    add_protein_usage_reactions,
    allocate_ec_for_catalyzed_reactions,
    build_rxn_enzyme_coupling,
    convert_to_irreversible,
    expand_model,
    invert_backwards_only_reactions,
    populate_enzyme_data,
    remove_pseudoreaction_gprs,
)

if TYPE_CHECKING:
    from ..adapter import ModelAdapter


def make_ec_model(
    model: cobra.Model,
    adapter: "ModelAdapter",
    *,
    gecko_light: bool = False,
    uniprot_db: Optional[UniprotDB] = None,
) -> tuple[EcModel, list[str]]:
    """Expand a conventional GEM into a basic enzyme-constrained model.

    Ported from GECKO MATLAB: src/geckomat/change_model/makeEcModel.m.

    Runs the full 12-stage pipeline:

    1. Remove GPRs from pseudoreactions.
    2. Invert reactions constrained to negative flux only.
    3. (Merged into stage 4.)
    4. Convert non-exchange reversible reactions to forward/reverse pair.
    5. Expand reactions with isozymes (OR in GPR) into one per isozyme.
    6. Allocate empty ec substructure (per-reaction slots).
    7. Populate ec with per-enzyme UniProt data.
    8. Build the reaction-enzyme coupling matrix.
    9. Add ``prot_<enzyme>`` pseudometabolites.
    10. Add the ``prot_pool`` pseudometabolite.
    11. Add ``usage_prot_<enzyme>`` reactions.
    12. Add the ``prot_pool_exchange`` reaction.

    kcat values are not set here. Call applyKcatConstraints (future)
    after populating ec.kcat with values from gather_kcats functions.

    Parameters
    ----------
    model
        A conventional cobra.Model. Mutated in place by stages 1-5 and
        then wrapped as an EcModel for stages 6-12.
    adapter
        Loaded ModelAdapter providing organism parameters and the
        location of the UniProt TSV.
    gecko_light
        If True, generate a simplified gecko-light model. Not yet
        implemented; raises NotImplementedError.
    uniprot_db
        Pre-loaded UniprotDB. If None, loaded from
        ``adapter.params.path / "data" / "uniprot.tsv"``. Provided as a
        parameter mainly for testing.

    Returns
    -------
    ec_model
        An EcModel with populated ec substructure, protein
        pseudometabolites, and pool machinery.
    no_uniprot
        List of model gene IDs that could not be matched to any UniProt
        entry. These reactions are flagged via rxn.notes but remain in
        the model.

    Raises
    ------
    NotImplementedError
        If ``gecko_light`` is True.
    FileNotFoundError
        If ``uniprot_db`` is None and no uniprot.tsv is found.
    ValueError
        Propagated from the individual stages on any consistency error
        (unresolvable compartment, no UniProt matches, etc.).
    """
    if gecko_light:
        raise NotImplementedError(
            "gecko_light mode is not yet implemented in geckopy. "
            "Only full ecModels are currently supported."
        )

    # Pre-flight: refuse to run twice on the same model. An existing
    # model.ec populated substructure is a strong signal the model has
    # already been converted.
    if isinstance(model, EcModel) and model.ec.n_rxns > 0:
        raise ValueError(
            "make_ec_model was called on a model that already has a "
            "populated ec substructure. Run it only on a conventional GEM."
        )

    # Stages 1-5: preprocess on a plain cobra.Model.
    remove_pseudoreaction_gprs(model, adapter)
    invert_backwards_only_reactions(model)
    convert_to_irreversible(model)
    expand_model(model)

    # Promote to EcModel for stages 6-12.
    ec_model = EcModel.from_cobra(model, adapter, gecko_light=gecko_light)

    # Stages 6-8: build the ec substructure.
    allocate_ec_for_catalyzed_reactions(ec_model)

    if uniprot_db is None:
        uniprot_path = adapter.params.path / "data" / "uniprot.tsv"
        uniprot_db = load_uniprot_tsv(uniprot_path)

    no_uniprot = populate_enzyme_data(ec_model, uniprot_db)
    build_rxn_enzyme_coupling(ec_model)

    # Stages 9-12: protein pool machinery.
    add_protein_pseudometabolites(ec_model)
    add_protein_pool_pseudometabolite(ec_model)
    add_protein_usage_reactions(ec_model)
    add_protein_pool_exchange_reaction(ec_model)

    return ec_model, no_uniprot
