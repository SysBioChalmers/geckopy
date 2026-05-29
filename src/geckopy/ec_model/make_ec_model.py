"""Top-level orchestrator: build an EcModel from a conventional GEM.

Ported from GECKO MATLAB: src/geckomat/change_model/makeEcModel.m.
"""
from __future__ import annotations

import logging
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


logger = logging.getLogger(__name__)


def make_ec_model(
    model: cobra.Model,
    adapter: "ModelAdapter",
    *,
    gecko_light: bool = False,
    uniprot_db: Optional[UniprotDB] = None,
) -> EcModel:
    """Expand a conventional GEM into a basic enzyme-constrained model.

    Ported from GECKO MATLAB: src/geckomat/change_model/makeEcModel.m.

    Runs the full 12-stage pipeline and returns the populated EcModel.
    Genes that could not be matched to a UniProt entry are reported as
    a logged warning and individually annotated on each affected
    reaction via ``rxn.notes["geckopy_warning"]`` (set by stage 7).

    MATLAB-COMPAT: GECKO MATLAB returns the list of unmatched genes as
    a separate output (``noUniprot``). geckopy logs a warning summary
    instead, since the per-reaction annotations are typically more
    useful for debugging than the flat list.

    kcat values are not set here. Call apply_kcat_constraints (after
    populating ec.kcat with values from gather_kcats functions, or via
    set_kcat_for_reactions).

    Parameters
    ----------
    model
        A conventional cobra.Model. Mutated in place by stages 1-5 and
        then wrapped as an EcModel for stages 6-12.
    adapter
        Loaded ModelAdapter providing organism parameters and the
        location of the UniProt TSV.
    gecko_light
        Not yet implemented; raises NotImplementedError.
    uniprot_db
        Pre-loaded UniprotDB. If None, loaded from
        ``adapter.params.path / "data" / "uniprot.tsv"``.

    Returns
    -------
    EcModel
        An EcModel with populated ec substructure, protein
        pseudometabolites, and pool machinery.

    Raises
    ------
    NotImplementedError
        If ``gecko_light`` is True.
    FileNotFoundError
        If ``uniprot_db`` is None and no uniprot.tsv is found.
    ValueError
        Propagated from individual stages on consistency errors.
    """
    if gecko_light:
        raise NotImplementedError(
            "gecko_light mode is not yet implemented in geckopy."
        )

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

    if no_uniprot:
        preview = ", ".join(no_uniprot[:5])
        more = "" if len(no_uniprot) <= 5 else f" (and {len(no_uniprot) - 5} more)"
        logger.warning(
            "%d gene(s) not found in UniProt and left enzyme-unconstrained: "
            "%s%s. Affected reactions are annotated via "
            "rxn.notes['geckopy_warning'].",
            len(no_uniprot),
            preview,
            more,
        )

    return ec_model
