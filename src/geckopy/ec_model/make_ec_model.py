"""Build an EcModel from a conventional GEM.

``make_ec_model`` is the top-level entry point for turning a plain
cobra metabolic model into an enzyme-constrained model. It runs
the 12-stage GECKO pipeline (preprocess, expand isozymes, look up
UniProt data, add protein pseudometabolites, set up the shared
pool, ...) and returns a populated ``EcModel``.

Note that kcat values are NOT set by this function. You usually
follow up with one of:

- ``geckopy.gather_kcats`` functions (fetch from BRENDA, DLKcat,
  custom files), then ``apply_kcat_constraints`` to push the
  values into the LP, or
- ``ec_model.pipeline.set_kcat.set_kcat_for_reactions`` for a
  small manual set.

The MATLAB tutorials in
``tutorials/full_ecModel/protocol.py`` show the full flow.

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
    """Build an enzyme-constrained model from a conventional GEM.

    Runs the 12-stage GECKO pipeline end to end and returns the
    populated ``EcModel``. The output has:

    - Reactions split per isozyme where applicable (``_EXP_<N>``
      suffix) and per direction for reversible ones (``_REV``).
    - One ``prot_<uniprot>`` pseudo-metabolite per enzyme.
    - One ``usage_prot_<uniprot>`` reaction per enzyme drawing
      from the shared protein pool.
    - A ``prot_pool`` pseudo-metabolite and a
      ``prot_pool_exchange`` reaction acting as the global
      protein budget.
    - A populated ``model.ec`` substructure (kcats, MW, sequences,
      coupling matrix).

    The kcat values themselves are NOT filled in here — they stay
    NaN until you populate them via the ``gather_kcats`` functions
    (BRENDA / DLKcat / custom files) and call
    ``apply_kcat_constraints``, or set them manually via
    ``set_kcat_for_reactions``.

    Parameters
    ----------
    model
        A conventional cobra.Model (the starting GEM). Mutated in
        place by stages 1-5 (preprocessing + isozyme expansion)
        and then wrapped as an EcModel for stages 6-12.
    adapter
        A loaded ModelAdapter. Carries organism parameters (taxonomy
        id, biomass reaction, sigma factor, ...) and the path to
        the UniProt cache.
    gecko_light
        Build a light ecModel instead of a full one. Not yet
        implemented in geckopy; raises NotImplementedError. See
        ``docs/gecko_light_status.md``.
    uniprot_db
        Pre-loaded UniprotDB. If None, the function looks for
        ``adapter.params.path / "data" / "uniprot.tsv"`` and loads
        it automatically.

    Returns
    -------
    EcModel
        The built ecModel. ``model.ec.kcat`` is all NaN; fill it in
        before solving.

    Raises
    ------
    NotImplementedError
        If ``gecko_light`` is True (not yet supported).
    FileNotFoundError
        If ``uniprot_db`` is None and no ``uniprot.tsv`` is found
        under the adapter's data folder.
    ValueError
        Propagated from individual pipeline stages on
        inconsistencies (e.g. genes referenced in GPRs that don't
        exist in the model).

    MATLAB-COMPAT: GECKO MATLAB returns the list of unmatched
    genes as a second output (``noUniprot``). geckopy logs a
    warning summary instead and annotates each affected reaction
    via ``rxn.notes["geckopy_warning"]`` — the annotations are
    usually more useful for debugging than the flat list.

    Ported from GECKO MATLAB: src/geckomat/change_model/makeEcModel.m.
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
