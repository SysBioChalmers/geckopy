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
from raven_toolbox.utils.sort import sort_identifiers

from ..databases import UniprotDB, load_uniprot_tsv
from .ec_model import EcModel
from .pipeline import (
    add_protein_pool_exchange_reaction,
    add_protein_pool_pseudometabolite,
    add_protein_pseudometabolites,
    add_protein_usage_reactions,
    allocate_ec_and_coupling_light,
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
    from ..databases.kegg_loader import KeggDB


logger = logging.getLogger(__name__)


def make_ec_model(
    model: cobra.Model,
    adapter: "ModelAdapter",
    *,
    gecko_light: bool = False,
    uniprot_db: Optional[UniprotDB] = None,
    kegg_db: "Optional[KeggDB]" = None,
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
    at 0 (0 marks "no kcat assigned") until you populate them
    via the ``gather_kcats`` functions (BRENDA / DLKcat / custom
    files) and call ``apply_kcat_constraints``, or set them manually
    via ``set_kcat_for_reactions``.

    Parameters
    ----------
    model
        A conventional cobra.Model (the starting GEM). Left
        unchanged: the function operates on an internal copy, so the
        caller's model is never mutated (matching the MATLAB GECKO
        ``makeEcModel`` value semantics).
    adapter
        A loaded ModelAdapter. Carries organism parameters (taxonomy
        id, biomass reaction, sigma factor, ...) and the path to
        the UniProt cache.
    gecko_light
        Build a light ecModel instead of a full one. The light layout
        skips the isozyme split (cobra reactions stay singular), omits
        the per-enzyme ``prot_<id>`` pseudometabolites and
        ``usage_prot_<id>`` reactions, and keeps only the shared
        ``prot_pool``. The per-isozyme bookkeeping moves into ``ec``
        instead: each cobra reaction with N isozymes produces N rows
        in ``ec.rxns`` distinguished by a 3-digit counter prefix
        (e.g. ``001_RXN1``, ``002_RXN1``). Suitable for large GEMs
        (HumanGEM, Recon3D) where the full layout's LP size becomes
        impractical. See ``docs/gecko_light_status.md``.
    uniprot_db
        Pre-loaded UniprotDB. If None, the function looks for
        ``adapter.params.path / "data" / "uniprot.tsv"`` and loads
        it automatically.
    kegg_db
        Pre-loaded KeggDB. Optional fallback source for protein
        information and EC numbers. When provided, stage 7 fills
        genes that UniProt missed (with ``ec.enzymes`` set to the
        UniProt accession carried in the KEGG row, or to the KEGG
        gene id when that accession is empty), and the EC-number
        lookup (``fill_eccodes_from_database``) consults KEGG when
        the UniProt result is empty or ends with ``-``. If None,
        KEGG is not used.

    Returns
    -------
    EcModel
        The built ecModel. ``model.ec.kcat`` is all 0; fill it in
        before solving.

    Raises
    ------
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
    if isinstance(model, EcModel) and model.ec.n_rxns > 0:
        raise ValueError(
            "make_ec_model was called on a model that already has a "
            "populated ec substructure. Run it only on a conventional GEM."
        )

    # Work on a copy: the pipeline stages below mutate the model in place
    # (splitting reactions, adding pseudometabolites, ...), but the caller's
    # input GEM must be left untouched.
    model = model.copy()

    # Stages 1-4: preprocess on a plain cobra.Model. Stage 5 (isozyme
    # expansion) is full-only; light keeps isozymes singular and tracks
    # them per-row in ec instead.
    remove_pseudoreaction_gprs(model, adapter)
    invert_backwards_only_reactions(model)
    convert_to_irreversible(model)
    if not gecko_light:
        expand_model(model)
        # Sort reactions, so that reversible and isozymic reactions are kept
        # near each other. MATLAB makeEcModel does this here (full models
        # only), before the protein pseudoreactions are appended, so the
        # metabolic block comes out alphabetically and the protein block
        # follows in enzyme order.
        sort_identifiers(model)

    # Promote to EcModel for stages 6-12.
    ec_model = EcModel.from_cobra(model, adapter, gecko_light=gecko_light)

    if uniprot_db is None:
        uniprot_path = adapter.params.path / "data" / "uniprot.tsv"
        uniprot_db = load_uniprot_tsv(uniprot_path)

    if gecko_light:
        # Light: stage 7 (per-enzyme data) runs first so the coupling
        # builder can index ec.genes when emitting the per-isozyme rows.
        # Stage 6 is folded into the light helper, which writes both
        # ec.rxns (one row per isozyme, with the 3-digit counter prefix)
        # and ec.rxn_enz_mat in one pass.
        no_uniprot = populate_enzyme_data(
            ec_model, uniprot_db, kegg_db=kegg_db,
        )
        allocate_ec_and_coupling_light(ec_model)
    else:
        # Full: stage 6 (allocate empty slots) runs first; stage 7
        # populates per-enzyme arrays; stage 8 fills the coupling matrix
        # from cobra GPRs (one AND-clause per reaction after the stage-5
        # expansion).
        allocate_ec_for_catalyzed_reactions(ec_model)
        no_uniprot = populate_enzyme_data(
            ec_model, uniprot_db, kegg_db=kegg_db,
        )
        build_rxn_enzyme_coupling(ec_model)

    # Stages 9-12: protein pool machinery. Light skips the per-enzyme
    # pseudometabolites + usage reactions (no per-enzyme bookkeeping in
    # the LP) but keeps the shared pool exactly as full models do.
    if not gecko_light:
        add_protein_pseudometabolites(ec_model)
    add_protein_pool_pseudometabolite(ec_model)
    if not gecko_light:
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

    # Catch any internal length/shape drift in the ec arrays before the
    # model is handed back (cheap, and surfaces build bugs early).
    ec_model.ec.validate()

    return ec_model
