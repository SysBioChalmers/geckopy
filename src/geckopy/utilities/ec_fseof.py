"""Flux-Scanning with Enforced Objective Function for an ecModel.

Thin wrapper around :func:`raven_toolbox.analysis.fseof.fseof` that adds
ec-specific concerns on top of raven's general-purpose FSEOF:

- resolves the biomass reaction from ``model.adapter`` when not passed
  explicitly;
- optionally runs a carbon-source consistency check (warns when the
  carbon-source lower bound is tighter than its uptake at biomass-max);
- filters ``usage_prot_*`` rows out of the returned ``scan`` /
  ``targets`` DataFrames — those reactions are protein-pool plumbing,
  never engineering targets.

The selection method and output shape come from raven-toolbox: each
reaction's flux is regressed against the enforced product flux across
the scan, and a reaction is a target if ``|correlation| >=
correlation_threshold`` and ``|slope| >= flux_eps``. Targets are
classified as ``amplify`` / ``knockdown`` / ``knockout``. See
:class:`raven_toolbox.analysis.fseof.FSEOFResult` for the full output
schema.

Ported from GECKO MATLAB: ``src/geckomat/utilities/ecFSEOF.m``,
re-designed around :func:`raven_toolbox.analysis.fseof.fseof`.

MATLAB-COMPAT: GECKO MATLAB uses strict monotonicity + top-25%-by-slope
to pick targets; geckopy delegates to raven-toolbox's regression-based
selection (``|correlation|`` against the enforced flux) with pFBA per
step. Regression is more robust to LP alternative optima — see
raven-toolbox's IMPROVEMENTS notes FS1-FS4 for the rationale.

MATLAB-COMPAT: GECKO MATLAB classifies targets as ``OE`` / ``KD`` /
``KO``; raven uses ``amplify`` / ``knockdown`` / ``knockout``. Same
semantics, different vocabulary.

MATLAB-COMPAT: GECKO MATLAB optionally writes the result tables to
TSV. geckopy returns the raven ``FSEOFResult``; callers can
``result.targets.to_csv(...)`` themselves.

MATLAB-COMPAT: GECKO MATLAB computes a per-gene ``essentiality`` column
by blocking each gene's ``usage_prot_<enzyme>`` reactions and re-solving.
geckopy drops this column. The same analysis is available via
:func:`cobra.flux_analysis.single_gene_deletion`, which works on
ecModels because GPRs already gate the catalysed reactions.

MATLAB-COMPAT: GECKO MATLAB splits targets into ``rxn_targets`` vs
``transport_targets`` (reactions whose metabolites span compartments).
geckopy returns one combined ``targets`` table; callers can split it
themselves via the ``subsystem`` column or by inspecting each reaction.
"""
from __future__ import annotations

import logging
from typing import TYPE_CHECKING, Optional

from raven_toolbox.analysis.fseof import FSEOFResult, fseof

if TYPE_CHECKING:
    from ..ec_model.ec_model import EcModel


logger = logging.getLogger(__name__)

_USAGE_PROT_PREFIX = "usage_prot_"


def ec_fseof(
    model: "EcModel",
    prod_target_rxn: str,
    cs_rxn: Optional[str] = None,
    *,
    n_steps: int = 10,
    bio_rxn: Optional[str] = None,
    max_fraction: float = 0.9,
    correlation_threshold: float = 0.9,
    flux_eps: float = 1e-6,
) -> FSEOFResult:
    """Run Flux-Scanning with Enforced Objective Function on an ecModel.

    Enforces an increasing flux toward ``prod_target_rxn`` while
    optimising biomass, runs pFBA at each step, regresses every
    reaction's flux against the enforced product flux, and returns
    reactions with strong slope-vs-correlation signal as engineering
    targets.

    Parameters
    ----------
    model
        EcModel with the protein-pool / usage-rxn machinery.
    prod_target_rxn
        Reaction id for the production target (typically a
        product-exchange).
    cs_rxn
        Reaction id for the main carbon source. When passed, a warning
        is emitted if its lower bound is tighter than its uptake at
        biomass-max (a common mis-configuration). ``None`` skips the
        check.
    n_steps
        Number of enforced-flux levels in the scan.
    bio_rxn
        Biomass reaction id. Defaults to ``adapter.params.bio_rxn`` when
        ``model.adapter`` is set.
    max_fraction
        Top of the scan range as a fraction of the theoretical maximum
        product flux. Default 0.9.
    correlation_threshold
        A reaction is a target when ``|corr(flux, enforced)| >=`` this
        value. Default 0.9.
    flux_eps
        Numerical floor for flat-flux detection.

    Returns
    -------
    FSEOFResult
        raven-toolbox's FSEOF output dataclass, with ``usage_prot_*``
        rows filtered out of both ``scan`` and ``targets``. See
        :class:`raven_toolbox.analysis.fseof.FSEOFResult` for the full
        field list (``scan``, ``enforced``, ``targets``, plus the
        ``amplification`` / ``knockout`` / ``gene_targets`` properties).
    """
    from ..adapter import resolve_param
    bio_rxn_id = resolve_param(
        model, bio_rxn, "bio_rxn",
        purpose="ec_fseof needs the biomass reaction id",
    )

    if cs_rxn is not None:
        _check_cs_consistency(model, bio_rxn_id, cs_rxn)

    result = fseof(
        model, prod_target_rxn,
        biomass_rxn=bio_rxn_id,
        n_steps=n_steps,
        max_fraction=max_fraction,
        correlation_threshold=correlation_threshold,
        flux_eps=flux_eps,
    )
    return _drop_usage_prot(result)


def _check_cs_consistency(
    model: "EcModel", bio_rxn_id: str, cs_rxn_id: str,
) -> None:
    """Warn when the carbon source can't deliver its biomass-max uptake.

    Matches the MATLAB ecFSEOF safety check: if the carbon-source
    reaction's lower bound is less negative than its biomass-max flux,
    the scan will silently flatten at the bound and the slopes mean
    less than they look.
    """
    cs_rxn_obj = model.reactions.get_by_id(cs_rxn_id)
    with model:
        model.objective = bio_rxn_id
        sol = model.optimize()
        if sol.status != "optimal":
            return
        cs_flux = float(sol.fluxes.get(cs_rxn_id, 0.0))
    if cs_rxn_obj.lower_bound < cs_flux:
        logger.warning(
            "ec_fseof: carbon source %r lower bound is %g but uptake "
            "at biomass-max is %g; consider tightening so the scan "
            "doesn't flatten at the bound.",
            cs_rxn_id, cs_rxn_obj.lower_bound, cs_flux,
        )


def _drop_usage_prot(result: FSEOFResult) -> FSEOFResult:
    """Strip ``usage_prot_*`` rows from raven's ``scan`` and ``targets``.

    These reactions exist solely to draw the per-enzyme prot mets out of
    the shared pool; they have no metabolic meaning of their own and
    can't be engineered. Filtering at the boundary keeps the rest of
    raven's output structure intact (the ``gene_targets`` property and
    the ``amplification`` / ``knockout`` slices recompute correctly off
    the filtered ``targets``).
    """
    scan = result.scan[
        ~result.scan.index.astype(str).str.startswith(_USAGE_PROT_PREFIX)
    ]
    targets = result.targets[
        ~result.targets["reaction"].str.startswith(_USAGE_PROT_PREFIX)
    ].reset_index(drop=True)
    return FSEOFResult(scan=scan, enforced=result.enforced, targets=targets)
