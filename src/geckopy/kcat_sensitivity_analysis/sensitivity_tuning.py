"""Iteratively raise the most-limiting kcats until a target growth rate
is reached.

Ported from GECKO MATLAB:
src/geckomat/kcat_sensitivity_analysis/sensitivityTuning.m.
"""
from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Optional

import numpy as np

from ..ec_model.pipeline.apply_kcat import apply_kcat_constraints

if TYPE_CHECKING:
    from ..ec_model.ec_model import EcModel


logger = logging.getLogger(__name__)


_USAGE_PREFIX = "usage_prot_"
_PROT_PREFIX = "prot_"
_TUNING_SOURCE = "sensitivityTuning"


@dataclass
class TunedKcatsResult:
    """Information about which kcats were tuned and by how much."""

    rxns: list[str] = field(default_factory=list)
    rxn_names: list[str] = field(default_factory=list)
    enzymes: list[str] = field(default_factory=list)
    old_kcat: np.ndarray = field(
        default_factory=lambda: np.empty(0, dtype=float)
    )
    new_kcat: np.ndarray = field(
        default_factory=lambda: np.empty(0, dtype=float)
    )
    source: list[str] = field(default_factory=list)


def sensitivity_tuning(
    model: "EcModel",
    *,
    desired_growth_rate: Optional[float] = None,
    fold_change: float = 10.0,
    prot_to_ignore: Optional[list[str]] = None,
    verbose: bool = True,
) -> TunedKcatsResult:
    """Iteratively bump the most-limiting kcat by ``fold_change`` until
    the model can reach ``desired_growth_rate``.

    Each iteration:

    1. Solves the LP at the current kcat values.
    2. If growth didn't increase since the last iteration, warns and
       breaks (likely a substrate / pool constraint, not a kcat one).
    3. If growth reached the target, returns.
    4. Otherwise picks the ``usage_prot_*`` reaction with the highest
       flux (= most-produced enzyme), excluding any whose UniProt ID
       is in ``prot_to_ignore``.
    5. Finds the catalysed reaction consuming the most of that
       enzyme's pseudometabolite and multiplies that reaction's
       ``ec.kcat`` by ``fold_change``.
    6. Annotates ``ec.notes`` with the original kcat and source
       (unless the source is already ``"sensitivityTuning"``).
    7. Re-applies kcat constraints for that single reaction.

    Ported from GECKO MATLAB:
    src/geckomat/kcat_sensitivity_analysis/sensitivityTuning.m.

    MATLAB-COMPAT: MATLAB usage rxns go reverse (most production =
    most negative flux); geckopy uses forward (most production =
    most positive flux). MATLAB picks ``min(drawFluxes)``, geckopy
    picks ``max(drawFluxes)``. Same enzyme either way.

    MATLAB-COMPAT: The MATLAB function has a separate gecko-light
    branch. geckopy does not yet support gecko-light pipelines;
    calling this on a gecko-light model raises
    ``NotImplementedError`` (consistent with
    ``constrain_enz_concs``).

    MATLAB-COMPAT: ``modelAdapter`` is dropped per established
    convention; ``desired_growth_rate`` defaults to
    ``model.adapter.params.gr_exp``.

    Parameters
    ----------
    model
        EcModel with the protein pool / usage rxn machinery
        installed and a populated ``ec.kcat``. Mutated in place.
    desired_growth_rate
        Target growth rate. Defaults to
        ``model.adapter.params.gr_exp``.
    fold_change
        Per-iteration kcat multiplier.
    prot_to_ignore
        UniProt IDs whose usage rxns are excluded from being picked
        as the limiting enzyme.
    verbose
        Whether per-iteration progress is logged at INFO.

    Returns
    -------
    TunedKcatsResult

    Raises
    ------
    ValueError
        If the initial LP is infeasible, or if
        ``desired_growth_rate`` is needed but no adapter / value is
        available.
    NotImplementedError
        If ``model.ec.gecko_light`` is True.
    """
    if model.ec.gecko_light:
        raise NotImplementedError(
            "sensitivity_tuning does not yet support gecko-light models."
        )

    if desired_growth_rate is None:
        if model.adapter is None or model.adapter.params.gr_exp is None:
            raise ValueError(
                "desired_growth_rate not provided and "
                "model.adapter.params.gr_exp is unavailable."
            )
        desired_growth_rate = float(model.adapter.params.gr_exp)

    if prot_to_ignore is None:
        prot_to_ignore = []

    bio_rxn_id = (
        model.adapter.params.bio_rxn if model.adapter is not None else None
    )
    if bio_rxn_id:
        model.objective = bio_rxn_id

    initial_sol = model.optimize()
    if (
        initial_sol.objective_value is None
        or np.isnan(initial_sol.objective_value)
    ):
        raise ValueError(
            "FBA of the input model is infeasible. Reduce the protein "
            "pool constraint with set_prot_pool_size and/or check the "
            "exchange constraints."
        )

    usage_rxn_ids = [
        r.id for r in model.reactions
        if r.id.startswith(_USAGE_PREFIX)
    ]
    ignore_usage_ids = {f"{_USAGE_PREFIX}{p}" for p in prot_to_ignore}

    last_growth = float("nan")
    tuned_ec_indices: list[int] = []  # indices into model.ec.rxns
    iteration = 1

    while True:
        sol = model.optimize()
        growth = float(sol.objective_value or 0.0)

        if growth == last_growth:
            logger.warning(
                "sensitivity_tuning: growth did not increase from %g; "
                "uptake or pool constraints may be limiting.", growth,
            )
            break
        last_growth = growth
        if verbose:
            logger.info(
                "sensitivity_tuning: iteration %d, growth = %g",
                iteration, growth,
            )
        if growth >= desired_growth_rate:
            break

        iteration += 1

        # Pick usage rxn with max flux (= most produced enzyme),
        # excluding any in prot_to_ignore.
        best_usage_id: Optional[str] = None
        best_usage_flux = -np.inf
        for uid in usage_rxn_ids:
            if uid in ignore_usage_ids:
                continue
            f = float(sol.fluxes.get(uid, 0.0))
            if f > best_usage_flux:
                best_usage_flux = f
                best_usage_id = uid

        if best_usage_id is None or best_usage_flux <= 0:
            logger.warning(
                "sensitivity_tuning: no usage rxn carries flux; cannot "
                "continue tuning."
            )
            break

        # The enzyme metabolite the usage rxn produces.
        enzyme_id = best_usage_id[len(_USAGE_PREFIX):]
        prot_met_id = f"{_PROT_PREFIX}{enzyme_id}"
        try:
            prot_met = model.metabolites.get_by_id(prot_met_id)
        except KeyError:
            logger.warning(
                "sensitivity_tuning: %s not found in model; skipping.",
                prot_met_id,
            )
            break

        # Find the catalysed rxn consuming the most of this protein.
        # Consumption: S[prot, rxn] < 0 with positive flux -> negative product.
        # We pick the rxn with the most-negative S[prot] * flux.
        target_rxn = None
        target_value = 0.0
        for rxn in prot_met.reactions:
            if rxn.id == best_usage_id:
                continue
            coeff = rxn.metabolites[prot_met]
            flux = float(sol.fluxes.get(rxn.id, 0.0))
            value = coeff * flux
            if value < target_value:
                target_value = value
                target_rxn = rxn

        if target_rxn is None:
            logger.warning(
                "sensitivity_tuning: no catalysed reaction consuming %s; "
                "cannot tune.", prot_met_id,
            )
            break

        # Find the ec entry for this catalysed reaction.
        if target_rxn.id not in model.ec.rxns:
            logger.warning(
                "sensitivity_tuning: rxn %s consumes %s but is not in "
                "model.ec.rxns; cannot tune.",
                target_rxn.id, prot_met_id,
            )
            break
        ec_idx = model.ec.rxns.index(target_rxn.id)

        # Annotate notes (idempotently).
        if model.ec.source[ec_idx] != _TUNING_SOURCE:
            old_note = model.ec.notes[ec_idx]
            new_note = (
                f"preTuneKcat={model.ec.kcat[ec_idx]} | "
                f"source:{model.ec.source[ec_idx]}"
            )
            model.ec.notes[ec_idx] = (
                f"{old_note}; {new_note}" if old_note else new_note
            )

        model.ec.kcat[ec_idx] = float(model.ec.kcat[ec_idx]) * fold_change
        model.ec.source[ec_idx] = _TUNING_SOURCE
        apply_kcat_constraints(model, update_rxns=[target_rxn.id])

        if ec_idx not in tuned_ec_indices:
            tuned_ec_indices.append(ec_idx)

    return _build_result(model, tuned_ec_indices)


def _build_result(
    model: "EcModel", ec_indices: list[int],
) -> TunedKcatsResult:
    """Snapshot pre/post kcat values for tuned reactions."""
    if not ec_indices:
        return TunedKcatsResult()

    rxns = [model.ec.rxns[i] for i in ec_indices]
    rxn_names = [
        model.reactions.get_by_id(r).name or r for r in rxns
    ]
    new_kcat = np.array(
        [float(model.ec.kcat[i]) for i in ec_indices], dtype=float,
    )

    # Recover old_kcat from notes (the "preTuneKcat=X" annotation).
    old_kcat: list[float] = []
    sources: list[str] = []
    for i in ec_indices:
        note = model.ec.notes[i]
        if "preTuneKcat=" in note:
            tail = note.split("preTuneKcat=")[-1]
            num_str = tail.split(" ")[0]
            try:
                old_kcat.append(float(num_str))
            except ValueError:
                old_kcat.append(float("nan"))
            if "source:" in note:
                src = note.split("source:")[-1]
                sources.append(src.strip())
            else:
                sources.append("")
        else:
            old_kcat.append(float("nan"))
            sources.append("")

    enzymes: list[str] = []
    rxn_enz_dense = model.ec.rxn_enz_mat.toarray()
    for i in ec_indices:
        enz_idxs = np.where(rxn_enz_dense[i] != 0)[0]
        enzymes.append(";".join(model.ec.enzymes[j] for j in enz_idxs))

    return TunedKcatsResult(
        rxns=rxns,
        rxn_names=rxn_names,
        enzymes=enzymes,
        old_kcat=np.array(old_kcat, dtype=float),
        new_kcat=new_kcat,
        source=sources,
    )
