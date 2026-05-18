"""Flux-Scanning with Enforced Objective Function for an ecModel.

Ported from GECKO MATLAB:
src/geckomat/utilities/ecFSEOF.m.
"""
from __future__ import annotations

import logging
from collections import defaultdict
from dataclasses import dataclass, field
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    from ..ec_model.ec_model import EcModel


logger = logging.getLogger(__name__)


from ..ec_model.constants import (
    PROT_PREFIX,
    USAGE_PREFIX,
)

_STANDARD_GPR = "standard"
_BIOMASS_INFEASIBLE_THRESHOLD = 1e-8
_TOP_QUANTILE = 0.75
_MAX_TARGET_FRACTION = 0.9


_RXN_TARGET_COLUMNS = [
    "rxn_id", "rxn_name", "slope", "gpr", "equation", "action",
]
_GENE_TARGET_COLUMNS = [
    "gene_id", "gene_name", "slope", "action", "essentiality",
]


@dataclass
class EcFseofResult:
    """Output of `ec_fseof`.

    Attributes
    ----------
    alpha
        Enforced production-target flux values used for the scan.
    v_matrix
        Per-target-rxn fluxes at each alpha. Index = rxn_id; columns
        are ``str(alpha[k])``. Excludes ``usage_prot_*`` rxns.
    rxn_targets
        DataFrame of reaction targets (excluding transport),
        sorted by descending slope. Columns: ``rxn_id``,
        ``rxn_name``, ``slope``, ``gpr``, ``equation``, ``action``.
    transport_targets
        Same shape as ``rxn_targets`` but for transport rxns
        (the same metabolite name appearing in multiple
        compartments).
    gene_targets
        DataFrame of gene targets, sorted by descending mean slope.
        Columns: ``gene_id``, ``gene_name``, ``slope``, ``action``,
        ``essentiality``.
    """

    alpha: np.ndarray = field(
        default_factory=lambda: np.empty(0, dtype=float)
    )
    v_matrix: pd.DataFrame = field(default_factory=pd.DataFrame)
    rxn_targets: pd.DataFrame = field(
        default_factory=lambda: pd.DataFrame(columns=_RXN_TARGET_COLUMNS)
    )
    transport_targets: pd.DataFrame = field(
        default_factory=lambda: pd.DataFrame(columns=_RXN_TARGET_COLUMNS)
    )
    gene_targets: pd.DataFrame = field(
        default_factory=lambda: pd.DataFrame(columns=_GENE_TARGET_COLUMNS)
    )


def ec_fseof(
    model: "EcModel",
    prod_target_rxn: str,
    cs_rxn: str,
    *,
    n_steps: int = 16,
) -> EcFseofResult:
    """Run Flux-Scanning with Enforced Objective Function (FSEOF).

    For each of `n_steps` alpha values between the initial production
    flux at biomass-max and 90% of the maximum theoretical production,
    fixes the production-target reaction to alpha and maximises
    biomass. Reactions whose flux changes monotonically with alpha
    are returned as engineering targets (OE, KD, KO).

    Ported from GECKO MATLAB:
    src/geckomat/utilities/ecFSEOF.m.

    MATLAB-COMPAT: GECKO MATLAB optionally writes the result tables
    to TSV files. geckopy returns the tables as DataFrames; callers
    can `.to_csv()` them themselves.

    MATLAB-COMPAT: GECKO MATLAB computes slope as
    ``(v[nSteps-1] - v[1]) / (maxTarget - maxTarget/nSteps - 1)``
    (1-indexed, with an oddly-formed denominator that subtracts 1
    from `maxTarget*(1 - 1/nSteps)`). geckopy uses the cleaner
    ``(v[-1] - v[0]) / (alpha[-1] - alpha[0])``.

    MATLAB-COMPAT: The MATLAB `standardizeGrRules` preprocessing is
    RAVEN-specific. cobrapy handles GPRs natively, so geckopy skips
    that step.

    Parameters
    ----------
    model
        EcModel with a populated ec substructure and the protein
        pool/usage rxn machinery.
    prod_target_rxn
        Reaction id for the production target (typically an
        exchange).
    cs_rxn
        Reaction id for the main carbon source uptake (used only for
        a sanity-check warning).
    n_steps
        Number of alpha values in the scan.

    Returns
    -------
    EcFseofResult
    """
    from ..adapter import resolve_adapter
    adapter = resolve_adapter(
        model,
        purpose="ec_fseof reads params.bio_rxn from the adapter",
    )
    if n_steps < 2:
        raise ValueError(f"n_steps must be >= 2, got {n_steps}")

    bio_rxn_id = adapter.params.bio_rxn
    bio_rxn = model.reactions.get_by_id(bio_rxn_id)
    prod_rxn = model.reactions.get_by_id(prod_target_rxn)
    cs_rxn_obj = model.reactions.get_by_id(cs_rxn)

    # Step 1: solve at biomass max, capture initial production flux
    # and check carbon source consistency.
    with model:
        model.objective = bio_rxn_id
        model.objective_direction = "max"
        sol = model.optimize()
        if sol.fluxes is None:
            raise ValueError("Initial biomass-max LP is infeasible.")
        ini_target = float(sol.fluxes[prod_target_rxn])
        cs_flux = float(sol.fluxes[cs_rxn])
        if cs_rxn_obj.lower_bound < cs_flux:
            logger.warning(
                "ec_fseof: carbon source lower bound is %g but "
                "uptake at biomass-max is %g; consider tightening.",
                cs_rxn_obj.lower_bound, cs_flux,
            )

    # Step 2: solve at production-target max; cap at 90%.
    with model:
        model.objective = prod_target_rxn
        model.objective_direction = "max"
        sol = model.optimize()
        if sol.fluxes is None:
            raise ValueError("Production-target-max LP is infeasible.")
        max_target = float(sol.fluxes[prod_target_rxn]) * _MAX_TARGET_FRACTION

    if max_target <= ini_target:
        raise ValueError(
            f"max_target ({max_target}) <= ini_target ({ini_target}); "
            f"FSEOF requires headroom. Check that prod_target_rxn can "
            f"actually carry more flux at lower biomass."
        )

    alpha = np.linspace(ini_target, max_target, n_steps)

    # Step 3: enforce alpha, maximise biomass for each step.
    rxn_ids = [r.id for r in model.reactions]
    n_rxns = len(rxn_ids)
    v_matrix = np.zeros((n_rxns, n_steps), dtype=float)
    for k, a in enumerate(alpha):
        with model:
            prod_rxn.lower_bound = float(a)
            prod_rxn.upper_bound = float(a)
            bio_rxn.lower_bound = 0.0
            model.objective = bio_rxn_id
            model.objective_direction = "max"
            sol = model.optimize()
            if sol.fluxes is None:
                v_matrix[:, k] = np.nan
            else:
                for i, rid in enumerate(rxn_ids):
                    v_matrix[i, k] = float(sol.fluxes[rid])

    # Step 4: filter rxns (have GPR, not "standard", any non-zero flux).
    has_gpr = np.array([
        bool(model.reactions[i].gene_reaction_rule)
        and _STANDARD_GPR not in model.reactions[i].gene_reaction_rule
        for i in range(n_rxns)
    ])
    abs_v = np.abs(v_matrix)
    nonzero = abs_v.max(axis=1) > 0
    keep_mask = has_gpr & nonzero
    if not keep_mask.any():
        return EcFseofResult(alpha=alpha)

    kept_indices = np.where(keep_mask)[0]
    kept_rxn_ids = [rxn_ids[i] for i in kept_indices]
    kept_v = v_matrix[kept_indices, :]
    kept_abs = abs_v[kept_indices, :]

    # Step 5: monotonic patterns + slope.
    actions: list[str] = []
    slopes: list[float] = []
    selected: list[bool] = []
    denom = abs(alpha[-1] - alpha[0])
    for i, rid in enumerate(kept_rxn_ids):
        diffs = np.diff(kept_abs[i])
        if np.all(diffs > 0):  # strictly ascending |flux|
            selected.append(True)
            actions.append("OE")
        elif np.all(diffs < 0):  # strictly descending |flux|
            selected.append(True)
            actions.append("KO" if kept_v[i, -1] == 0 else "KD")
        else:
            selected.append(False)
            actions.append("")
        slopes.append(
            abs(kept_v[i, -1] - kept_v[i, 0]) / denom if denom > 0 else 0.0
        )

    sel_arr = np.array(selected)
    if not sel_arr.any():
        return EcFseofResult(alpha=alpha)

    sel_indices = np.where(sel_arr)[0]
    target_rxn_ids = [kept_rxn_ids[i] for i in sel_indices]
    target_v = kept_v[sel_indices, :]
    target_slopes = np.array([slopes[i] for i in sel_indices])
    target_actions = [actions[i] for i in sel_indices]

    # Step 6: top 25% by slope.
    threshold = float(np.quantile(target_slopes, _TOP_QUANTILE))
    quant_mask = target_slopes > threshold
    if not quant_mask.any():
        return EcFseofResult(alpha=alpha)

    target_rxn_ids = [r for r, k in zip(target_rxn_ids, quant_mask) if k]
    target_v = target_v[quant_mask, :]
    target_slopes = target_slopes[quant_mask]
    target_actions = [a for a, k in zip(target_actions, quant_mask) if k]

    # Sort by descending slope.
    order = np.argsort(-target_slopes)
    target_rxn_ids = [target_rxn_ids[i] for i in order]
    target_v = target_v[order, :]
    target_slopes = target_slopes[order]
    target_actions = [target_actions[i] for i in order]

    # Step 7: gene-level analysis.
    gene_to_rxns: dict[str, list[int]] = defaultdict(list)
    for rxn_idx_in_targets, rid in enumerate(target_rxn_ids):
        rxn = model.reactions.get_by_id(rid)
        for gene in rxn.genes:
            if gene.id == _STANDARD_GPR:
                continue
            gene_to_rxns[gene.id].append(rxn_idx_in_targets)

    gene_targets_rows: list[dict] = []
    for gene_id, rxn_indices in gene_to_rxns.items():
        actions_for_gene = sorted({target_actions[i] for i in rxn_indices})
        essentiality = _check_essentiality(model, gene_id, bio_rxn_id)

        if len(actions_for_gene) > 1:
            if "OE" in actions_for_gene:
                action = "OE"
            elif essentiality != "essential":
                action = "KO"
            else:
                action = "KD"
        else:
            action = actions_for_gene[0]

        slope_values = [target_slopes[i] for i in rxn_indices]
        gene = model.genes.get_by_id(gene_id)
        gene_name = gene.name or ""
        gene_targets_rows.append({
            "gene_id": gene_id,
            "gene_name": gene_name,
            "slope": float(np.mean(slope_values)),
            "action": action,
            "essentiality": essentiality,
        })

    gene_targets = pd.DataFrame(
        gene_targets_rows,
        columns=["gene_id", "gene_name", "slope", "action", "essentiality"],
    ).sort_values("slope", ascending=False).reset_index(drop=True)

    # Step 8: build rxn_targets and transport_targets DataFrames.
    keep_for_rxn_table = [
        not rid.startswith(USAGE_PREFIX) for rid in target_rxn_ids
    ]
    final_rxn_ids = [
        r for r, k in zip(target_rxn_ids, keep_for_rxn_table) if k
    ]
    final_v = target_v[keep_for_rxn_table, :]
    final_slopes = target_slopes[keep_for_rxn_table]
    final_actions = [
        a for a, k in zip(target_actions, keep_for_rxn_table) if k
    ]

    rxn_rows: list[dict] = []
    transport_flags: list[bool] = []
    for rid, slope, action in zip(final_rxn_ids, final_slopes, final_actions):
        rxn = model.reactions.get_by_id(rid)
        non_prot_met_names = [
            m.name for m in rxn.metabolites
            if not m.id.startswith(PROT_PREFIX)
        ]
        is_transport = (
            len(non_prot_met_names) != len(set(non_prot_met_names))
        )
        transport_flags.append(is_transport)
        rxn_rows.append({
            "rxn_id": rid,
            "rxn_name": rxn.name or "",
            "slope": float(slope),
            "gpr": rxn.gene_reaction_rule,
            "equation": rxn.build_reaction_string(),
            "action": action,
        })

    cols = ["rxn_id", "rxn_name", "slope", "gpr", "equation", "action"]
    all_rxn_df = pd.DataFrame(rxn_rows, columns=cols)
    transport_df = all_rxn_df[transport_flags].reset_index(drop=True)
    rxn_targets_df = all_rxn_df[
        [not f for f in transport_flags]
    ].reset_index(drop=True)

    # v_matrix DataFrame (excluding usage_prot_* rxns).
    v_matrix_df = pd.DataFrame(
        final_v,
        index=pd.Index(final_rxn_ids, name="rxn_id"),
        columns=[str(a) for a in alpha],
    )

    return EcFseofResult(
        alpha=alpha,
        v_matrix=v_matrix_df,
        rxn_targets=rxn_targets_df,
        transport_targets=transport_df,
        gene_targets=gene_targets,
    )


def _check_essentiality(
    model: "EcModel", gene_id: str, bio_rxn_id: str,
) -> str:
    """Block all usage_prot_<X> rxns for enzymes that map to this
    gene; if biomass collapses, the gene is essential."""
    if gene_id not in model.ec.genes:
        return ""
    enz_indices = [
        i for i, g in enumerate(model.ec.genes) if g == gene_id
    ]
    usage_rxn_ids = [
        f"{USAGE_PREFIX}{model.ec.enzymes[i]}" for i in enz_indices
    ]
    cobra_rxn_ids = {r.id for r in model.reactions}
    usage_rxn_ids = [r for r in usage_rxn_ids if r in cobra_rxn_ids]
    if not usage_rxn_ids:
        return ""
    with model:
        for rid in usage_rxn_ids:
            r = model.reactions.get_by_id(rid)
            r.lower_bound = 0.0
            r.upper_bound = 0.0
        model.objective = bio_rxn_id
        sol = model.optimize()
        if sol.fluxes is None:
            return "essential"
        bio_flux = float(sol.fluxes.get(bio_rxn_id, 0.0))
        if bio_flux < _BIOMASS_INFEASIBLE_THRESHOLD:
            return "essential"
    return ""
