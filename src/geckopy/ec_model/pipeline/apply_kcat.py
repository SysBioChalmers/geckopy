"""Apply kcat-derived stoichiometric constraints to an ecModel.

Ported from GECKO MATLAB: src/geckomat/change_model/applyKcatConstraints.m.
Only the full (non-light) formulation is supported; the gecko-light
branch raises NotImplementedError.
"""
from __future__ import annotations

import warnings
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from ..ec_model import EcModel


from ..constants import POOL_ID, PROT_PREFIX


def apply_kcat_constraints(
    model: "EcModel",
    update_rxns: list[str] | None = None,
) -> None:
    """Translate ec.kcat values into stoichiometric coefficients.

    Ported from GECKO MATLAB: src/geckomat/change_model/applyKcatConstraints.m
    (full-model branch only).

    For each reaction in ``ec.rxns`` with a valid kcat and at least one
    associated enzyme, writes a negative coefficient at
    ``S[prot_<enzyme>, rxn]`` of magnitude::

        subunits * MW / (kcat * 3600)

    where ``subunits`` is the corresponding entry of ``ec.rxn_enz_mat``
    (1 for a single subunit, higher for multi-copy complexes), ``MW``
    is the enzyme molecular weight from ``ec.mw`` (in Da), and the
    factor 3600 converts kcat from 1/s (its conventional unit and the
    unit used in ``ec.kcat``) to 1/h so it aligns with cobrapy's
    mmol/gDW/h flux convention.

    Entries with kcat == 0 or kcat == NaN are treated as "no kcat
    assigned"; their stoichiometric coefficients are set to 0 (i.e.,
    the enzyme places no constraint on that reaction).

    MATLAB-COMPAT: GECKO MATLAB clears existing coefficients before
    checking kcat validity, so flipping a kcat to NaN or 0 and
    re-applying clears the prior constraint. geckopy matches this
    semantics exactly. (An earlier draft of geckopy did the check
    first and skipped clearing on all-invalid; this was changed to
    match MATLAB.) No MATLAB-side change required.

    The function is idempotent: existing non-zero coefficients at
    ``S[prot_<enzyme>, rxn]`` are first cleared for the reactions being
    updated, then freshly written from current ec.kcat. Running it
    twice yields the same result as running it once. Setting a kcat
    back to NaN/0 and re-applying clears the previously written
    coefficient.

    Parameters
    ----------
    model
        An EcModel produced by ``make_ec_model``, with at least some
        kcat values assigned to ``model.ec.kcat``. Mutated in place.
    update_rxns
        Reaction IDs (from ``model.ec.rxns``) to update. If None,
        updates all entries.

    Raises
    ------
    NotImplementedError
        If the model is gecko-light (not yet supported).
    ValueError
        If ``update_rxns`` contains IDs not present in ``ec.rxns``.

    Warns
    -----
    UserWarning
        If every kcat to be applied is 0 or NaN. Old coefficients are
        still cleared in this case; nothing is written afterward.
    """
    if model.ec.gecko_light:
        raise NotImplementedError(
            "apply_kcat_constraints: gecko-light formulation is not yet "
            "implemented."
        )

    ec_rxn_ids = model.ec.rxns

    # Resolve update_rxns to a numpy array of ec-row indices.
    if update_rxns is None:
        selected_idx = np.arange(len(ec_rxn_ids))
    else:
        id_to_idx = {rxn_id: i for i, rxn_id in enumerate(ec_rxn_ids)}
        unknown = [r for r in update_rxns if r not in id_to_idx]
        if unknown:
            raise ValueError(
                f"update_rxns contains IDs not present in ec.rxns: {unknown[:5]}"
            )
        selected_idx = np.array([id_to_idx[r] for r in update_rxns], dtype=int)

    if selected_idx.size == 0:
        return

    # Step 1: clear existing prot_<enzyme> coefficients on every reaction
    # in the update set. We never touch prot_pool here. This step always
    # runs regardless of kcat validity, so that flipping a kcat to NaN
    # and re-applying genuinely clears the old constraint.
    prot_met_ids = {
        m.id for m in model.metabolites
        if m.id.startswith(PROT_PREFIX) and m.id != POOL_ID
    }

    for idx in selected_idx:
        rxn_id = ec_rxn_ids[idx]
        rxn = model.reactions.get_by_id(rxn_id)
        to_clear = {
            m: 0.0
            for m in rxn.metabolites
            if m.id in prot_met_ids
        }
        if to_clear:
            rxn.add_metabolites(to_clear, combine=False)

    # Step 2: identify which selected entries have a real kcat (>0).
    kcats_selected = model.ec.kcat[selected_idx]
    valid = kcats_selected > 0

    if not valid.any():
        warnings.warn(
            "apply_kcat_constraints: ec.kcat has no real entries for the "
            "selected reactions; existing constraints were cleared and no "
            "new ones written.",
            UserWarning,
            stacklevel=2,
        )
        return

    # Step 3: write fresh coefficients for valid kcats.
    enz_idx_to_met_id = {
        i: f"{PROT_PREFIX}{enz}"
        for i, enz in enumerate(model.ec.enzymes)
    }
    rxn_enz_mat_csr = model.ec.rxn_enz_mat.tocsr()

    for local_i, idx in enumerate(selected_idx):
        if not valid[local_i]:
            continue

        kcat = float(kcats_selected[local_i])
        rxn_id = ec_rxn_ids[idx]
        rxn = model.reactions.get_by_id(rxn_id)

        row = rxn_enz_mat_csr.getrow(idx)
        enzyme_indices = row.indices
        subunit_counts = row.data

        if enzyme_indices.size == 0:
            continue

        updates: dict = {}
        for enz_idx, subunits in zip(enzyme_indices, subunit_counts):
            mw = float(model.ec.mw[enz_idx])
            coef = -subunits * mw / (kcat * 3600.0)
            met = model.metabolites.get_by_id(enz_idx_to_met_id[enz_idx])
            updates[met] = coef

        if updates:
            rxn.add_metabolites(updates, combine=False)
