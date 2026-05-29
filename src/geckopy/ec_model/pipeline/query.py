"""Query helpers over the ec substructure.

These functions do not mutate the model; they pull subsets of the ec
data for inspection, reporting, or downstream pipeline use.
"""
from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    from ..ec_model import EcModel


def get_reactions_from_enzyme(
    model: "EcModel", protein_id: str
) -> pd.DataFrame:
    """Return all reactions catalyzed by a given enzyme, as a DataFrame.

    Ported from GECKO MATLAB:
    src/geckomat/change_model/getReactionsFromEnzyme.m.

    MATLAB-COMPAT: MATLAB returns five separate outputs (rxns, kcat,
    idx, rxnNames, grRules). geckopy returns a single pandas DataFrame
    with columns ``rxn_id``, ``kcat``, ``name``, ``gpr``. The
    ``idx`` output is dropped because it is trivially derivable as
    ``ec.rxns.index(rxn_id)``. To get just the reaction IDs:
    ``df['rxn_id'].tolist()``. To sort by kcat:
    ``df.sort_values('kcat')``.

    MATLAB-COMPAT: MATLAB matches case-insensitively. geckopy is
    case-sensitive throughout.

    MATLAB-COMPAT: MATLAB returns empty outputs when the protein ID is
    unknown. geckopy raises ValueError.

    Parameters
    ----------
    model
        An EcModel with populated ec.enzymes and ec.rxn_enz_mat.
    protein_id
        UniProt accession (or whatever IDs are in ec.enzymes),
        case-sensitive.

    Returns
    -------
    pandas.DataFrame
        One row per catalyzed reaction, with columns:

        - ``rxn_id``: the reaction ID.
        - ``kcat``: the kcat value (1/s) from ec.kcat (0 if unset).
        - ``name``: the reaction name from the cobra model.
        - ``gpr``: the gene-protein-reaction rule.

        Empty DataFrame with the right columns if the enzyme catalyzes
        no reactions in the current model.

    Raises
    ------
    ValueError
        If ``protein_id`` is not found in ec.enzymes.
    """
    try:
        prot_idx = model.ec.enzymes.index(protein_id)
    except ValueError:
        raise ValueError(
            f"protein_id '{protein_id}' not found in ec.enzymes."
        ) from None

    col = model.ec.rxn_enz_mat.tocsc().getcol(prot_idx)
    ec_rxn_indices = np.array(col.nonzero()[0], dtype=int)

    rxn_ids = [model.ec.rxns[i] for i in ec_rxn_indices]
    kcats = model.ec.kcat[ec_rxn_indices].tolist()

    names: list[str] = []
    gprs: list[str] = []
    for rxn_id in rxn_ids:
        rxn = model.reactions.get_by_id(rxn_id)
        names.append(rxn.name or "")
        gprs.append(rxn.gene_reaction_rule)

    return pd.DataFrame({
        "rxn_id": rxn_ids,
        "kcat": kcats,
        "name": names,
        "gpr": gprs,
    })
