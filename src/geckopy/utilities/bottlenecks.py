"""Identify enzymes limiting the current objective.

Ported from geckopy_old/geckopy/flux_analysis.py:275-281
(get_protein_bottlenecks), adapted to use the Enzyme proxy.
"""
from __future__ import annotations

from typing import TYPE_CHECKING

import pandas as pd

if TYPE_CHECKING:
    from ..ec_model.ec_model import EcModel


def get_enzyme_bottlenecks(model: "EcModel", *, top: int = 10) -> pd.DataFrame:
    """Top-N enzymes by absolute shadow price.

    Solves the model, ranks enzymes by ``|shadow_price|``, returns a
    DataFrame indexed by uniprot id with columns
    ``[gene, shadow_price, flux, cap_usage, upper_bound]``.

    Raises ``RuntimeError`` if the solve is not optimal.

    Ported from geckopy_old/geckopy/flux_analysis.py:275-281.
    """
    sol = model.optimize()
    if sol.status != "optimal":
        raise RuntimeError(
            f"Solver status {sol.status!r}; cannot rank bottlenecks"
        )

    rows = []
    for enz in model.enzymes:
        rows.append({
            "uniprot": enz.id,
            "gene": enz.gene,
            "shadow_price": enz.shadow_price,
            "flux": enz.flux,
            "cap_usage": enz.cap_usage,
            "upper_bound": enz.upper_bound,
        })
    df = pd.DataFrame(rows).set_index("uniprot")
    df = df.reindex(
        df["shadow_price"].abs().sort_values(ascending=False).index
    )
    return df.head(top)
