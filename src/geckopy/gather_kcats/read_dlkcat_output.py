"""Read a DLKcat predictions TSV into a kcat_list DataFrame.

Ported from GECKO MATLAB:
src/geckomat/gather_kcats/readDLKcatOutput.m.
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING

import pandas as pd

if TYPE_CHECKING:
    from ..ec_model.ec_model import EcModel

logger = logging.getLogger(__name__)


_OUTPUT_COLUMNS = [
    "rxn_id",
    "source",
    "eccode",
    "substrates",
    "genes",
    "kcat",
    "wildcard_level",
    "origin",
]


def read_dlkcat_output(
    model: "EcModel",
    file_path: str | Path,
) -> pd.DataFrame:
    """Parse DLKcat's output TSV into a kcat_list DataFrame.

    Ported from GECKO MATLAB:
    src/geckomat/gather_kcats/readDLKcatOutput.m.

    The input file is the 6-column TSV that DLKcat writes (see
    ``run_dlkcat``), with one header line and columns
    ``rxn_id, gene, substrate, smiles, sequence, kcat``. Rows whose
    kcat column is non-numeric (``"NA"``, empty, etc.) are silently
    dropped before validation against the model.

    The output DataFrame matches the schema produced by
    ``fuzzy_kcat_matching`` so downstream functions
    (``apply_kcat_list``, ``merge_dlkcat_and_fuzzy_kcats``) can
    consume both interchangeably. Columns not produced by DLKcat
    (``eccode``, ``wildcard_level``, ``origin``) are filled with
    ``""`` / ``<NA>``.

    MATLAB-COMPAT: GECKO MATLAB takes a ``modelAdapter`` and defaults
    the path to ``adapter.params.path/data/DLKcat.tsv``. geckopy
    requires ``file_path`` explicitly; the caller resolves it.

    MATLAB-COMPAT: GECKO MATLAB's substrate-name match against
    ``model.metNames`` is case-sensitive (``ismember`` default).
    geckopy uses case-insensitive matching, which is more lenient
    and avoids false-positive failures when SBML loaders differ in
    capitalization. Tracked as a MATLAB-side improvement in
    ``docs/future_improvements.md``.

    Parameters
    ----------
    model
        EcModel with populated ``ec.rxns`` and ``metabolites``.
    file_path
        Path to the DLKcat output TSV.

    Returns
    -------
    pandas.DataFrame
        Rows where the kcat field was numeric, with the columns
        listed above. ``source`` is ``"DLKcat"`` for every row.

    Raises
    ------
    FileNotFoundError
        If ``file_path`` does not exist.
    ValueError
        If the file contains no numeric kcat values; if any
        substrate name from the file is not in
        ``model.metabolites.<name>`` (case-insensitive); or if any
        ``rxn_id`` is not in ``model.ec.rxns``.
    """
    file_path = Path(file_path)
    if not file_path.is_file():
        raise FileNotFoundError(f"DLKcat output file not found: {file_path}")

    df = pd.read_csv(
        file_path,
        sep="\t",
        header=0,
        names=["rxn_id", "gene", "substrate", "smiles", "sequence", "kcat"],
        dtype=str,
        keep_default_na=False,
    )

    numeric_kcat = pd.to_numeric(df["kcat"], errors="coerce")
    valid_mask = numeric_kcat.notna()

    if not valid_mask.any():
        raise ValueError(
            f"{file_path} contains no numeric kcat values. Did you forget "
            f"to run run_dlkcat() first?"
        )

    model_met_names_lower = {m.name.lower() for m in model.metabolites}
    file_subs_lower = df["substrate"].str.lower()
    unknown_subs_mask = ~file_subs_lower.isin(model_met_names_lower)
    if unknown_subs_mask.any():
        unknown = sorted(df.loc[unknown_subs_mask, "substrate"].unique())
        preview = unknown[:5]
        raise ValueError(
            f"DLKcat output references {len(unknown)} substrate name(s) "
            f"not present in model.metabolites (case-insensitive). "
            f"Examples: {preview}. The output was likely generated from a "
            f"different ecModel."
        )

    ec_rxn_ids = set(model.ec.rxns)
    unknown_rxns_mask = ~df["rxn_id"].isin(ec_rxn_ids)
    if unknown_rxns_mask.any():
        unknown = sorted(df.loc[unknown_rxns_mask, "rxn_id"].unique())
        preview = unknown[:5]
        raise ValueError(
            f"DLKcat output references {len(unknown)} reaction ID(s) "
            f"not present in model.ec.rxns. Examples: {preview}. The "
            f"output was likely generated from a different ecModel."
        )

    kept = df[valid_mask].reset_index(drop=True)
    kept_kcats = numeric_kcat[valid_mask].reset_index(drop=True).astype(float)
    n = len(kept)

    out = pd.DataFrame({
        "rxn_id": kept["rxn_id"].astype(str).tolist(),
        "source": ["DLKcat"] * n,
        "eccode": [""] * n,
        "substrates": [[s] for s in kept["substrate"].astype(str)],
        "genes": [[g] for g in kept["gene"].astype(str)],
        "kcat": kept_kcats.tolist(),
        "wildcard_level": pd.array([pd.NA] * n, dtype="Int64"),
        "origin": pd.array([pd.NA] * n, dtype="Int64"),
    })

    dropped = len(df) - n
    if dropped > 0:
        logger.info(
            "read_dlkcat_output: read %d row(s); dropped %d row(s) with "
            "non-numeric kcat.",
            n, dropped,
        )
    else:
        logger.info("read_dlkcat_output: read %d row(s) from %s.", n, file_path)

    return out
