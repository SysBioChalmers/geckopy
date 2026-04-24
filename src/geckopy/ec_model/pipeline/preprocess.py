"""Structural preprocessing of the GEM before enzyme extension.

Corresponds to stages 1 to 4 of makeEcModel in GECKO MATLAB. These
stages do not touch the `ec` substructure; they only reshape the
underlying model so subsequent stages operate on a clean irreversible form.
"""
from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import cobra
    from ...adapter import ModelAdapter


def remove_pseudoreaction_gprs(
    model: "cobra.Model", adapter: "ModelAdapter"
) -> list[str]:
    """Clear gene-protein-reaction rules from pseudoreactions.

    Corresponds to stage 1 of GECKO MATLAB `makeEcModel`. A reaction is
    treated as a pseudoreaction if either:

    - its name contains the substring ``pseudoreaction`` (case-sensitive,
      matching MATLAB's default `contains` behavior); or
    - its ID is listed in the first column of ``data/pseudoRxns.tsv``
      inside the adapter's project folder, if that file exists.

    Parameters
    ----------
    model
        A cobra.Model that will be mutated in place: the gene reaction
        rule of each matched reaction is cleared (empty string). This
        automatically updates the reaction's `.genes` set via cobrapy's
        GPR handling.
    adapter
        Provides the project path used to locate the optional TSV.

    Returns
    -------
    list of str
        Sorted IDs of reactions whose GPR was cleared. Intended for
        logging and test assertions.
    """
    pseudo_ids: set[str] = set()

    # By reaction name (MATLAB: contains(model.rxnNames,'pseudoreaction'))
    for rxn in model.reactions:
        if rxn.name and "pseudoreaction" in rxn.name:
            pseudo_ids.add(rxn.id)

    # By optional TSV (first column = reaction ID)
    tsv_path = adapter.params.path / "data" / "pseudoRxns.tsv"
    if tsv_path.is_file():
        with open(tsv_path, "r", encoding="utf-8") as f:
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if parts and parts[0]:
                    pseudo_ids.add(parts[0])

    # Clear GPRs for matched reactions that exist in the model.
    model_rxn_ids = {r.id for r in model.reactions}
    cleared: list[str] = []
    for rxn_id in pseudo_ids & model_rxn_ids:
        rxn = model.reactions.get_by_id(rxn_id)
        rxn.gene_reaction_rule = ""
        cleared.append(rxn_id)

    return sorted(cleared)