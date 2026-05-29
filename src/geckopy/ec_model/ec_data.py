"""Enzyme-constraint data container for ecModels.

Mirrors the `model.ec` substructure in GECKO MATLAB. Two groups of fields:

Per-reaction fields (length N_rxns, one entry per catalyzed reaction or
isozyme-specific reaction):
    rxns, kcat, source, notes, eccodes

Per-enzyme fields (length N_enzymes, one entry per unique enzyme):
    genes, enzymes, mw, sequence, concs

Coupling between them (shape N_rxns x N_enzymes):
    rxn_enz_mat
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np
from scipy import sparse


@dataclass
class EcData:
    """Enzyme-constraint data attached to an EcModel.

    Full vs light models differ in the semantics of `rxns`:

    Full model
        Model reactions are duplicated per isozyme by expandModel, and
        `rxns` contains one entry per catalyzed reaction. Length of
        `rxns` equals the number of such reactions in the host model.

    Light model
        Model reactions are not duplicated. `rxns` contains duplicate
        entries per isozyme, ordered the same as model.reactions.

    Per-enzyme fields and rxn_enz_mat have identical semantics in both.
    """

    gecko_light: bool = False

    # Per-reaction fields (length N_rxns)
    rxns: list[str] = field(default_factory=list)
    kcat: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=float))
    source: list[str] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)
    eccodes: list[str] = field(default_factory=list)

    # Per-enzyme fields (length N_enzymes)
    genes: list[str] = field(default_factory=list)
    enzymes: list[str] = field(default_factory=list)
    mw: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=float))
    sequence: list[str] = field(default_factory=list)
    concs: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=float))

    # Coupling (N_rxns x N_enzymes)
    rxn_enz_mat: sparse.csr_matrix = field(
        default_factory=lambda: sparse.csr_matrix((0, 0), dtype=float)
    )

    @property
    def n_rxns(self) -> int:
        return len(self.rxns)

    @property
    def n_enzymes(self) -> int:
        return len(self.enzymes)

    def validate(self) -> None:
        """Raise ValueError if internal field lengths are inconsistent."""
        n_r, n_e = self.n_rxns, self.n_enzymes

        rxn_lengths = {
            "kcat": len(self.kcat),
            "source": len(self.source),
            "notes": len(self.notes),
            "eccodes": len(self.eccodes),
        }
        for name, length in rxn_lengths.items():
            if length != n_r:
                raise ValueError(
                    f"ec.{name} has length {length}, expected {n_r} "
                    f"(matching ec.rxns)"
                )

        enz_lengths = {
            "enzymes": len(self.enzymes),
            "mw": len(self.mw),
            "sequence": len(self.sequence),
            "concs": len(self.concs),
        }
        for name, length in enz_lengths.items():
            if length != n_e:
                raise ValueError(
                    f"ec.{name} has length {length}, expected {n_e} "
                    f"(matching ec.genes)"
                )

        if self.rxn_enz_mat.shape != (n_r, n_e):
            raise ValueError(
                f"ec.rxn_enz_mat has shape {self.rxn_enz_mat.shape}, "
                f"expected ({n_r}, {n_e})"
            )

    @staticmethod
    def empty_for_reactions(n_rxns: int, n_enzymes: int = 0, *,
                            gecko_light: bool = False) -> "EcData":
        """Create an EcData with empty strings and NaN arrays of given size.

        Used by makeEcModel step 6 to preallocate the structure before
        populating it.
        """
        return EcData(
            gecko_light=gecko_light,
            rxns=[""] * n_rxns,
            kcat=np.full(n_rxns, np.nan, dtype=float),
            source=[""] * n_rxns,
            notes=[""] * n_rxns,
            eccodes=[""] * n_rxns,
            genes=[""] * n_enzymes,
            enzymes=[""] * n_enzymes,
            mw=np.full(n_enzymes, np.nan, dtype=float),
            sequence=[""] * n_enzymes,
            concs=np.full(n_enzymes, np.nan, dtype=float),
            rxn_enz_mat=sparse.lil_matrix((n_rxns, n_enzymes), dtype=float).tocsr(),
        )