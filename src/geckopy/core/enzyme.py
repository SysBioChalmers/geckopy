"""Define the Enzyme class for enzyme-constrained models."""

from copy import deepcopy
from typing import TYPE_CHECKING, FrozenSet, Optional

from cobra.core.object import Object


if TYPE_CHECKING:
    from geckopy.core.ec_model import ECModel


class Enzyme(Object):
    """Represents an enzyme (protein) in an enzyme-constrained metabolic model.

    An Enzyme corresponds to a single protein identified by a database
    accession (e.g. UniProt ID). It stores biophysical properties used to
    build enzyme-capacity constraints in an ECModel.

    Parameters
    ----------
    id : str, optional
        Database accession for the protein (e.g. UniProt ID such as "P00549").
    name : str, optional
        Human-readable protein name.
    mw : float, optional
        Molecular weight in g/mol.
    sequence : str, optional
        Amino acid sequence (single-letter code).
    concentration : float, optional
        Measured enzyme concentration in mmol/gDCW (e.g. from proteomics).
        None if no measurement is available.

    Attributes
    ----------
    mw : float or None
        Molecular weight [g/mol].
    sequence : str or None
        Amino acid sequence.
    concentration : float or None
        Measured concentration [mmol/gDCW]. None when not constrained.
    """

    def __init__(
        self,
        id: Optional[str] = None,
        name: str = "",
        mw: Optional[float] = None,
        sequence: Optional[str] = None,
        concentration: Optional[float] = None,
    ) -> None:
        """Initialize an Enzyme.

        Parameters
        ----------
        id : str, optional
            Database accession (e.g. UniProt ID).
        name : str, optional
            Human-readable protein name.
        mw : float, optional
            Molecular weight in g/mol.
        sequence : str, optional
            Amino acid sequence (single-letter code).
        concentration : float, optional
            Measured enzyme concentration in mmol/gDCW.
        """
        super().__init__(id=id, name=name)
        self.mw = mw
        self.sequence = sequence
        self.concentration = concentration
        self._model: Optional["ECModel"] = None
        self._reaction: set = set()

    @property
    def reactions(self) -> FrozenSet:
        """Return a frozenset of reactions constrained by this enzyme.

        Returns
        -------
        FrozenSet
            Reactions whose flux is limited by the capacity of this enzyme.
        """
        return frozenset(self._reaction)

    @property
    def model(self) -> Optional["ECModel"]:
        """Return the ECModel this enzyme belongs to.

        Returns
        -------
        ECModel or None
        """
        return self._model

    @property
    def pseudometabolite_id(self) -> Optional[str]:
        """Return the ID of the enzyme pseudometabolite (``prot_<id>``).

        Returns
        -------
        str or None
        """
        if self.id is None:
            return None
        return f"prot_{self.id}"

    @property
    def usage_reaction_id(self) -> Optional[str]:
        """Return the ID of the enzyme usage pseudoreaction (``usage_prot_<id>``).

        Returns
        -------
        str or None
        """
        if self.id is None:
            return None
        return f"usage_prot_{self.id}"

    def __getstate__(self) -> dict:
        """Return serializable state, dropping the back-reference to the model."""
        state = super().__getstate__()
        state["_reaction"] = set()
        return state

    def copy(self) -> "Enzyme":
        """Return a deep copy of this enzyme (not associated with any model).

        Returns
        -------
        Enzyme
        """
        return deepcopy(self)

    def __repr__(self) -> str:
        mw_str = f", mw={self.mw}" if self.mw is not None else ""
        return f"<Enzyme {self.id}{mw_str} at {id(self):#x}>"
