"""Define ECData and ECModel for enzyme-constrained metabolic models.

An enzyme-constrained model (ECModel) extends :class:`cobra.Model` with
enzyme capacity constraints derived from catalytic rate constants (kcat) and
enzyme molecular weights. This implementation follows the conventions of the
GECKO Toolbox (Sánchez et al., 2017; Domenzain et al., 2022).

Two model flavors are supported:

* **Full model** – each enzyme gets its own pseudometabolite (``prot_<ID>``)
  and usage pseudoreaction (``usage_prot_<ID>``), all drawing from a shared
  ``prot_pool`` metabolite.
* **Light model** (``geckoLight=True``) – only the single most efficient
  isozyme per reaction contributes to the shared ``prot_pool`` directly;
  no individual enzyme pseudoreactions are added.

References
----------
Sánchez et al. (2017) Mol Syst Biol 13:935.
Domenzain et al. (2022) Nat Commun 13:3136.
"""

import logging
from copy import deepcopy
from typing import Dict, Iterable, List, Optional, Union

import numpy as np
from cobra import Model
from cobra.core.dictlist import DictList
from cobra.core.metabolite import Metabolite
from cobra.core.reaction import Reaction

from geckopy.core.enzyme import Enzyme


logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

#: Prefix applied to enzyme pseudometabolite IDs (e.g. ``prot_P00549``).
ENZYME_METABOLITE_PREFIX = "prot_"
#: ID of the shared protein pool pseudometabolite.
POOL_METABOLITE_ID = "prot_pool"
#: Prefix applied to enzyme usage pseudoreaction IDs (e.g. ``usage_prot_P00549``).
ENZYME_REACTION_PREFIX = "usage_prot_"
#: ID of the protein pool exchange pseudoreaction.
POOL_REACTION_ID = "prot_pool_exchange"


# ---------------------------------------------------------------------------
# ECData
# ---------------------------------------------------------------------------


class ECData:
    """Container for enzyme-constraint data, mirroring GECKO's ``ec`` struct.

    All parallel arrays share the same length *n* where *n* is the number of
    reaction–enzyme associations (one enzyme may appear in several reactions
    and one reaction may be catalyzed by several isozymes).

    Parameters
    ----------
    geckoLight : bool, optional
        When True the model uses the light formulation (single pool constraint
        per reaction, no individual enzyme pseudoreactions). Default False.

    Attributes
    ----------
    geckoLight : bool
        Whether the light model formulation is used.
    rxns : list of str
        Reaction IDs for each association entry.
    kcat : numpy.ndarray of float
        Catalytic rate constants in s⁻¹.
    source : list of str
        Provenance of each kcat value (e.g. ``"BRENDA"``, ``"DLKcat"``).
    notes : list of str
        Free-text annotations per entry.
    eccodes : list of str
        EC numbers (e.g. ``"1.2.3.4"``) per entry.
    genes : list of str
        Gene identifiers corresponding to each enzyme entry.
    enzymes : list of str
        Enzyme/UniProt accessions for each entry.
    mw : numpy.ndarray of float
        Molecular weights in g/mol for each entry.
    sequence : list of str
        Amino acid sequences for each entry.
    concs : numpy.ndarray of float
        Measured enzyme concentrations in mmol/gDCW (NaN when unavailable).
    rxnEnzMat : numpy.ndarray or None
        2-D array of shape *(n_reactions, n_enzymes)* containing subunit copy
        numbers. Rows correspond to unique reactions, columns to unique
        enzymes. Populated by complex-stoichiometry data; None when not set.
    """

    def __init__(self, geckoLight: bool = False) -> None:
        """Initialize an empty ECData container."""
        self.geckoLight: bool = geckoLight

        self.rxns: List[str] = []
        self.kcat: np.ndarray = np.empty(0, dtype=float)
        self.source: List[str] = []
        self.notes: List[str] = []
        self.eccodes: List[str] = []
        self.genes: List[str] = []
        self.enzymes: List[str] = []
        self.mw: np.ndarray = np.empty(0, dtype=float)
        self.sequence: List[str] = []
        self.concs: np.ndarray = np.full(0, np.nan, dtype=float)

        self.rxnEnzMat: Optional[np.ndarray] = None

    @property
    def n(self) -> int:
        """Number of reaction–enzyme association entries."""
        return len(self.rxns)

    def to_dict(self) -> dict:
        """Return a plain-Python-dict representation of this ECData.

        Returns
        -------
        dict
            Keys match attribute names; numpy arrays are converted to lists.
        """
        return {
            "geckoLight": self.geckoLight,
            "rxns": list(self.rxns),
            "kcat": self.kcat.tolist(),
            "source": list(self.source),
            "notes": list(self.notes),
            "eccodes": list(self.eccodes),
            "genes": list(self.genes),
            "enzymes": list(self.enzymes),
            "mw": self.mw.tolist(),
            "sequence": list(self.sequence),
            "concs": self.concs.tolist(),
            "rxnEnzMat": (
                self.rxnEnzMat.tolist() if self.rxnEnzMat is not None else None
            ),
        }

    @classmethod
    def from_dict(cls, data: dict) -> "ECData":
        """Reconstruct an ECData from a plain-Python dict.

        Parameters
        ----------
        data : dict
            As produced by :meth:`to_dict`.

        Returns
        -------
        ECData
        """
        ec = cls(geckoLight=data.get("geckoLight", False))
        ec.rxns = list(data.get("rxns", []))
        ec.kcat = np.array(data.get("kcat", []), dtype=float)
        ec.source = list(data.get("source", []))
        ec.notes = list(data.get("notes", []))
        ec.eccodes = list(data.get("eccodes", []))
        ec.genes = list(data.get("genes", []))
        ec.enzymes = list(data.get("enzymes", []))
        ec.mw = np.array(data.get("mw", []), dtype=float)
        ec.sequence = list(data.get("sequence", []))
        concs_raw = data.get("concs", [])
        ec.concs = np.array(concs_raw, dtype=float) if concs_raw else np.full(0, np.nan)
        rxn_enz = data.get("rxnEnzMat")
        ec.rxnEnzMat = np.array(rxn_enz, dtype=float) if rxn_enz is not None else None
        return ec

    def copy(self) -> "ECData":
        """Return a deep copy of this ECData."""
        return deepcopy(self)

    def __repr__(self) -> str:
        return (
            f"<ECData geckoLight={self.geckoLight}, "
            f"n={self.n} reaction-enzyme associations>"
        )


# ---------------------------------------------------------------------------
# ECModel
# ---------------------------------------------------------------------------


class ECModel(Model):
    """Enzyme-constrained genome-scale metabolic model.

    Extends :class:`cobra.Model` with:

    * an :attr:`enzymes` DictList of :class:`~geckopy.core.Enzyme` objects,
    * an :attr:`ec` attribute (:class:`ECData`) that stores all enzyme-constraint
      metadata (kcat values, molecular weights, gene associations, etc.),
    * convenience properties to access the enzyme pseudometabolites, usage
      pseudoreactions, and the protein pool metabolite/reaction present in the
      underlying stoichiometric matrix.

    Parameters
    ----------
    id_or_model : str or cobra.Model or ECModel, optional
        * **string** – used as the model ID; an empty model is created.
        * :class:`cobra.Model` – wraps an existing model; :attr:`ec` will be empty.
        * :class:`ECModel` – full copy (including EC data).
        * None – empty model with id ``None``.
    name : str, optional
        Human-readable model description.
    geckoLight : bool, optional
        Whether to use the light model formulation (default False).

    Attributes
    ----------
    enzymes : DictList of Enzyme
        All enzymes registered in this model, keyed by UniProt accession.
    ec : ECData
        Enzyme-constraint metadata.
    """

    ENZYME_METABOLITE_PREFIX: str = ENZYME_METABOLITE_PREFIX
    POOL_METABOLITE_ID: str = POOL_METABOLITE_ID
    ENZYME_REACTION_PREFIX: str = ENZYME_REACTION_PREFIX
    POOL_REACTION_ID: str = POOL_REACTION_ID

    def __init__(
        self,
        id_or_model: Union[str, Model, "ECModel", None] = None,
        name: Optional[str] = None,
        geckoLight: bool = False,
    ) -> None:
        """Initialize an ECModel."""
        if isinstance(id_or_model, ECModel):
            super().__init__(id_or_model=id_or_model, name=name)
            self.enzymes: DictList[Enzyme] = deepcopy(id_or_model.enzymes)
            self.ec: ECData = id_or_model.ec.copy()
            for enzyme in self.enzymes:
                enzyme._model = self
        elif isinstance(id_or_model, Model):
            super().__init__(id_or_model=id_or_model, name=name)
            self.enzymes = DictList()
            self.ec = ECData(geckoLight=geckoLight)
        else:
            super().__init__(id_or_model=id_or_model, name=name)
            self.enzymes = DictList()
            self.ec = ECData(geckoLight=geckoLight)

    # ------------------------------------------------------------------
    # Enzyme management
    # ------------------------------------------------------------------

    def add_enzymes(self, enzymes: Iterable[Enzyme]) -> None:
        """Add enzyme objects to the model.

        Parameters
        ----------
        enzymes : iterable of Enzyme
            Enzyme objects to add. Enzymes already present (same id) are
            skipped with a warning.
        """
        for enzyme in list(enzymes):
            if not isinstance(enzyme, Enzyme):
                raise TypeError(f"Expected Enzyme, got {type(enzyme)}")
            if self.enzymes.has_id(enzyme.id):
                logger.warning(
                    "Enzyme %s is already in the model, skipping.", enzyme.id
                )
                continue
            enzyme._model = self
            self.enzymes.append(enzyme)

    def remove_enzymes(
        self, enzymes: Iterable[Enzyme], prune: bool = False
    ) -> None:
        """Remove enzyme objects from the model.

        Parameters
        ----------
        enzymes : iterable of Enzyme
            Enzyme objects to remove.
        prune : bool, optional
            If True, also remove the associated pseudometabolite and usage
            pseudoreaction from the model (if they exist). Default False.
        """
        for enzyme in list(enzymes):
            if not self.enzymes.has_id(enzyme.id):
                logger.warning(
                    "Enzyme %s is not in the model, skipping.", enzyme.id
                )
                continue
            if prune:
                self._remove_enzyme_reactions(enzyme)
            enzyme._model = None
            self.enzymes.remove(enzyme)

    def _remove_enzyme_reactions(self, enzyme: Enzyme) -> None:
        """Remove the pseudometabolite and usage reaction for an enzyme."""
        usage_rxn_id = enzyme.usage_reaction_id
        if usage_rxn_id and self.reactions.has_id(usage_rxn_id):
            self.remove_reactions([self.reactions.get_by_id(usage_rxn_id)])

        prot_met_id = enzyme.pseudometabolite_id
        if prot_met_id and self.metabolites.has_id(prot_met_id):
            self.remove_metabolites([self.metabolites.get_by_id(prot_met_id)])

    # ------------------------------------------------------------------
    # Convenience properties
    # ------------------------------------------------------------------

    @property
    def enzyme_metabolites(self) -> List[Metabolite]:
        """Return all individual enzyme pseudometabolites (``prot_*``).

        Excludes the shared protein pool metabolite; use :attr:`pool_metabolite`
        to access that.

        Returns
        -------
        list of Metabolite
        """
        prefix = self.ENZYME_METABOLITE_PREFIX
        pool_id = self.POOL_METABOLITE_ID
        return [
            m for m in self.metabolites
            if m.id.startswith(prefix) and m.id != pool_id
        ]

    @property
    def pool_metabolite(self) -> Optional[Metabolite]:
        """Return the protein pool pseudometabolite, or None if absent.

        Returns
        -------
        Metabolite or None
        """
        pool_id = self.POOL_METABOLITE_ID
        return (
            self.metabolites.get_by_id(pool_id)
            if self.metabolites.has_id(pool_id)
            else None
        )

    @property
    def enzyme_reactions(self) -> List[Reaction]:
        """Return all enzyme usage pseudoreactions (``usage_prot_*``).

        Returns
        -------
        list of Reaction
        """
        prefix = self.ENZYME_REACTION_PREFIX
        return [r for r in self.reactions if r.id.startswith(prefix)]

    @property
    def pool_reaction(self) -> Optional[Reaction]:
        """Return the protein pool exchange pseudoreaction, or None if absent.

        Returns
        -------
        Reaction or None
        """
        pool_rxn_id = self.POOL_REACTION_ID
        return (
            self.reactions.get_by_id(pool_rxn_id)
            if self.reactions.has_id(pool_rxn_id)
            else None
        )

    @property
    def constrain_pool(self) -> Optional[float]:
        """Return the lower bound of the protein pool exchange reaction.

        Returns None when the pool reaction is not present.

        Returns
        -------
        float or None
        """
        rxn = self.pool_reaction
        return rxn.lower_bound if rxn is not None else None

    @constrain_pool.setter
    def constrain_pool(self, value: float) -> None:
        """Set the lower bound of the protein pool exchange reaction.

        Parameters
        ----------
        value : float
            Lower bound (negative, since the exchange reaction consumes the pool).
        """
        rxn = self.pool_reaction
        if rxn is None:
            raise ValueError(
                f"Protein pool reaction '{self.POOL_REACTION_ID}' is not "
                "present in the model. Add it before setting the pool constraint."
            )
        rxn.lower_bound = value

    # ------------------------------------------------------------------
    # Serialization
    # ------------------------------------------------------------------

    def __setstate__(self, state: dict) -> None:
        """Restore model state and re-wire enzyme back-references."""
        super().__setstate__(state)
        if not hasattr(self, "enzymes"):
            self.enzymes = DictList()
        if not hasattr(self, "ec"):
            self.ec = ECData()
        for enzyme in self.enzymes:
            enzyme._model = self

    def copy(self) -> "ECModel":
        """Return a copy of this ECModel."""
        return ECModel(id_or_model=self)

    def __repr__(self) -> str:
        return (
            f"<ECModel {self.id} at {id(self):#x} "
            f"({len(self.reactions)} reactions, "
            f"{len(self.metabolites)} metabolites, "
            f"{len(self.enzymes)} enzymes)>"
        )
