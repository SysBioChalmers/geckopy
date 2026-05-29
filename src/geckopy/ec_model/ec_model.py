"""ecModel: cobra.Model extended with enzyme-constraint data."""
from __future__ import annotations

from typing import TYPE_CHECKING, Optional, Union

import cobra

from .ec_data import EcData

if TYPE_CHECKING:
    from ..adapter import ModelAdapter


class EcModel(cobra.Model):
    """A cobra.Model with attached enzyme-constraint data.

    Stores the GECKO 3 `ec` substructure as `self.ec` and holds a
    reference to the ModelAdapter that was used to build it, so
    downstream functions only need the model object.

    Attributes
    ----------
    ec : EcData
        Enzyme-constraint data (kcats, enzymes, coupling matrix).
    adapter : ModelAdapter or None
        The adapter used when building the ecModel. May be None for
        models loaded from disk that do not carry adapter state.
    """

    def __init__(
        self,
        id_or_model: Optional[Union[str, cobra.Model]] = None,
        name: Optional[str] = None,
        *,
        adapter: Optional["ModelAdapter"] = None,
        gecko_light: bool = False,
    ):
        super().__init__(id_or_model, name)
        self.ec: EcData = EcData(gecko_light=gecko_light)
        self.adapter: Optional["ModelAdapter"] = adapter

    @classmethod
    def from_cobra(
        cls,
        model: cobra.Model,
        adapter: "ModelAdapter",
        *,
        gecko_light: bool = False,
    ) -> "EcModel":
        """Wrap an existing cobra.Model as an EcModel.

        Does not yet run makeEcModel; that is a separate pipeline step
        that mutates the wrapped model and populates `self.ec`.
        """
        ec_model = cls(model, adapter=adapter, gecko_light=gecko_light)
        return ec_model

    def validate_ec(self) -> None:
        """Validate internal consistency of the ec substructure."""
        self.ec.validate()
        unknown = set(self.ec.rxns) - {r.id for r in self.reactions}
        if unknown:
            preview = sorted(unknown)[:5]
            raise ValueError(
                f"{len(unknown)} reaction IDs in ec.rxns are not present "
                f"in the model (examples: {preview})"
            )