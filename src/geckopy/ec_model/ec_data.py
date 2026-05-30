"""Re-export of :class:`raven_python.io.EcData`.

``EcData`` is the typed enzyme-constraint substructure attached to every
ecModel. The dataclass and its YAML schema (the ``ec-rxns`` / ``ec-enzymes``
/ ``gecko_light`` top-level sections, the ``EcData.empty`` preallocation
factory, the ``EcData.validate`` shape-integrity check) live in
raven-python, in line with MATLAB RAVEN where ``model.ec`` is a RAVEN-owned
struct and GECKO is just one of its consumers. geckopy re-exports the
class so existing ``from geckopy import EcData`` imports keep working.

Algorithmic operations on an ``EcData`` (kcat sourcing, isozyme
expansion, sensitivity tuning, ...) stay in geckopy — that's where the
GECKO pipeline lives. The split mirrors the RAVEN/GECKO MATLAB layering.
"""
from raven_python.io import EcData

__all__ = ["EcData"]
