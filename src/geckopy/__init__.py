"""geckopy – enzyme-constrained metabolic modeling in Python.

Provides :class:`~geckopy.core.ECModel`, :class:`~geckopy.core.ECData`, and
:class:`~geckopy.core.Enzyme` as the core data structures for building and
analysing enzyme-constrained genome-scale models (ecGEMs), following the
conventions of the GECKO Toolbox (Sánchez et al. 2017; Domenzain et al. 2022).
"""

from geckopy.core.ec_model import ECData, ECModel
from geckopy.core.enzyme import Enzyme


__all__ = [
    "ECData",
    "ECModel",
    "Enzyme",
]
