"""Kcat sensitivity analysis and tuning against experimental data.

Mirrors GECKO MATLAB's `src/kcat_tuning/`
directory.
"""
from .corrections import (
    Correction,
    annotate_from_model,
    corrections,
    corrections_tsv,
)
from .find_max_value import find_max_value
from .sensitivity_tuning import TunedKcatsResult, sensitivity_tuning
from .sigma_fitter import SigmaFitterResult, fit_sigma, sigma_fitter
from .truncate_values import truncate_values

try:
    from . import evotune
except ImportError as exc:  # pragma: no cover - exercised
    # only when the optional `cma` dependency is absent.
    # Bound to a module-level name: Python clears an `except ... as`
    # target when the block exits, so the class below cannot close over
    # `exc` itself.
    _evotune_import_error = exc

    class _MissingEvotune:
        """Stand-in for the `evotune` subpackage when its
        `cma` dependency isn't installed. Raises only on first attribute
        access, not on import, so plain `geckopy.kcat_tuning`
        stays usable without the optional extra."""

        def __getattr__(self, name: str):
            raise ImportError(
                "Evotune: CMA-ES kcat tuning requires the optional "
                "'evotune' extra: pip install "
                "geckopy[evotune]"
            ) from _evotune_import_error

    evotune = _MissingEvotune()  # type: ignore[assignment]

__all__ = [
    "Correction",
    "SigmaFitterResult",
    "TunedKcatsResult",
    "annotate_from_model",
    "corrections",
    "corrections_tsv",
    "evotune",
    "find_max_value",
    "fit_sigma",
    "sensitivity_tuning",
    "sigma_fitter",  # deprecated; alias of fit_sigma
    "truncate_values",
]
