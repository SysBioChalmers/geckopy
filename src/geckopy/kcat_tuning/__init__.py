"""Kcat sensitivity analysis and tuning against experimental data.

Mirrors GECKO MATLAB's `src/kcat_tuning/`
directory.
"""
from .find_max_value import find_max_value
from .sensitivity_tuning import TunedKcatsResult, sensitivity_tuning
from .sigma_fitter import SigmaFitterResult, fit_sigma, sigma_fitter
from .truncate_values import truncate_values

try:
    from . import evolutionary_tuning
except ImportError as exc:  # pragma: no cover - exercised
    # only when the optional `cma` dependency is absent.
    # Bound to a module-level name: Python clears an `except ... as`
    # target when the block exits, so the class below cannot close over
    # `exc` itself.
    _evolutionary_tuning_import_error = exc

    class _MissingEvolutionaryTuning:
        """Stand-in for the `evolutionary_tuning` subpackage when its
        `cma` dependency isn't installed. Raises only on first attribute
        access, not on import, so plain `geckopy.kcat_tuning`
        stays usable without the optional extra."""

        def __getattr__(self, name: str):
            raise ImportError(
                "Evolutionary (CMA-ES) kcat tuning requires the optional "
                "'evolutionary-tuning' extra: pip install "
                "geckopy[evolutionary-tuning]"
            ) from _evolutionary_tuning_import_error

    evolutionary_tuning = _MissingEvolutionaryTuning()  # type: ignore[assignment]

__all__ = [
    "SigmaFitterResult",
    "TunedKcatsResult",
    "evolutionary_tuning",
    "find_max_value",
    "fit_sigma",
    "sensitivity_tuning",
    "sigma_fitter",  # deprecated; alias of fit_sigma
    "truncate_values",
]
