"""Sensitivity analysis on enzyme kcat values.

Mirrors GECKO MATLAB's `src/geckomat/kcat_sensitivity_analysis/`
directory.
"""
from .find_max_value import find_max_value
from .sensitivity_tuning import TunedKcatsResult, sensitivity_tuning
from .sigma_fitter import SigmaFitterResult, fit_sigma, sigma_fitter
from .truncate_values import truncate_values

try:
    from . import bayesian
except ImportError as exc:  # pragma: no cover - exercised
    # only when the optional `pyabc` dependency is absent.
    # Bound to a module-level name: Python clears an `except ... as`
    # target when the block exits, so the class below cannot close over
    # `exc` itself.
    _bayesian_import_error = exc

    class _MissingBayesian:
        """Stand-in for the `bayesian` subpackage when its `pyabc`
        dependency isn't installed. Raises only on first attribute
        access, not on import, so plain `geckopy.kcat_sensitivity_analysis`
        stays usable without the optional extra."""

        def __getattr__(self, name: str):
            raise ImportError(
                "Bayesian kcat tuning requires the optional 'bayesian' "
                "extra: pip install geckopy[bayesian]"
            ) from _bayesian_import_error

    bayesian = _MissingBayesian()  # type: ignore[assignment]

__all__ = [
    "SigmaFitterResult",
    "TunedKcatsResult",
    "bayesian",
    "find_max_value",
    "fit_sigma",
    "sensitivity_tuning",
    "sigma_fitter",  # deprecated; alias of fit_sigma
    "truncate_values",
]
