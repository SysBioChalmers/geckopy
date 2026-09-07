"""Load experimental data for Bayesian kcat tuning.

Ported from GECKO MATLAB:
src/geckomat/kcat_sensitivity_analysis/Bayesian/loadBayesianData.m.

MATLAB's ``loadBayesianData`` loads ``bayesianFluxData.tsv`` and
``bayesianMaxGrowth.tsv`` through the *same* ``loadFluxData`` parser used
for regular flux data (verified against ``develop4``'s current source,
not assumed) -- both files share the identical
``Condition Ptot grRate <met> (<rxn>) ... bayesianRMSEweight source``
column layout. geckopy already has that parser as
:func:`~geckopy.databases.flux_data.load_flux_data`, so both files reuse
it unchanged; only ``bayesianZeroExch.tsv`` (a single ``Rxns`` id column)
needs a loader of its own.

Unlike MATLAB -- which treats a missing/unreadable file as `[]` via a
try/catch inside ``loadFluxData`` -- all three files are treated
uniformly here as optional: a missing file yields ``None`` (or an empty
list for ``zero_flux``), never an exception. This also sidesteps a bug
flagged in the MATLAB review (``REVIEW.md`` #3): there, only the
``fluxData`` field was null-guarded before ``.biomass`` was attached,
while ``maxGrate`` wasn't, so a missing ``bayesianMaxGrowth.tsv`` there
silently produced a struct with just a ``.biomass`` field instead of a
clean "absent" signal.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Optional

from ...databases.flux_data import FluxData, load_flux_data

if TYPE_CHECKING:
    from ...adapter import ModelAdapter

_FLUX_DATA_FILENAME = "bayesianFluxData.tsv"
_MAX_GROWTH_FILENAME = "bayesianMaxGrowth.tsv"
_ZERO_EXCH_FILENAME = "bayesianZeroExch.tsv"


@dataclass
class BayesianData:
    """Experimental data used by Bayesian kcat tuning.

    Attributes
    ----------
    flux_data
        Per-condition growth/exchange-flux measurements, from
        ``bayesianFluxData.tsv``. ``None`` if that file isn't present.
    max_grate
        Per-condition maximum-growth-rate measurements (one active
        carbon source per condition, at ``-1000`` i.e. unconstrained),
        from ``bayesianMaxGrowth.tsv``. ``None`` if that file isn't
        present.
    zero_flux
        Reaction IDs assumed to carry zero flux, from
        ``bayesianZeroExch.tsv``. Empty if that file isn't present.
    """

    flux_data: Optional[FluxData] = None
    max_grate: Optional[FluxData] = None
    zero_flux: list[str] = field(default_factory=list)


def load_bayesian_data(adapter: "ModelAdapter") -> BayesianData:
    """Load the three Bayesian-kcat-tuning experimental data files.

    Each file is optional and resolved as
    ``<adapter.params.path>/data/<filename>``; a missing file yields
    ``None`` (``flux_data``/``max_grate``) or an empty list
    (``zero_flux``) rather than raising.

    Parameters
    ----------
    adapter
        A ``ModelAdapter`` whose ``params.path`` locates the ``data/``
        folder.

    Returns
    -------
    BayesianData
    """
    data_dir = Path(adapter.params.path) / "data"
    return BayesianData(
        flux_data=_load_optional_flux_data(data_dir / _FLUX_DATA_FILENAME),
        max_grate=_load_optional_flux_data(data_dir / _MAX_GROWTH_FILENAME),
        zero_flux=_load_zero_exch(data_dir / _ZERO_EXCH_FILENAME),
    )


def _load_optional_flux_data(path: Path) -> Optional[FluxData]:
    if not path.is_file():
        return None
    return load_flux_data(path)


def _load_zero_exch(path: Path) -> list[str]:
    """Parse ``bayesianZeroExch.tsv``: a single ``Rxns`` header, then
    one reaction ID per line."""
    if not path.is_file():
        return []
    lines = path.read_text(encoding="utf-8").splitlines()
    rxn_ids = [line.strip() for line in lines[1:] if line.strip()]
    return rxn_ids
