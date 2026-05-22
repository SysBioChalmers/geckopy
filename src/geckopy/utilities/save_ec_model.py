"""Save an ecModel to disk.

The inverse of ``load_ec_model``. Writes the canonical geckopy
YAML format documented in ``docs/yaml_format.md`` for
``.yml`` / ``.yaml`` filenames, and SBML for ``.xml`` / ``.sbml``
(via ``geckopy.io.sbml.write_sbml_ec_model``).

The YAML writer uses ``cobra.io.dict.model_to_dict`` for the
cobra-shaped portion (metabolites, reactions, genes,
compartments) and **cobra's own YAML serializer** to write it, so
the file is byte-for-byte the cobrapy YAML format (``!!omap``
tagged). Four GECKO-specific top-level keys are added alongside:
``ec-rxns``, ``ec-enzymes``, ``gecko_light``, ``metaData``.
cobra-aware tools ignore those extra keys, so the file also loads
as a plain cobra model. Empty / NaN fields are omitted to keep the
file compact; the loader fills them back in.

MATLAB-COMPAT: MATLAB ``saveEcModel`` mutates the input model
(``ecModel.description = ['Enzyme-constrained model of '
ecModel.id]``). geckopy writes the description into the output
document only; the input is untouched.

MATLAB-COMPAT: MATLAB ``saveEcModel`` falls back to
``ModelAdapterManager.getDefault()`` when no adapter is supplied.
geckopy has no global default adapter — the caller passes one
explicitly, relies on ``model.adapter``, or uses an absolute
filename.

Ported from GECKO MATLAB: src/geckomat/utilities/saveEcModel.m.
"""
from __future__ import annotations

import math
from datetime import date
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import TYPE_CHECKING, Optional, Union

from cobra.io.dict import model_to_dict
from cobra.io.yaml import yaml as _cobra_yaml

from ..ec_model.ec_data import EcData
from ..ec_model.ec_model import EcModel

if TYPE_CHECKING:
    from ..adapter import ModelAdapter


_DEFAULT_FILENAME = "ecModel.yml"
_YAML_SUFFIXES = {".yml", ".yaml"}
_SBML_SUFFIXES = {".xml", ".sbml"}


def save_ec_model(
    model: EcModel,
    filename: Optional[Union[str, Path]] = None,
    adapter: Optional["ModelAdapter"] = None,
) -> Path:
    """Save an ecModel to disk in YAML or SBML format.

    Dispatches on the file extension:

    - ``.yml`` / ``.yaml`` -> geckopy canonical YAML (this module).
    - ``.xml`` / ``.sbml`` -> SBML (via
      ``geckopy.io.sbml.write_sbml_ec_model``).

    Parameters
    ----------
    model
        The ecModel to save. Must have populated ``model.ec``
        (both ``ec.rxns`` and ``ec.enzymes`` non-empty). Empty
        ecModels are rejected because they're almost always a bug
        (you probably forgot to run ``make_ec_model`` first).
    filename
        Destination path. If ``None``, defaults to ``ecModel.yml``.
        Relative paths are resolved against
        ``<adapter.params.path>/models/<filename>``. Absolute
        paths are used as-is.
    adapter
        Used to resolve a relative ``filename``. If not supplied
        as an argument, falls back to ``model.adapter``. Only
        required when ``filename`` is relative AND
        ``model.adapter`` is ``None``.

    Returns
    -------
    Path
        The absolute path of the file that was written.

    Raises
    ------
    ValueError
        If ``model.ec.rxns`` or ``model.ec.enzymes`` is empty, if
        the file extension isn't recognised, or if the filename is
        relative and no adapter is available.
    """
    if model.ec is None or len(model.ec.rxns) == 0 or len(model.ec.enzymes) == 0:
        raise ValueError(
            "model.ec is empty (zero rxns or zero enzymes); "
            "save_ec_model requires a populated ecModel. Run "
            "`make_ec_model` first to populate model.ec."
        )

    resolved_adapter = adapter if adapter is not None else getattr(
        model, "adapter", None,
    )
    path = _resolve_path(filename, resolved_adapter)

    if path.suffix.lower() in _SBML_SUFFIXES:
        from ..io.sbml import write_sbml_ec_model
        path.parent.mkdir(parents=True, exist_ok=True)
        write_sbml_ec_model(model, path)
        return path

    # cobra-shaped portion, exactly as cobrapy serialises it.
    doc = model_to_dict(model)
    # GECKO-specific extras as plain native types (cobra's YAML
    # serialiser renders every mapping as !!omap, so the extras match
    # the cobra portion's style). Kept out of model_to_dict so the
    # cobra part stays byte-identical to cobra's own output.
    doc["gecko_light"] = bool(model.ec.gecko_light)
    doc["metaData"] = _to_native(_build_metadata(model))
    doc["ec-rxns"] = _to_native(_build_ec_rxns_list(model.ec))
    doc["ec-enzymes"] = _to_native(_build_ec_enzymes_list(model.ec))

    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        _cobra_yaml.dump(doc, f)
    return path


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #

def _resolve_path(
    filename: Optional[Union[str, Path]],
    adapter: Optional["ModelAdapter"],
) -> Path:
    if filename is None:
        filename = _DEFAULT_FILENAME
    path = Path(filename)
    if (
        path.suffix.lower() not in _YAML_SUFFIXES
        and path.suffix.lower() not in _SBML_SUFFIXES
    ):
        raise ValueError(
            "ecModels are saved in YAML or SBML format only "
            f"(.yml/.yaml/.xml/.sbml). Got: {path.suffix!r}."
        )
    if not path.is_absolute():
        if adapter is None:
            raise ValueError(
                "filename is relative; supply an `adapter` (or set "
                "`model.adapter`) whose `params.path` will be used "
                "to resolve it under `<path>/models/<filename>`."
            )
        path = adapter.params.path / "models" / path
    return path


def _build_metadata(model: EcModel) -> dict:
    """Provenance fields written under the top-level `metaData` key.

    Currently populated with the save date, geckopy version, and a
    description string mirroring MATLAB's
    ``Enzyme-constrained model of <id>``.

    Future enrichment (Q2-b): author / email / organization /
    taxonomy could be sourced from `adapter.params` once those
    fields are added to `ModelParameters`. Tracked in
    `docs/future_improvements.md`.
    """
    try:
        geckopy_version = version("geckopy")
    except PackageNotFoundError:
        geckopy_version = "unknown"
    return {
        "date": date.today().isoformat(),
        "geckopy_version": geckopy_version,
        "description": f"Enzyme-constrained model of {model.id}",
    }


def _build_ec_rxns_list(ec: EcData) -> list:
    """Translate per-rxn ec fields + rxn_enz_mat rows into the
    canonical YAML list-of-mappings form.

    Empty `source` / `notes` / `eccodes` strings are omitted. `kcat`
    is always written: a real turnover number when set, otherwise
    ``0`` (the "no kcat assigned" sentinel, matching MATLAB GECKO).
    """
    coo = ec.rxn_enz_mat.tocoo()
    per_row_enzymes: list[dict[str, float]] = [
        {} for _ in range(ec.n_rxns)
    ]
    for i, j, v in zip(coo.row, coo.col, coo.data):
        per_row_enzymes[int(i)][ec.enzymes[int(j)]] = float(v)

    out = []
    for i in range(ec.n_rxns):
        entry: dict = {"id": ec.rxns[i], "kcat": float(ec.kcat[i])}
        if ec.source[i]:
            entry["source"] = ec.source[i]
        if ec.notes[i]:
            entry["notes"] = ec.notes[i]
        if ec.eccodes[i]:
            entry["eccodes"] = _eccodes_to_yaml(ec.eccodes[i])
        entry["enzymes"] = per_row_enzymes[i]
        out.append(entry)
    return out


def _build_ec_enzymes_list(ec: EcData) -> list:
    """Translate per-enzyme ec fields into the canonical YAML
    list-of-mappings form.

    NaN `mw` / `concs` and empty `sequence` are omitted; the loader
    fills them back in as NaN / empty string.
    """
    out = []
    for j in range(ec.n_enzymes):
        entry: dict = {"genes": ec.genes[j], "enzymes": ec.enzymes[j]}
        if not math.isnan(ec.mw[j]):
            entry["mw"] = float(ec.mw[j])
        if ec.sequence[j]:
            entry["sequence"] = ec.sequence[j]
        if not math.isnan(ec.concs[j]):
            entry["concs"] = float(ec.concs[j])
        out.append(entry)
    return out


def _to_native(value):
    """Recursively coerce ruamel.yaml scalar wrappers (ScalarInt,
    ScalarFloat, ...) and numpy scalars to native Python types.

    cobra-py's YAML loader stores ruamel scalar types on attributes
    like ``Metabolite.charge``; passing those through to the YAML
    safe-dumper raises RepresenterError. Coerce at the boundary so
    the on-disk file uses plain Python primitives.
    """
    if isinstance(value, dict):
        return {_to_native(k): _to_native(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_to_native(v) for v in value]
    if isinstance(value, bool):
        return bool(value)
    if isinstance(value, int):
        return int(value)
    if isinstance(value, float):
        return float(value)
    return value


def _eccodes_to_yaml(eccodes: str):
    """Convert an internal `;`-joined eccodes string back to the
    YAML form: a scalar string for one EC, a list for multiple."""
    parts = [p for p in eccodes.split(";") if p]
    if len(parts) <= 1:
        return eccodes
    return parts
