"""Load an ecModel from a YAML file.

Ported from GECKO MATLAB:
src/geckomat/utilities/loadEcModel.m.

The on-disk format is documented in `docs/yaml_format.md`. It is
designed to be a strict superset of cobra-py's canonical YAML schema:
the cobra-shaped portion (`metabolites`, `reactions`, `genes`,
`compartments`, ...) is loaded directly via
`cobra.io.dict.model_from_dict`, while the GECKO-specific top-level
keys (`ec-rxns`, `ec-enzymes`, `gecko_light`, `metaData`) are handled
by this module.

This loader does NOT accept the legacy MATLAB / RAVEN format
(`!!omap`-tagged, outer sequence wrapper, `metaData`-nested `id`).
That format will be handled by an external one-off conversion script
once MATLAB-side `writeYAMLmodel` is updated; see
`docs/future_improvements.md` and `docs/yaml_format.md` for the
canonical schema and the migration table.

MATLAB-COMPAT: `loadEcModel.m` validates `endsWith(filename,
{'yml','yaml'})` and errors otherwise, then has an `elseif
endsWith(filename, {'xml','sbml'})` branch which is unreachable.
geckopy keeps the strict YAML check (mirror of MATLAB's documented
intent) and drops the dead SBML branch.

MATLAB-COMPAT: MATLAB `loadEcModel` falls back to
`ModelAdapterManager.getDefault()` when no adapter is supplied.
geckopy has no global default adapter; the caller must either
pass an adapter explicitly OR pass an absolute path that does not
need adapter-relative resolution.
"""
from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Optional, Union

import numpy as np
from cobra.io.dict import model_from_dict
from ruamel.yaml import YAML
from scipy import sparse

from ..ec_model.ec_data import EcData
from ..ec_model.ec_model import EcModel

if TYPE_CHECKING:
    from ..adapter import ModelAdapter


_DEFAULT_FILENAME = "ecModel.yml"
_YAML_SUFFIXES = {".yml", ".yaml"}
_SBML_SUFFIXES = {".xml", ".sbml"}


def load_ec_model(
    filename: Optional[Union[str, Path]] = None,
    adapter: Optional["ModelAdapter"] = None,
) -> EcModel:
    """Load an ecModel YAML file from disk.

    Parameters
    ----------
    filename
        Path to the ecModel YAML. If `None`, defaults to
        `ecModel.yml`. If relative, resolved against
        `<adapter.params.path>/models/<filename>`. If absolute,
        used as-is and `adapter` may be `None`.
    adapter
        ModelAdapter used to resolve a relative `filename`. The
        loaded `EcModel` carries this adapter as `model.adapter`.
        Required when `filename` is relative.

    Returns
    -------
    EcModel
        The loaded ecModel with `model.ec` populated. The
        `gecko_light` flag is read from the file (top-level
        `gecko_light:` key) and propagated to the EcData.

    Raises
    ------
    ValueError
        If the file extension is not `.yml` or `.yaml`, if the
        filename is relative and no adapter is given, or if the
        YAML lacks the GECKO-specific `ec-rxns` / `ec-enzymes`
        keys (i.e. it is not an ecModel).
    FileNotFoundError
        If the resolved file does not exist.
    """
    path = _resolve_path(filename, adapter)
    if path.suffix.lower() in _SBML_SUFFIXES:
        from ..io.sbml import read_sbml_ec_model
        return read_sbml_ec_model(path, adapter=adapter)
    data = _read_yaml(path)

    ec_rxns_raw = data.pop("ec-rxns", None)
    ec_enzymes_raw = data.pop("ec-enzymes", None)
    if ec_rxns_raw is None or ec_enzymes_raw is None:
        raise ValueError(
            f"{path}: YAML lacks `ec-rxns` and/or `ec-enzymes` "
            "top-level keys; this is not a geckopy ecModel YAML"
        )

    gecko_light = bool(data.pop("gecko_light", False))
    data.pop("metaData", None)  # preserved on disk; ignored on load

    cobra_model = model_from_dict(data)
    ec_data = _build_ec_data(
        ec_rxns_raw, ec_enzymes_raw, gecko_light=gecko_light,
    )

    ec_model = EcModel.from_cobra(
        cobra_model, adapter=adapter, gecko_light=gecko_light,
    )
    ec_model.ec = ec_data
    return ec_model


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
            "ecModels are distributed in YAML or SBML format "
            f"(.yml/.yaml/.xml/.sbml). Got: {path.suffix!r}."
        )
    if not path.is_absolute():
        if adapter is None:
            raise ValueError(
                "filename is relative; supply an `adapter` whose "
                "`params.path` will be used to resolve it under "
                "`<path>/models/<filename>`."
            )
        path = adapter.params.path / "models" / path
    if not path.is_file():
        raise FileNotFoundError(f"ecModel file not found: {path}")
    return path


def _read_yaml(path: Path) -> dict:
    yaml = YAML(typ="safe")
    with open(path, "r", encoding="utf-8") as f:
        data = yaml.load(f)
    if not isinstance(data, dict):
        raise ValueError(
            f"{path}: top-level YAML must be a mapping, got "
            f"{type(data).__name__}. Legacy MATLAB / RAVEN format "
            "(outer sequence wrapper) is not supported; see "
            "docs/yaml_format.md."
        )
    return data


def _build_ec_data(
    ec_rxns_raw: list,
    ec_enzymes_raw: list,
    *,
    gecko_light: bool,
) -> EcData:
    """Construct `EcData` from the parsed `ec-rxns` and `ec-enzymes`
    YAML lists."""
    n_e = len(ec_enzymes_raw)
    genes = [str(e["genes"]) for e in ec_enzymes_raw]
    enzymes = [str(e["enzymes"]) for e in ec_enzymes_raw]
    mw = np.array(
        [float(e.get("mw", np.nan)) for e in ec_enzymes_raw], dtype=float,
    )
    sequence = [str(e.get("sequence", "")) for e in ec_enzymes_raw]
    concs = np.array(
        [float(e.get("concs", np.nan)) for e in ec_enzymes_raw], dtype=float,
    )

    enz_index = {eid: i for i, eid in enumerate(enzymes)}

    n_r = len(ec_rxns_raw)
    rxns = [str(r["id"]) for r in ec_rxns_raw]
    kcat = np.array(
        [float(r.get("kcat", np.nan)) for r in ec_rxns_raw], dtype=float,
    )
    source = [str(r.get("source", "")) for r in ec_rxns_raw]
    notes = [str(r.get("notes", "")) for r in ec_rxns_raw]
    eccodes = [_canonicalize_eccodes(r.get("eccodes", "")) for r in ec_rxns_raw]

    mat = sparse.lil_matrix((n_r, n_e), dtype=float)
    for i, r in enumerate(ec_rxns_raw):
        for enz_id, stoich in (r.get("enzymes") or {}).items():
            j = enz_index.get(str(enz_id))
            if j is None:
                raise ValueError(
                    f"ec-rxns[{i}] (id={r.get('id')!r}) references "
                    f"enzyme {enz_id!r} that is not present in "
                    "ec-enzymes."
                )
            mat[i, j] = float(stoich)

    return EcData(
        gecko_light=gecko_light,
        rxns=rxns,
        kcat=kcat,
        source=source,
        notes=notes,
        eccodes=eccodes,
        genes=genes,
        enzymes=enzymes,
        mw=mw,
        sequence=sequence,
        concs=concs,
        rxn_enz_mat=mat.tocsr(),
    )


def _canonicalize_eccodes(value) -> str:
    """Coerce an EC-codes field to a single `;`-joined string.

    The schema accepts either a scalar string (`"1.1.1.1"`) or a
    list of strings (`["1.1.1.1", "1.1.99.40"]`); both round-trip
    to the same internal representation.
    """
    if value is None:
        return ""
    if isinstance(value, str):
        return value
    if isinstance(value, (list, tuple)):
        return ";".join(str(v) for v in value)
    return str(value)
