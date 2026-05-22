"""Load an ecModel from a YAML file.

Reads the cobrapy YAML format (``!!omap`` tagged, exactly what
``cobra.io.save_yaml_model`` writes) plus the GECKO-specific
top-level keys (``ec-rxns``, ``ec-enzymes``, ``gecko_light``,
``metaData``). The cobra-shaped portion is rebuilt by
``cobra.io.dict.model_from_dict``; this module reads the ec keys
on top.

The loader also dispatches SBML files (`.xml` / `.sbml`) to
``geckopy.io.sbml.read_sbml_ec_model``, so the same call works
for both formats.

Legacy MATLAB / RAVEN ecModels load too: ``_normalize_legacy_layout``
lifts ``id`` / ``name`` / ``version`` out of ``metaData`` and moves
per-metabolite ``smiles`` into ``annotation`` before handing the
document to cobra. (The legacy ``!!omap`` tags themselves parse
straight into mappings, so no special handling is needed for those.)

ecModels written before MATLAB GECKO switched ``usage_prot_*`` and
``prot_pool_exchange`` to the forward direction (positive flux)
used the opposite sign convention. The loader detects that case
on load and flips the affected reactions in place, so downstream
geckopy code always sees the current convention.

ec.kcat uses ``0`` as the "no kcat assigned" sentinel (matching
MATLAB GECKO).

MATLAB-COMPAT: the MATLAB ``loadEcModel`` validates
``endsWith(filename, {'yml','yaml'})`` and then has a dead
``elseif endsWith(filename,{'xml','sbml'})`` branch. geckopy
implements both paths properly; SBML is no longer dead code.

MATLAB-COMPAT: MATLAB ``loadEcModel`` falls back to
``ModelAdapterManager.getDefault()`` when no adapter is supplied.
geckopy has no global default adapter — the caller passes an
adapter explicitly, or uses an absolute filename.

Ported from GECKO MATLAB: src/geckomat/utilities/loadEcModel.m.
"""
from __future__ import annotations

import warnings
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
    """Load an ecModel from a YAML or SBML file.

    Dispatches on the file extension:

    - ``.yml`` / ``.yaml`` -> read as a geckopy canonical YAML
      ecModel (this module).
    - ``.xml`` / ``.sbml`` -> delegate to
      ``geckopy.io.sbml.read_sbml_ec_model``.

    Parameters
    ----------
    filename
        Path to the ecModel file. If ``None``, defaults to
        ``ecModel.yml``. If the path is relative, it is resolved
        against ``<adapter.params.path>/models/<filename>``. If
        absolute, it is used as-is and ``adapter`` may be
        ``None``.
    adapter
        A ``ModelAdapter`` used to resolve a relative ``filename``,
        and attached to the loaded model as ``model.adapter``.
        Required when ``filename`` is relative.

    Returns
    -------
    EcModel
        The loaded ecModel with ``model.ec`` populated. The
        ``gecko_light`` flag is read from the top-level
        ``gecko_light:`` key in the file (defaults to ``False``).

    Raises
    ------
    ValueError
        If the file extension isn't recognised, if the filename
        is relative and no adapter is given, or if the YAML lacks
        the GECKO-specific ``ec-rxns`` / ``ec-enzymes`` top-level
        keys (in which case it is not an ecModel).
    FileNotFoundError
        If the resolved file doesn't exist.
    """
    path = _resolve_path(filename, adapter)
    if path.suffix.lower() in _SBML_SUFFIXES:
        from ..io.sbml import read_sbml_ec_model
        return read_sbml_ec_model(path, adapter=adapter)
    data = _read_yaml(path)
    data = _normalize_legacy_layout(data)

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
    _flip_legacy_prot_direction(cobra_model)
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
    # cobra `!!omap` and plain mappings both parse to a dict. A very old
    # RAVEN file written as a bare `---` sequence of single-key maps (no
    # `!!omap`) parses to a list; merge it into one mapping.
    if isinstance(data, list):
        merged: dict = {}
        for item in data:
            if isinstance(item, dict):
                merged.update(item)
        data = merged
    if not isinstance(data, dict):
        raise ValueError(
            f"{path}: top-level YAML must be a mapping or a sequence of "
            f"single-key mappings, got {type(data).__name__}. "
            "See docs/yaml_format.md."
        )
    return data


def _normalize_legacy_layout(data: dict) -> dict:
    """Bring a legacy MATLAB / RAVEN ecModel YAML to the cobra layout.

    Two legacy quirks are normalised (both no-ops on a current file):

    - ``id`` / ``name`` / ``version`` nested under ``metaData`` are
      lifted to the top level (where cobra expects them).
    - a per-metabolite top-level ``smiles`` key is moved into
      ``annotation['smiles']`` (as a one-element list).
    """
    meta = data.get("metaData")
    if isinstance(meta, dict):
        for key in ("id", "name", "version"):
            if meta.get(key) and not data.get(key):
                data[key] = meta[key]

    mets = data.get("metabolites")
    if isinstance(mets, list):
        for met in mets:
            if isinstance(met, dict) and "smiles" in met:
                smiles = met.pop("smiles")
                annotation = met.get("annotation")
                if not isinstance(annotation, dict):
                    annotation = {}
                    met["annotation"] = annotation
                if "smiles" not in annotation and smiles:
                    annotation["smiles"] = (
                        smiles if isinstance(smiles, list) else [smiles]
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
    # 0 = "no kcat assigned"; real turnover numbers are always positive.
    kcat = np.array(
        [float(r.get("kcat", 0.0)) for r in ec_rxns_raw], dtype=float,
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


def _flip_legacy_prot_direction(model) -> None:
    """Detect and flip pre-forward-direction protein reactions in place.

    Older MATLAB GECKO ecModels defined ``usage_prot_*`` and
    ``prot_pool_exchange`` as "reverse" reactions: their flux was
    negative, and the stoichiometry signs were correspondingly
    swapped. The current convention (in both geckopy and recent
    MATLAB GECKO) treats them as ordinary forward reactions, with
    positive flux. When a loaded model still uses the older
    convention we flip the affected reactions in place so the rest
    of geckopy never has to handle two shapes.

    The signature we look for is any ``usage_prot_*`` or
    ``prot_pool_exchange`` reaction whose lower bound is negative.
    """
    flipped: list[str] = []
    for rxn in model.reactions:
        if not (
            rxn.id.startswith("usage_prot_")
            or rxn.id == "prot_pool_exchange"
        ):
            continue
        if rxn.lower_bound >= -1e-9:
            continue
        # Swap stoichiometry: subtract 2*current to land at -current.
        rxn.add_metabolites(
            {met: -2.0 * coef for met, coef in rxn.metabolites.items()},
            combine=True,
        )
        lb, ub = rxn.lower_bound, rxn.upper_bound
        rxn.lower_bound = -ub
        rxn.upper_bound = -lb
        flipped.append(rxn.id)
    if flipped:
        warnings.warn(
            f"ecModel uses the older reverse-direction convention "
            f"for {len(flipped)} protein usage/pool reaction(s); "
            "flipping to the current forward convention.",
            stacklevel=3,
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
    return ";".join(str(v) for v in value)
