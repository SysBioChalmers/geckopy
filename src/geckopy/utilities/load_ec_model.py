"""Load an ecModel from a YAML file.

Thin wrapper around :func:`raven_python.io.read_yaml_model`. raven-python
owns the YAML schema for the ec sections (see
``raven_python.io.ec_data``): it reads the cobra-shaped portion plus the
``ec-rxns`` / ``ec-enzymes`` / ``gecko_light`` top-level sections,
populates ``model.ec`` as a typed :class:`~raven_python.io.EcData`, and
normalises three legacy MATLAB GECKO quirks (top-level ``smiles`` on
metabolites, reverse-direction ``usage_prot_*`` / ``prot_pool_exchange``,
bare-``-`` document root).

geckopy keeps three application-level concerns:

1. **File-extension validation** — only ``.yml`` / ``.yaml`` is accepted.
   SBML ecModel I/O has been removed; raise a helpful ValueError if the
   caller passes ``.xml`` / ``.sbml``.
2. **Adapter-aware path resolution** — relative filenames resolve under
   ``<adapter.params.path>/models/<filename>``.
3. **Not-an-ecModel guard + adapter attachment** — fast-fail with a
   clear message if the YAML lacks the ec sections, and attach the
   adapter to ``model.adapter`` so downstream geckopy code finds it.

MATLAB-COMPAT: MATLAB ``loadEcModel`` falls back to
``ModelAdapterManager.getDefault()`` when no adapter is supplied.
geckopy has no global default adapter — pass an adapter explicitly, or
use an absolute filename.

Ported from GECKO MATLAB: src/geckomat/utilities/loadEcModel.m.
"""
from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Optional, Union

from raven_python.io import EcData, read_yaml_model

from ..ec_model.ec_model import EcModel

if TYPE_CHECKING:
    from ..adapter import ModelAdapter


_DEFAULT_FILENAME = "ecModel.yml"
_YAML_SUFFIXES = {".yml", ".yaml"}


def load_ec_model(
    filename: Optional[Union[str, Path]] = None,
    adapter: Optional["ModelAdapter"] = None,
) -> EcModel:
    """Load an ecModel from a YAML file.

    Parameters
    ----------
    filename
        Path to the ecModel file. If ``None``, defaults to ``ecModel.yml``.
        Relative paths resolve against
        ``<adapter.params.path>/models/<filename>``; absolute paths are
        used as-is and ``adapter`` may be ``None``.
    adapter
        A ``ModelAdapter`` used to resolve a relative ``filename`` and
        attached to the loaded model as ``model.adapter``. Required when
        ``filename`` is relative.

    Returns
    -------
    EcModel
        The loaded ecModel with a populated ``model.ec``. The
        ``gecko_light`` flag, ec arrays, and coupling matrix are all
        read from the file.

    Raises
    ------
    ValueError
        If the extension isn't YAML, the filename is relative without an
        adapter, or the YAML lacks the GECKO ``ec-rxns`` / ``ec-enzymes``
        sections (in which case it isn't an ecModel).
    FileNotFoundError
        If the resolved file doesn't exist.
    """
    path = _resolve_path(filename, adapter)

    cobra_model = read_yaml_model(path)
    ec = getattr(cobra_model, "ec", None)
    if not isinstance(ec, EcData) or ec.n_rxns == 0 or ec.n_enzymes == 0:
        raise ValueError(
            f"{path}: YAML lacks the GECKO `ec-rxns` and/or `ec-enzymes` "
            "top-level sections — this is not a geckopy ecModel YAML."
        )

    ec_model = EcModel.from_cobra(
        cobra_model, adapter=adapter, gecko_light=ec.gecko_light,
    )
    ec_model.ec = ec
    return ec_model


def _resolve_path(
    filename: Optional[Union[str, Path]],
    adapter: Optional["ModelAdapter"],
) -> Path:
    if filename is None:
        filename = _DEFAULT_FILENAME
    path = Path(filename)
    if path.suffix.lower() not in _YAML_SUFFIXES:
        raise ValueError(
            "ecModels are loaded from YAML only (.yml/.yaml). Got: "
            f"{path.suffix!r}. SBML ecModel I/O has been removed."
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
