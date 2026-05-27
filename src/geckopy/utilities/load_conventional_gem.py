"""Load the conventional (non-ec) starting GEM for a project.

Every geckopy project has a "starting GEM" — a plain cobra model
without any GECKO additions. This is the model that
``make_ec_model`` extends into an ecModel.

This helper reads that file from the path set in the adapter
(``adapter.params.conv_gem``). The extension picks the reader:
``.yml``/``.yaml`` via ``load_yaml_model``, ``.json`` via
``load_json_model``, ``.mat`` via ``load_matlab_model``, and
anything else via ``read_sbml_model``.

Ported from GECKO MATLAB:
src/geckomat/utilities/loadConventionalGEM.m.
"""
from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import cobra

if TYPE_CHECKING:
    from ..adapter import ModelAdapter


_YAML_SUFFIXES = {".yml", ".yaml"}
_JSON_SUFFIXES = {".json"}
_MAT_SUFFIXES = {".mat"}


def load_conventional_gem(adapter: "ModelAdapter") -> cobra.Model:
    """Load the starting GEM file pointed to by ``adapter``.

    Parameters
    ----------
    adapter
        A ``ModelAdapter`` whose ``params.conv_gem`` field gives
        the path to the GEM file. The path is used as-is; the
        adapter's loader has already resolved it relative to the
        adapter folder.

    Returns
    -------
    cobra.Model
        The loaded model. Use it as the input to
        ``make_ec_model``.

    Raises
    ------
    FileNotFoundError
        If ``params.conv_gem`` doesn't point to an existing file.
    """
    path = Path(adapter.params.conv_gem)
    if not path.is_file():
        raise FileNotFoundError(f"conv_gem file not found: {path}")
    suffix = path.suffix.lower()
    if suffix in _YAML_SUFFIXES:
        return cobra.io.load_yaml_model(str(path))
    if suffix in _JSON_SUFFIXES:
        return cobra.io.load_json_model(str(path))
    if suffix in _MAT_SUFFIXES:
        return cobra.io.load_matlab_model(str(path))
    return cobra.io.read_sbml_model(str(path))
