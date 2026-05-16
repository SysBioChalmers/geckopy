"""Load the conventional (non-ec) GEM referenced by an adapter.

Ported from GECKO MATLAB:
src/geckomat/utilities/loadConventionalGEM.m.

Reads the file at ``adapter.params.conv_gem``, dispatching on
extension: YAML via :func:`cobra.io.load_yaml_model`, anything
else via :func:`cobra.io.read_sbml_model`.
"""
from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import cobra

if TYPE_CHECKING:
    from ..adapter import ModelAdapter


_YAML_SUFFIXES = {".yml", ".yaml"}


def load_conventional_gem(adapter: "ModelAdapter") -> cobra.Model:
    """Load the conventional GEM file pointed to by ``adapter``.

    Parameters
    ----------
    adapter
        ModelAdapter whose ``params.conv_gem`` points at the GEM
        file. Resolved as-is (the adapter's loader has already
        absolutised the path relative to the adapter folder).

    Returns
    -------
    cobra.Model
        The loaded model.

    Raises
    ------
    FileNotFoundError
        If ``params.conv_gem`` does not point to an existing file.
    """
    path = Path(adapter.params.conv_gem)
    if not path.is_file():
        raise FileNotFoundError(f"conv_gem file not found: {path}")
    if path.suffix.lower() in _YAML_SUFFIXES:
        return cobra.io.load_yaml_model(str(path))
    return cobra.io.read_sbml_model(str(path))
