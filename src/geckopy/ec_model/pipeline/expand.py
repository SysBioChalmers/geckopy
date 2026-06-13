"""Expand reactions with isozymes into one reaction per isozyme.

Corresponds to stage 5 of makeEcModel in GECKO MATLAB. Only applies to full
ecModels (geckoLight skips this stage).

This implementation was originally written for geckopy and then adopted as the
canonical home in raven-toolbox (which uses cobrapy's GPR AST instead of
RAVEN's string manipulation, matching the geckopy version). geckopy re-exports
it from raven-toolbox under the same name so the rest of the pipeline is
unaffected.
"""
from __future__ import annotations

# Re-export from raven-toolbox — see module docstring for the migration history.
from raven_toolbox.manipulation.expand import expand_model

__all__ = ["expand_model"]
