"""Generate commented TOML templates from the adapter schema.

The pydantic schema in `params.py` is the single source of truth for
parameter names, types, defaults, and descriptions. This module walks
that schema and emits a TOML file where every field is either:
uncommented with a placeholder (for required fields), or commented out
showing its default value (for optional fields). Field descriptions
are emitted as TOML comments above each field.
"""
from __future__ import annotations

import inspect
from pathlib import Path
from typing import Any, Literal, get_args, get_origin

from pydantic import BaseModel

from .params import ModelParameters

_INJECTED_FIELDS = {"path"}  # loader injects at runtime, never in template
_ADVANCED_SECTIONS = {"bayesian"}
_REQUIRED_PLACEHOLDERS = {
    "conv_gem": '"models/REPLACE_ME.xml"',
    "org_name": '"REPLACE_ME"',
}


def _format_toml_value(value: Any) -> str | None:
    """Format a Python value as a TOML literal, or None if unrepresentable."""
    if value is None:
        return None
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, str):
        return f'"{value}"'
    if isinstance(value, Path):
        return f'"{value.as_posix()}"'
    if isinstance(value, (int, float)):
        return repr(value)
    if isinstance(value, list):
        parts = [_format_toml_value(v) for v in value]
        if any(p is None for p in parts):
            return None
        return "[" + ", ".join(parts) + "]"
    return None


def _annotation_hint(annotation: Any) -> str:
    """Extra comment hint; currently only for Literal types."""
    if get_origin(annotation) is Literal:
        choices = ", ".join(repr(a) for a in get_args(annotation))
        return f"One of: {choices}"
    return ""


def _emit_scalar(
    name: str, *,
    description: str | None,
    annotation: Any,
    default_value: Any,
    required: bool,
    placeholder: str | None = None,
) -> list[str]:
    lines: list[str] = []
    if description:
        lines.append(f"# {description}")
    hint = _annotation_hint(annotation)
    if hint:
        lines.append(f"# {hint}")

    if required:
        assert placeholder is not None
        lines.append(f"{name} = {placeholder}")
    else:
        formatted = _format_toml_value(default_value)
        if formatted is None:
            lines.append(f"# {name} =")
        else:
            lines.append(f"# {name} = {formatted}")
    lines.append("")
    return lines


def _is_nested_model(annotation: Any) -> bool:
    return inspect.isclass(annotation) and issubclass(annotation, BaseModel)


def generate_template_toml(
    *, advanced: bool = False, project_name: str = "my-ecmodel",
) -> str:
    """Generate a commented TOML template for a new ecModel project.

    Parameters
    ----------
    advanced
        Include the advanced sections (currently just Bayesian
        kcat-tuning hyperparameters); omitted by default since most
        projects never need to touch them.
    project_name
        Name written into the template's header comment.

    Returns
    -------
    str
        The full TOML template text, ending in a newline.
    """
    lines: list[str] = [
        f"# Adapter configuration for: {project_name}",
        "# Uncomment and edit lines below to override defaults.",
        "# Required fields are already uncommented.",
        "",
    ]

    top_level: list[tuple[str, Any]] = []
    sections: list[tuple[str, type[BaseModel]]] = []

    for name, field in ModelParameters.model_fields.items():
        if name in _INJECTED_FIELDS:
            continue
        if _is_nested_model(field.annotation):
            sections.append((name, field.annotation))
        else:
            top_level.append((name, field))

    # TOML requires top-level keys before any [section].
    for name, field in top_level:
        if field.is_required():
            placeholder = _REQUIRED_PLACEHOLDERS.get(name, '"..."')
            lines.extend(_emit_scalar(
                name, description=field.description,
                annotation=field.annotation, default_value=None,
                required=True, placeholder=placeholder,
            ))
        else:
            lines.extend(_emit_scalar(
                name, description=field.description,
                annotation=field.annotation, default_value=field.default,
                required=False,
            ))

    has_omitted_advanced = False
    for name, model_cls in sections:
        if name in _ADVANCED_SECTIONS and not advanced:
            has_omitted_advanced = True
            continue
        lines.append(f"[{name}]")
        instance = model_cls()
        for sub_name, sub_field in model_cls.model_fields.items():
            lines.extend(_emit_scalar(
                sub_name, description=sub_field.description,
                annotation=sub_field.annotation,
                default_value=getattr(instance, sub_name),
                required=False,
            ))

    if has_omitted_advanced:
        lines.extend([
            "# Note: advanced sections (Bayesian kcat-tuning hyperparameters)",
            "# are omitted by default. Regenerate with --advanced to include them.",
        ])

    return "\n".join(lines) + "\n"
