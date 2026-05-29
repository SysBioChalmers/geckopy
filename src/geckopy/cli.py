"""Command-line interface for geckopy."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from geckopy.adapter.template import generate_template_toml


ADAPTER_PY_STUB = '''\
"""Custom adapter for this ecModel.

Uncomment and customize the class below to override default behavior,
typically when your model needs a non-standard way to identify
spontaneous reactions or map gene IDs. Then import and use it:

    from myModel.adapter import MyAdapter
    adapter = MyAdapter.from_folder("myModel")

If you only need to change parameter values, edit model_adapter.toml
and use the base ModelAdapter directly.
"""
# from geckopy import ModelAdapter
#
#
# class MyAdapter(ModelAdapter):
#     def get_spontaneous_reactions(self, model):
#         # Return a list of reaction IDs for spontaneous reactions.
#         return []
'''


def cmd_init(args: argparse.Namespace) -> int:
    folder = Path(args.folder).resolve()
    if folder.exists() and any(folder.iterdir()):
        print(f"Error: {folder} already exists and is not empty.", file=sys.stderr)
        return 1

    folder.mkdir(parents=True, exist_ok=True)
    for sub in ("data", "models", "output"):
        (folder / sub).mkdir(exist_ok=True)

    toml_text = generate_template_toml(
        advanced=args.advanced, project_name=folder.name
    )
    (folder / "model_adapter.toml").write_text(toml_text)
    (folder / "adapter.py").write_text(ADAPTER_PY_STUB)

    print(f"Created ecModel project at {folder}")
    print(
        f"Next: edit {folder / 'model_adapter.toml'} to set conv_gem and "
        f"org_name, then place your GEM SBML file in {folder / 'models'}."
    )
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="geckopy", description="Command-line interface for geckopy."
    )
    sub = parser.add_subparsers(dest="command", required=True)

    init_parser = sub.add_parser(
        "init", help="Scaffold a new ecModel project folder"
    )
    init_parser.add_argument(
        "folder",
        help="Target folder (will be created; must not already exist or be empty)",
    )
    init_parser.add_argument(
        "--advanced", action="store_true",
        help="Include advanced sections (Bayesian kcat-tuning hyperparameters)",
    )
    init_parser.set_defaults(func=cmd_init)

    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
