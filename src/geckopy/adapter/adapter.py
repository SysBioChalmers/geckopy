"""Project-level configuration for a geckopy ecModel.

A ``ModelAdapter`` bundles everything that's specific to *this*
ecModel project: the organism's name and taxonomy id, the path
to the starting GEM, the path to UniProt and other cached data,
the biomass reaction, the assumed enzyme saturation factor, and so
on. Downstream geckopy functions take an adapter as input rather
than dozens of separate keyword arguments.

A project lives in a folder containing:

- ``model_adapter.toml`` — the parameters (see
  ``geckopy.adapter.params.ModelParameters`` for the schema).
- ``models/`` — the starting GEM and any ecModel outputs.
- ``data/`` — cached databases (UniProt, BRENDA, ComplexPortal,
  proteomics, etc.).

Then:

.. code-block:: python

    adapter = ModelAdapter.from_folder("my_project")

For most users that's all you need. Subclass ``ModelAdapter``
only when an organism needs custom logic — typically a custom
``get_spontaneous_reactions`` or a custom gene-id mapping.
"""
from __future__ import annotations

import tomllib
from pathlib import Path
from typing import TYPE_CHECKING

from .params import ModelParameters

if TYPE_CHECKING:
    import cobra


class ModelAdapter:
    """Holds the organism-specific parameters and overridable
    helpers for one ecModel project.

    Most users don't subclass this directly: they make a folder
    with a ``model_adapter.toml`` and call
    ``ModelAdapter.from_folder(folder)``.

    Subclass only when a particular organism needs custom logic
    (e.g. an unusual way to identify spontaneous reactions, or a
    bespoke gene-id mapping that the default UniProt-table-driven
    one can't handle).
    """

    def __init__(self, params: ModelParameters):
        self.params = params

    @classmethod
    def from_folder(cls, folder: str | Path) -> "ModelAdapter":
        """Load an adapter from a folder containing model_adapter.toml.

        Parameters
        ----------
        folder
            Path to the ecModel project folder. Must contain
            `model_adapter.toml`.

        Returns
        -------
        A ModelAdapter (or subclass) instance with parameters loaded
        from the TOML file.
        """
        folder = Path(folder).resolve()
        config_path = folder / "model_adapter.toml"
        if not config_path.is_file():
            raise FileNotFoundError(
                f"No model_adapter.toml found in {folder}. "
                "See docs for the expected folder structure."
            )

        with open(config_path, "rb") as f:
            data = tomllib.load(f)

        # The path field is injected from the folder location, not the TOML,
        # so that the TOML file is location-independent.
        data["path"] = folder
        params = ModelParameters.model_validate(data)

        # Resolve conv_gem relative to the folder if it is not absolute.
        if not params.conv_gem.is_absolute():
            params.conv_gem = folder / params.conv_gem

        # Automatic discovery of a subclass in `adapter.py` is intentionally
        # not enabled by default. Users should import their subclass
        # explicitly. To enable auto-discovery, uncomment the block below.
        #
        # adapter_py = folder / "adapter.py"
        # if adapter_py.is_file() and cls is ModelAdapter:
        #     import importlib.util
        #     spec = importlib.util.spec_from_file_location(
        #         f"_adapter_{folder.name}", adapter_py
        #     )
        #     module = importlib.util.module_from_spec(spec)
        #     spec.loader.exec_module(module)
        #     subclasses = [
        #         v for v in vars(module).values()
        #         if isinstance(v, type)
        #         and issubclass(v, ModelAdapter)
        #         and v is not ModelAdapter
        #     ]
        #     if len(subclasses) == 1:
        #         return subclasses[0](params)
        #     elif len(subclasses) > 1:
        #         raise RuntimeError(
        #             f"Multiple ModelAdapter subclasses found in {adapter_py}. "
        #             "Import the one you want explicitly."
        #         )

        return cls(params)

    # Overridable methods. The default implementations match the MATLAB
    # base class behavior. Subclass and override as needed.

    def get_spontaneous_reactions(self, model: "cobra.Model") -> list[str]:
        """Return IDs of reactions that proceed without enzyme catalysis.

        The default returns an empty list. Override in a subclass to
        identify spontaneous reactions for your organism.
        """
        return []

    def get_uniprot_compatible_genes(self, in_genes: list[str]) -> list[str]:
        """Map model gene IDs to IDs usable in UniProt queries.

        The default returns the input unchanged. Override when the model
        uses gene identifiers that need transformation before UniProt
        lookup (for example, stripping a prefix).
        """
        return list(in_genes)

    def get_uniprot_ids_from_table(self, model_genes: list[str]) -> list[str]:
        """Map model gene IDs to UniProt IDs via data/uniprotConversion.tsv.

        If the conversion table is absent, returns the input unchanged.
        The table is a two-column TSV with header row: model gene ID,
        UniProt ID.
        """
        table_path = self.params.path / "data" / "uniprotConversion.tsv"
        if not table_path.is_file():
            return list(model_genes)

        mapping: dict[str, str] = {}
        with open(table_path, "r", encoding="utf-8") as f:
            next(f)  # skip header
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 2:
                    mapping[parts[0]] = parts[1]

        return [mapping.get(g, g) for g in model_genes]

    def get_brenda_db_folder(self) -> Path:
        """Return the folder where the BRENDA database files are stored.

        Default: the BRENDA snapshot bundled with geckopy
        (`<geckopy>/data/brenda/`), refreshed via `geckopy brenda-refresh`.
        Override to point at a project-local copy.
        """
        return Path(__file__).resolve().parent.parent / "data" / "brenda"

    def get_phyl_dist_path(self) -> Path:
        """Return the path to the phylogenetic distance structure file.

        Default: `<project_path>/data/PhylDist.mat`.
        """
        return self.params.path / "data" / "PhylDist.mat"
