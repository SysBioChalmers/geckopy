"""Parameter schemas for ecModel adapters.

These pydantic models mirror the structure of the MATLAB adapter's
`obj.params` struct, with snake_case field names per Python convention.
The schema is the single source of truth for what a valid adapter
configuration looks like, and pydantic validates TOML files against it.
"""
from pathlib import Path
from typing import Optional
from pydantic import BaseModel, Field, ConfigDict


class KeggParams(BaseModel):
    """KEGG database lookup parameters."""
    model_config = ConfigDict(extra="forbid")

    id: str = Field(
        default="sce",
        description="KEGG organism code, see https://www.genome.jp/kegg/catalog/org_list.html",
    )
    gene_id_field: str = Field(
        default="kegg",
        description="Field in the KEGG entry that matches the model's gene IDs",
    )


class UniprotParams(BaseModel):
    """UniProt lookup parameters."""
    model_config = ConfigDict(extra="forbid")

    tax_id: Optional[str] = Field(
        default=None,
        description="NCBI taxonomic ID for UniProt query",
    )
    reviewed: bool = Field(
        default=False,
        description="Restrict UniProt query to reviewed (Swiss-Prot) entries only",
    )
    gene_id_field: str = Field(
        default="gene_names",
        description="UniProt field used to match model gene IDs",
    )


class ComplexParams(BaseModel):
    """Complex Portal lookup parameters."""
    model_config = ConfigDict(extra="forbid")

    taxonomic_id: Optional[int] = Field(
        default=None,
        description="Taxonomic ID for Complex Portal query",
    )


class ModelParameters(BaseModel):
    """Top-level parameters for an ecModel adapter.

    All paths are resolved relative to `path` (the adapter folder) unless
    given as absolute paths. `path` itself is set automatically by
    ModelAdapter.from_folder().
    """
    model_config = ConfigDict(extra="forbid")

    path: Path = Field(description="Root folder of this ecModel project")
    conv_gem: Path = Field(description="Path to the conventional GEM (SBML) file")

    org_name: str = Field(description="Scientific name of the organism")

    sigma: float = Field(default=0.5, description="Average enzyme saturation factor")
    p_tot: float = Field(default=0.5, description="Total protein content [g protein/gDw]")
    f: float = Field(default=0.5, description="Fraction of enzymes in the model [g/g]")
    gr_exp: float = Field(default=0.41, description="Reference growth rate [1/h]")

    kegg: KeggParams = Field(default_factory=KeggParams)
    uniprot: UniprotParams = Field(default_factory=UniprotParams)
    complex: ComplexParams = Field(default_factory=ComplexParams)