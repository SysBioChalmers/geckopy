"""Parameter schemas for ecModel adapters."""
from pathlib import Path
from typing import Literal, Optional

from pydantic import BaseModel, ConfigDict, Field


class KeggParams(BaseModel):
    model_config = ConfigDict(extra="forbid")

    id: str = Field(default="sce", description="KEGG organism code")
    gene_id: str = Field(
        default="kegg",
        description=(
            "Field in the KEGG entry matching model gene IDs. "
            "'kegg' uses the default KEGG entry identifier; other options "
            "include 'NCBI-GeneID', 'UniProt', 'Ensembl'."
        ),
    )


class UniprotParams(BaseModel):
    model_config = ConfigDict(extra="forbid")

    type: Literal["taxonomy", "proteome"] = Field(
        default="taxonomy",
        description="Whether 'id' is a taxonomy ID or a proteome ID",
    )
    id: str = Field(
        default="559292",
        description="Taxonomy ID (e.g. '559292') or proteome ID, matching 'type'",
    )
    gene_id_field: str = Field(
        default="gene_oln",
        description="UniProt return field matching model gene IDs",
    )
    reviewed: bool = Field(
        default=False,
        description="Restrict queries to reviewed (Swiss-Prot) entries",
    )


class ComplexParams(BaseModel):
    model_config = ConfigDict(extra="forbid")

    taxonomic_id: Optional[int] = Field(
        default=None, description="Taxonomy ID for Complex Portal query"
    )


class BayesianParams(BaseModel):
    """Hyperparameters for Bayesian kcat fitting (ABC-SMC)."""
    model_config = ConfigDict(extra="forbid")

    sigma0_log_default: float = 0.5
    kcat_sources: list[str] = Field(
        default_factory=lambda: ["dlkcat", "brenda", "custom"]
    )
    sigma0_log_source: list[float] = Field(default_factory=lambda: [0.4, 0.2, 0.1])

    shrink_thr_default: float = 1.5
    shrink_thr_source: list[float] = Field(default_factory=lambda: [1.5, 3.5, 5.5])
    variance_cap_default: float = 10.0
    variance_cap_source: list[float] = Field(default_factory=lambda: [10.0, 4.0, 2.0])

    force_prior_thr_default: float = -1.0
    force_prior_thr_source: list[float] = Field(
        default_factory=lambda: [-1.0, 4.0, 8.0]
    )
    sparsity_threshold: float = 0.3

    schedule_generations: list[int] = Field(default_factory=lambda: [1, 2, 9, 15])
    schedule_samples: list[int] = Field(default_factory=lambda: [1000, 800, 600, 400])

    target_accept: float = 10.0
    min_keep: float = 0.3
    max_keep: float = 0.6

    rmse_threshold: float = 0.2
    max_generations: int = 150


class OkpParams(BaseModel):
    """OpenKineticsPredictor settings (used by submit/fetch_open_kinetics_predictor).

    The API key is intentionally NOT here (it is a secret): provide it as a
    function argument, the OKP_API_KEY environment variable, or the file
    ``<path>/data/okpApiKey.txt``.
    """
    model_config = ConfigDict(extra="forbid")

    method: str = Field(
        default="CataPro",
        description=(
            "Predictor method: CataPro, CatPred, DLKcat, EITLEM, "
            "KinForm-H, KinForm-L, UniKP (see GET /api/v1/methods/)."
        ),
    )
    handle_long_sequences: str = Field(
        default="truncate",
        description="How OKP handles sequences exceeding a method's max length.",
    )
    include_similarity_columns: bool = Field(
        default=True,
        description="Append per-row similarity-to-training-data columns.",
    )
    canonicalize_substrates: bool = Field(
        default=True,
        description="Canonicalize substrate SMILES server-side before prediction.",
    )


class ModelParameters(BaseModel):
    """Top-level parameters for an ecModel adapter."""
    model_config = ConfigDict(extra="forbid")

    path: Path = Field(description="Root folder of this ecModel project")
    conv_gem: Path = Field(description="Path to the conventional GEM SBML file")

    org_name: str = Field(description="Scientific name of the organism")

    sigma: float = Field(default=0.5, description="Average enzyme saturation factor")
    p_tot: float = Field(default=0.5, description="Total protein content [g/gDw]")
    f: float = Field(default=0.5, description="Fraction of enzymes in model [g/g]")
    gr_exp: float = Field(default=0.41, description="Reference growth rate [1/h]")

    c_source: str = Field(
        default="", description="Reaction ID for preferred carbon source exchange"
    )
    bio_rxn: str = Field(
        default="", description="Reaction ID for the biomass pseudoreaction"
    )
    enzyme_comp: str = Field(
        default="cytoplasm",
        description="Compartment where protein pseudometabolites are placed",
    )

    kegg: KeggParams = Field(default_factory=KeggParams)
    uniprot: UniprotParams = Field(default_factory=UniprotParams)
    complex: ComplexParams = Field(default_factory=ComplexParams)
    bayesian: BayesianParams = Field(default_factory=BayesianParams)
    okp: OkpParams = Field(default_factory=OkpParams)
