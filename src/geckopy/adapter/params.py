"""Parameter schemas for ecModel adapters."""
from pathlib import Path
from typing import Literal, Optional

from pydantic import BaseModel, ConfigDict, Field, model_validator


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


class SourceGroupRule(BaseModel):
    """One trust-tier entry in ``BayesianParams.source_groups``.

    ``ec.source`` holds literal strings like ``"dlkcat"``, ``"brenda"``,
    and -- for OpenKineticsPredictor kcats -- the raw predictor method
    name (e.g. ``"CataPro"``), not a generic ``"okp"`` tag. A rule
    matches a given ``ec.source`` value if it's listed in ``sources``,
    or (when ``match_okp`` is True) if it equals the project's
    configured OKP method (``OkpParams.method``).
    """
    model_config = ConfigDict(extra="forbid")

    sources: list[str] = Field(
        default_factory=list,
        description="Literal ec.source strings belonging to this group.",
    )
    match_okp: bool = Field(
        default=False,
        description=(
            "Also match ec.source values equal to the project's "
            "configured OpenKineticsPredictor method (OkpParams.method)."
        ),
    )


class BayesianParams(BaseModel):
    """Hyperparameters for Bayesian kcat fitting (ABC-SMC).

    Sources not matched by any ``source_groups`` entry fall back to the
    ``*_default`` fields (matching MATLAB's ``noKcatSource`` behaviour:
    ``sigma0logDefault`` etc. apply to every kcat first, then only ones
    with a recognised source get overridden by their group's value).
    """
    model_config = ConfigDict(extra="forbid")

    sigma0_log_default: float = 0.5
    source_groups: dict[str, SourceGroupRule] = Field(
        default_factory=lambda: {
            "dlkcat": SourceGroupRule(sources=["dlkcat"]),
            "brenda": SourceGroupRule(sources=["brenda"]),
            "custom": SourceGroupRule(sources=["custom"]),
        }
    )
    sigma0_log_source: dict[str, float] = Field(
        default_factory=lambda: {"dlkcat": 0.4, "brenda": 0.2, "custom": 0.1}
    )

    shrink_thr_default: float = 1.5
    shrink_thr_source: dict[str, float] = Field(
        default_factory=lambda: {"dlkcat": 1.5, "brenda": 3.5, "custom": 5.5}
    )

    force_prior_thr_default: float = -1.0
    force_prior_thr_source: dict[str, float] = Field(
        default_factory=lambda: {"dlkcat": -1.0, "brenda": 4.0, "custom": 8.0}
    )
    sparsity_threshold: float = 0.5

    schedule_generations: list[int] = Field(default_factory=lambda: [1, 2, 9, 15])
    schedule_samples: list[int] = Field(default_factory=lambda: [1000, 800, 600, 400])

    min_keep: float = 0.3
    max_keep: float = 0.6

    rmse_threshold: float = 0.2
    max_generations: int = 150

    max_growth_weight: float = 1.0

    # Proposal-width adaptation. Off by default: MATLAB has no such
    # feedback, so enabling it is a deliberate departure from the
    # faithful path (see docs/internal/matlab_replication_results.md).
    adapt_proposal_width: bool = False
    target_accept_rate: float = 0.15
    proposal_adaptation_rate: float = 2.0
    proposal_scale_bounds: tuple[float, float] = (0.02, 2.0)

    @model_validator(mode="after")
    def _check_group_keys_and_schedule_lengths(self) -> "BayesianParams":
        """The per-source dicts must have exactly ``source_groups``'
        keys, and the ABC-SMC schedule lists must line up with each
        other; otherwise a downstream lookup silently falls back to a
        default or a positional zip silently mismatches."""
        group_names = set(self.source_groups)
        for name in (
            "sigma0_log_source",
            "shrink_thr_source",
            "force_prior_thr_source",
        ):
            keys = set(getattr(self, name))
            if keys != group_names:
                raise ValueError(
                    f"BayesianParams.{name} keys {sorted(keys)} must match "
                    f"source_groups keys {sorted(group_names)}."
                )
        if len(self.schedule_generations) != len(self.schedule_samples):
            raise ValueError(
                "BayesianParams.schedule_generations and schedule_samples "
                "must have equal length."
            )
        return self


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

    # No range constraints: callers intentionally pass out-of-range values
    # (e.g. an inflated f) for sensitivity experiments and test fixtures.
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

    # kcat-aggregation defaults applied when the matching function is called
    # without an explicit aggregate/criteria argument. One field per function
    # so each can keep its MATLAB-GECKO default while still being centrally
    # overridable on the adapter. See docs/kcat_aggregation.md for the
    # empirical rationale for offering `median` as an alternative to `max`.
    kcat_aggregate_brenda: Literal["max", "median"] = Field(
        default="max",
        description=(
            "Default aggregation for `fuzzy_kcat_matching` when no `aggregate` "
            "is passed: collapse the BRENDA rows matched at one search level "
            "to a single kcat. `max` matches MATLAB GECKO; `median` is more "
            "robust to engineered-mutant/non-physiological-substrate tails."
        ),
    )
    kcat_aggregate_candidates: Literal["max", "min", "median", "mean"] = Field(
        default="max",
        description=(
            "Default aggregation for `apply_kcat_list` when no `criteria` is "
            "passed: collapse multiple candidate kcats for the same reaction "
            "(e.g. several EC tokens per reaction) to one. `max` matches "
            "MATLAB GECKO."
        ),
    )
    kcat_aggregate_isozymes: Literal["max", "median", "mean"] = Field(
        default="mean",
        description=(
            "Default aggregation for `fill_kcats_from_isozymes` when no "
            "`aggregate` is passed: combine known kcats across sibling "
            "isozymes to fill a missing entry. `mean` matches MATLAB GECKO."
        ),
    )

    kegg: KeggParams = Field(default_factory=KeggParams)
    uniprot: UniprotParams = Field(default_factory=UniprotParams)
    complex: ComplexParams = Field(default_factory=ComplexParams)
    bayesian: BayesianParams = Field(default_factory=BayesianParams)
    okp: OkpParams = Field(default_factory=OkpParams)
