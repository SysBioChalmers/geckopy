# %% [markdown]
# # full_ecModel tutorial
#
# Python port of the MATLAB tutorial that accompanies the GECKO 3
# Nature Protocols paper (https://doi.org/10.1038/s41596-023-00931-7).
# Demonstrates the reconstruction and analysis of a *full* ecModel
# using yeast-GEM as the starting point. STEP numbers match the
# MATLAB tutorial and the Nature Protocols paper.
#
# This script covers Stages 0-3 (model construction + tuning).
# Stages 4-5 (proteomics integration and simulation/analysis) will
# be added in subsequent phases.
#
# **Do not use the ecModel produced here outside this tutorial.**
# A maintained ecYeastGEM is distributed via
# https://github.com/SysBioChalmers/yeast-GEM (release 9.2.0+).

# %% [markdown]
# ## Setup

# %%
from pathlib import Path

import cobra

from geckopy import EcModel, ModelAdapter, make_ec_model
from geckopy.databases import (
    load_brenda_data,
    load_dlkcat_ignore_lists,
    load_phyl_dist,
    load_uniprot_tsv,
)
from geckopy.databases import find_met_smiles
from geckopy.ec_model.pipeline.apply_complex_data import apply_complex_data
from geckopy.ec_model.pipeline.apply_custom_kcats import apply_custom_kcats
from geckopy.ec_model.pipeline.apply_kcat import apply_kcat_constraints
from geckopy.ec_model.pipeline.fill_kcats import get_kcat_across_isozymes
from geckopy.ec_model.pipeline.protein_pool import set_prot_pool_size
from geckopy.ec_model.pipeline.set_kcat import set_kcat_for_reactions
from geckopy.gather_kcats.fuzzy_kcat_matching import fuzzy_kcat_matching
from geckopy.gather_kcats.get_standard_kcat import get_standard_kcat
from geckopy.gather_kcats.merge_dlkcat_and_fuzzy_kcats import (
    merge_dlkcat_and_fuzzy_kcats,
)
from geckopy.gather_kcats.read_dlkcat_output import read_dlkcat_output
from geckopy.gather_kcats.select_kcat_value import select_kcat_value
from geckopy.gather_kcats.write_dlkcat_input import write_dlkcat_input
from geckopy.get_enzyme_data.ec_from_database import get_ec_from_database
from geckopy.get_enzyme_data.ec_from_gem import get_ec_from_gem
from geckopy.utilities import (
    load_conventional_gem,
    load_ec_model,
    save_ec_model,
)

# %% [markdown]
# ## STAGE 0: Preparation
#
# **STEP 1-7** Project structure and adapter parameters are already
# in place: see `model_adapter.toml` for the organism-specific
# parameters and the `data/` and `models/` subfolders.

# %% [markdown]
# ## STAGE 1: Expansion from a starting GEM to an ecModel structure
#
# **STEP 8** Set the model adapter. geckopy has no global default
# adapter; the loaded adapter must be passed explicitly to functions
# that need it, or attached to `model.adapter`.

# %%
adapter = ModelAdapter.from_folder(Path(__file__).parent)
params = adapter.params
print(f"Organism: {params.org_name}")
print(f"Biomass reaction: {params.bio_rxn}")

# %% [markdown]
# **STEP 9** Load conventional yeast-GEM. `load_conventional_gem`
# reads the file at `adapter.params.conv_gem`.

# %%
model = load_conventional_gem(adapter)
print(f"Conventional GEM: {len(model.reactions)} reactions, "
      f"{len(model.metabolites)} metabolites, {len(model.genes)} genes")

# %% [markdown]
# **STEP 10-11** Build the ecModel. We make a *full* GECKO ecModel
# (`gecko_light=False`); see `tutorials/light_ecModel` for the light
# variant. UniProt data is loaded from `data/uniprot.tsv`.

# %%
uniprot_db = load_uniprot_tsv(params.path / "data" / "uniprot.tsv")
ec_model = make_ec_model(model, adapter, uniprot_db=uniprot_db)
print(f"ecModel: {len(ec_model.reactions)} reactions, "
      f"{len(ec_model.ec.enzymes)} enzymes")

# %% [markdown]
# **STEP 12-13** Annotate with complex data. The Complex Portal
# JSON is shipped at `data/ComplexPortal.json`; the downloader
# `geckopy.databases.get_complex_data` would refresh it.

# %%
apply_complex_data(ec_model, path=params.path / "data" / "ComplexPortal.json")

# %% [markdown]
# **STEP 14** Save Stage 1 ecModel.

# %%
save_ec_model(ec_model, "ecYeastGEM_stage1.yml", adapter=adapter)

# %% [markdown]
# ## STAGE 2: Integration of kcat into the ecModel structure
#
# To resume from disk, uncomment:
# ```python
# ec_model = load_ec_model("ecYeastGEM_stage1.yml", adapter=adapter)
# ```

# %% [markdown]
# **STEP 16-17** Gather EC numbers. First take what's annotated in
# the GEM, then fill the rest from UniProt. The MATLAB tutorial
# notes the yeast-GEM EC annotations are not thoroughly curated and
# overwrites all of them with database-derived values; we do the
# same.

# %%
get_ec_from_database(ec_model, uniprot_db)

# %% [markdown]
# **STEP 18-19** Gather kcat values from BRENDA via fuzzy matching.
# The BRENDA dumps live at `data/max_KCAT.txt` etc.; PhylDist.mat is
# the KEGG phylogenetic-distance file.
#
# Note: the PhylDist.mat shipped here is a stub. For
# better-than-tutorial results, build a real PhylDist via KEGG (a
# `get_phyl_dist` function is on the porting roadmap).

# %%
brenda = load_brenda_data(params.path / "data")
phyl_dist = load_phyl_dist(params.path / "data" / "PhylDist.mat")
kcat_list_fuzzy = fuzzy_kcat_matching(ec_model, brenda, phyl_dist)
print(f"Fuzzy BRENDA matches: {len(kcat_list_fuzzy)} rows")

# %% [markdown]
# **STEP 20-22** Gather metabolite SMILES. `find_met_smiles` queries
# PubChem; with a populated cache (`data/smilesDB.tsv`) no network
# is needed.

# %%
find_met_smiles(ec_model, cache_path=params.path / "data" / "smilesDB.tsv")

# %% [markdown]
# **STEP 23** Prepare DLKcat input file. The `data/DLKcat.tsv`
# shipped here already has predicted kcat values; uncomment to
# regenerate (this discards existing predictions).
#
# ```python
# ignore_lists = load_dlkcat_ignore_lists(params.path / "data")
# write_dlkcat_input(
#     ec_model, params.path / "data" / "DLKcat.tsv",
#     ignore_lists, overwrite=True,
# )
# ```

# %% [markdown]
# **STEP 24** Run DLKcat. (External tool; not invoked here. The
# shipped `data/DLKcat.tsv` already has output from a previous run.)

# %% [markdown]
# **STEP 25** Load DLKcat output.

# %%
kcat_list_dlkcat = read_dlkcat_output(
    ec_model, params.path / "data" / "DLKcat.tsv",
)
print(f"DLKcat predictions: {len(kcat_list_dlkcat)} rows")

# %% [markdown]
# **STEP 26** Combine kcat from BRENDA and DLKcat.

# %%
kcat_list_merged = merge_dlkcat_and_fuzzy_kcats(
    kcat_list_dlkcat, kcat_list_fuzzy,
)

# %% [markdown]
# **STEP 27** Populate `ec_model.ec.kcat` from the merged list.

# %%
select_kcat_value(ec_model, kcat_list_merged)

# %% [markdown]
# **STEP 28** Apply manually-curated kcat values from
# `data/customKcats.tsv`.

# %%
apply_custom_kcats(ec_model, path=params.path / "data" / "customKcats.tsv")

# %% [markdown]
# **STEP 29** Propagate kcat values across isozymes (sibling
# enzymes that catalyse the same reaction).

# %%
get_kcat_across_isozymes(ec_model)

# %% [markdown]
# **STEP 30** Get standard kcat. Assigns a protein cost to
# reactions without a gene association, except for exchange /
# transport / pseudoreactions and any reactions listed in
# `data/pseudoRxns.tsv`.

# %%
get_standard_kcat(ec_model, uniprot_db)

# %% [markdown]
# **STEP 31** Apply kcat constraints to the stoichiometric matrix.
# `apply_kcat_constraints` translates `ec.kcat`, `ec.mw`, and
# `ec.rxn_enz_mat` into the protein pseudo-substrate stoichiometry
# in `model.S`. Re-run any time `ec.kcat`, `ec.rxn_enz_mat`, or
# `ec.mw` changes.

# %%
apply_kcat_constraints(ec_model)

# %% [markdown]
# **STEP 32** Set the upper bound of the protein pool exchange.
# `Ptot * f * sigma` from the adapter; can be overridden per call.

# %%
set_prot_pool_size(
    ec_model,
    p_tot=params.p_tot, f=params.f, sigma=params.sigma,
)

# %%
save_ec_model(ec_model, "ecYeastGEM_stage2.yml", adapter=adapter)

# %% [markdown]
# ## STAGE 3: Model tuning
#
# **STEP 33-38** Test the maximum growth rate with an unconstrained
# glucose uptake. In MATLAB this uses RAVEN's `setParam` /
# `solveLP`; in geckopy we set bounds on `cobra.Reaction` directly
# and use `model.optimize()`.

# %%
ec_model.reactions.get_by_id(params.c_source).lower_bound = -1000
ec_model.objective = params.bio_rxn
sol = ec_model.optimize()
growth_rate = sol.fluxes[params.bio_rxn]
print(f"Growth rate: {growth_rate:.4f} /hour")
print(f"(Below the 0.41 /hour reference for S. cerevisiae)")

# %% [markdown]
# **STEP 43-44** Sensitivity tuning. Iteratively bumps the most
# limiting kcat by `fold_change` until the model can reach
# `desired_growth_rate` (defaults to `adapter.params.gr_exp`).
# Returns a result with the tuned kcats.

# %%
from geckopy.kcat_sensitivity_analysis import sensitivity_tuning

tuning_result = sensitivity_tuning(ec_model)
final_growth = ec_model.optimize().fluxes[params.bio_rxn]
print(f"Tuned {len(tuning_result.rxns)} kcats; "
      f"final growth rate: {final_growth:.4f}")

# %% [markdown]
# Note: `bayesianSensitivityTuning` (the ABC-SMC variant introduced
# in GECKO 3.3.0) is not yet ported. The MATLAB tutorial also
# skips it for the protocol walkthrough; the maintained ecYeastGEM
# in the yeast-GEM repository uses it instead.

# %% [markdown]
# **STEP 45-51** Curate kcat values based on the tuning result.
# As an example the MATLAB tutorial increases the kcat of r_0079
# (5'-phosphoribosylformyl glycinamidine synthetase) from 0.05 to a
# computed 5.34 /sec, derived from a paper-reported specific
# activity. We replicate that here.

# %%
import numpy as np

enz_idx = ec_model.ec.enzymes.index("P38972")
enz_mw = ec_model.ec.mw[enz_idx]
sa = 2.15  # umol/min/mg protein
kcat_per_sec = sa / 1000 / 60 * enz_mw
print(f"Computed kcat for r_0079: {kcat_per_sec:.4f} /sec")
set_kcat_for_reactions(ec_model, ["r_0079"], kcat_per_sec)

# %%
save_ec_model(ec_model, "ecYeastGEM_stage3.yml", adapter=adapter)

# %% [markdown]
# **STEP 52** Save the curated ecModel without proteomics
# integration (the canonical "ecYeastGEM" file).

# %%
save_ec_model(ec_model, "ecYeastGEM.yml", adapter=adapter)

print("Stage 3 complete. Stages 4-5 (proteomics integration and "
      "simulation / analysis) follow in a separate phase.")
