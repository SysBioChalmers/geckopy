"""Anaerobic switch for the yeast-GEM ecModel.

Port of `GECKO/tutorials/full_ecModel/code/anaerobicModel_GECKO.m`,
which `YeastGEMAdapter.makeModelAnaerobic` calls for every Bayesian
tuning condition whose oxygen exchange is measured as exactly 0.

The changes are organism-specific, so they live with the tutorial
rather than in geckopy, mirroring where GECKO keeps them. Pass
:func:`make_anaerobic` to `bayesian_kcat_tuning(make_anaerobic=...)`;
it is a module-level function so it survives pickling to the scoring
pool's workers.

Every change goes through the cobra API while the caller holds the
model's context (`with model:`), so all of them revert when the
condition's simulation ends.
"""
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from geckopy import EcModel

#: Uptakes opened for anaerobic growth: sterols and fatty acids cannot
#: be synthesised without oxygen, and the vitamins feed NAD(P)H and CoA
#: synthesis.
_OPEN_UPTAKE = (
    "r_1757",  # ergosterol
    "r_1915",  # lanosterol
    "r_2106",  # zymosterol
    "r_2134",  # 14-demethyllanosterol
    "r_1994",  # palmitoleate
    "r_2189",  # oleate
    "r_1967",  # nicotinate
    "r_1548",  # (R)-pantothenate
)

#: Blocked outright. MDH2 and IDP2 are repressed on glucose and absent
#: from the proteome; the sterol uptake is blocked because it recycles
#: NADH to ergosterol.
_BLOCKED = (
    "r_2137",  # ergosta-5,7,22,24(28)-tetraen-3beta-ol
    "r_0714",  # MDH2, malate dehydrogenase (cytoplasmic)
    "r_0659",  # IDP2, isocitrate dehydrogenase (NADP)
)

#: FADH2 turnover from disulphide bond formation via Ero1, added to the
#: biomass pseudoreaction.
_FADH2_PROD = 0.08


def make_anaerobic(model: "EcModel") -> None:
    """Constrain a yeast-GEM ecModel to anaerobic conditions, in place."""
    # Heme a leaves the cofactor pseudoreaction: its synthesis needs O2.
    cofactor = model.reactions.get_by_id("r_4598")
    heme = model.metabolites.get_by_id("s_3714")
    if heme in cofactor.metabolites:
        cofactor.add_metabolites({heme: -cofactor.metabolites[heme]})

    model.reactions.get_by_id("r_1992").lower_bound = 0.0  # O2 uptake
    for rxn_id in _OPEN_UPTAKE:
        model.reactions.get_by_id(rxn_id).lower_bound = -1000.0
    for rxn_id in _BLOCKED:
        # MATLAB's setParam('eq', ...) matches the reaction ID exactly,
        # so the `_REV` half of a split reversible pair is left alone.
        model.reactions.get_by_id(rxn_id).bounds = (0.0, 0.0)

    biomass = model.reactions.get_by_id("r_4041")
    biomass.add_metabolites({
        model.metabolites.get_by_id("s_0689"): _FADH2_PROD,       # FADH2[c]
        model.metabolites.get_by_id("s_0687"): -_FADH2_PROD,      # FAD[c]
        model.metabolites.get_by_id("s_0794"): -2 * _FADH2_PROD,  # H+[c]
    })
