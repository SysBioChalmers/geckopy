"""SBML reader and writer for ecModels.

SBML is the standard exchange format for constraint-based
modelling. Writing an ecModel as SBML means other tools (MATLAB
GECKO, Memote, escher, BiGG, ...) can read at least the cobra
portion. Tools that don't understand the GECKO extensions still
see a self-consistent FBA model — the enzyme mass-balance
constraints are explicit in the stoichiometry.

The encoding (the "MW_KCAT" format, GECKO 3's default) is:

- Each enzyme is a regular metabolite named ``prot_<uniprot>``.
- These metabolites are grouped in an SBML ``Group`` named
  ``Protein`` so an SBML-aware reader can find them.
- Molecular weight is carried in the species notes (cobra's
  ``metabolite.notes`` dict).
- Proteomics concentration goes in the species
  ``initialAmount``.
- The stoichiometric coefficient of ``prot_<id>`` on a catalysed
  reaction is ``-MW / (3600 * kcat)`` per subunit, matching the
  GECKO 3 MATLAB convention exactly so no rescaling is needed at
  round-trip.

Ported from the legacy geckopy package described in Carrasco et al.
(2023, https://doi.org/10.1128/spectrum.01705-23), file
geckopy/io/sbml.py. Simplified to a single encoding (no
``EcStoichiometry`` enum, no dual KCAT mode).
"""
from __future__ import annotations

import logging
import re
from pathlib import Path
from typing import TYPE_CHECKING, Optional

import cobra
import libsbml
import numpy as np
from scipy import sparse

from ..ec_model.constants import POOL_EXCHANGE_ID, POOL_ID, PROT_PREFIX, USAGE_PREFIX

if TYPE_CHECKING:
    from ..adapter.adapter import ModelAdapter
    from ..ec_model.ec_model import EcModel

logger = logging.getLogger(__name__)

_PROTEIN_GROUP = "Protein"
_COBRA_MET_PREFIX = "M_"


# --------------------------------------------------------------------------- #
# Write
# --------------------------------------------------------------------------- #

def write_sbml_ec_model(model: "EcModel", filename: str | Path) -> None:
    """Write an EcModel as SBML with enzymes in a ``Protein`` Group.

    Coefficients on ``prot_<id>`` metabolites are scaled to
    ``-MW/(3600*kcat)`` so tools reading the file without geckopy
    get a self-consistent FBA model where enzyme balances are
    explicit. The Protein group + per-species MW notes are the
    GECKO-specific extensions and are silently ignored by such
    tools.

    Ported from the legacy geckopy package (Carrasco et al., 2023,
    https://doi.org/10.1128/spectrum.01705-23),
    geckopy/io/sbml.py:799-1190
    (write_sbml_ec_model), simplified to MW_KCAT encoding only.
    """
    filename = str(filename)

    # 1. Tag prot_<id> metabolites with their MW (cobra serialises
    #    metabolite.notes into SBML <notes>; this carries through
    #    the round-trip cleanly without fighting libsbml.appendNotes).
    annotated = _annotate_mw(model)

    # 2. Write via cobra (handles all standard SBML). The in-memory
    #    coefficient on prot_<id> mets is already MW_KCAT-encoded
    #    (-subunits * MW / (kcat * 3600)), matching the on-disk
    #    convention used by MATLAB GECKO, so no rescaling needed.
    cobra.io.write_sbml_model(annotated, filename)

    # 3. Inject Protein group + MW notes via libsbml.
    reader = libsbml.SBMLReader()
    doc = reader.readSBML(filename)
    sbml_model = doc.getModel()

    # Enable Groups extension (cobra usually loads it).
    if not doc.isPackageEnabled("groups"):
        doc.enablePackage(
            "http://www.sbml.org/sbml/level3/version1/groups/version1",
            "groups",
            True,
        )
    groups_plugin = sbml_model.getPlugin("groups")
    if groups_plugin is None:
        raise RuntimeError(
            "libsbml groups package not available; cannot write Protein group"
        )

    group = groups_plugin.createGroup()
    group.setId(_PROTEIN_GROUP)
    group.setName(_PROTEIN_GROUP)
    group.setKind("collection")

    for u in model.ec.enzymes:
        species_id = f"{_COBRA_MET_PREFIX}{PROT_PREFIX}{u}"
        species = sbml_model.getSpecies(species_id)
        if species is None:
            continue
        idx = model.ec.enzymes.index(u)
        conc = float(model.ec.concs[idx])
        species.setInitialAmount(0.0 if np.isnan(conc) else conc)
        # MW is written into the species notes via cobra's metabolite.notes
        # dict in _scale_for_sbml (avoids libsbml.appendNotes returning -5).
        member = group.createMember()
        member.setIdRef(species_id)

    libsbml.writeSBMLToFile(doc, filename)


def _annotate_mw(model: "EcModel") -> "EcModel":
    """Return a deep-copy of ``model`` with each prot_<id>
    metabolite's MW written into its ``notes`` dict.

    cobra serialises ``metabolite.notes`` into SBML ``<notes>`` as
    ``<p>key: value</p>`` entries, which both cobra (on read) and
    geckopy's reader pick up as ``met.notes["mw"]``. This avoids
    the libsbml ``appendNotes`` quirk (returns status -5 even for
    well-formed XHTML).
    """
    import copy as _copy

    out = model.copy()
    out.ec = _copy.deepcopy(model.ec)
    for i, enz in enumerate(model.ec.enzymes):
        mw = float(model.ec.mw[i])
        if not np.isfinite(mw):
            continue
        met_id = f"{PROT_PREFIX}{enz}"
        if met_id in {m.id for m in out.metabolites}:
            out.metabolites.get_by_id(met_id).notes["mw"] = str(mw)
    return out


# --------------------------------------------------------------------------- #
# Read
# --------------------------------------------------------------------------- #

def read_sbml_ec_model(
    filename: str | Path,
    *,
    adapter: Optional["ModelAdapter"] = None,
) -> "EcModel":
    """Read an SBML file written with MW_KCAT encoding back into an EcModel.

    Reconstructs ``model.ec`` from Protein group membership, species
    notes (MW), and reaction stoichiometry. Inverse of
    ``write_sbml_ec_model``; pure-FBA SBMLs without the Protein
    group come back as an EcModel with an empty ``ec`` substructure.

    Parameters
    ----------
    filename
        Path to the SBML file.
    adapter
        ModelAdapter to attach to the returned EcModel. Required;
        downstream functions need it for organism parameters.

    Ported from the legacy geckopy package (Carrasco et al., 2023,
    https://doi.org/10.1128/spectrum.01705-23),
    geckopy/io/sbml.py:95-700
    (read_sbml_ec_model), simplified to MW_KCAT only.
    """
    from ..ec_model.ec_data import EcData
    from ..ec_model.ec_model import EcModel

    filename = str(filename)
    if adapter is None:
        raise ValueError("read_sbml_ec_model requires an adapter")

    # 1. Standard cobra read.
    cobra_model = cobra.io.read_sbml_model(filename)

    # 2. Walk libsbml document for the Protein group.
    reader = libsbml.SBMLReader()
    doc = reader.readSBML(filename)
    sbml_model = doc.getModel()
    groups_plugin = sbml_model.getPlugin("groups")

    enzyme_ids: list[str] = []
    enzyme_mw: dict[str, float] = {}
    enzyme_conc: dict[str, float] = {}
    cobra_met_ids = {m.id for m in cobra_model.metabolites}
    if groups_plugin is not None:
        for gi in range(groups_plugin.getNumGroups()):
            g = groups_plugin.getGroup(gi)
            if g.getId() != _PROTEIN_GROUP:
                continue
            for mi in range(g.getNumMembers()):
                member = g.getMember(mi)
                sid = member.getIdRef()
                species = sbml_model.getSpecies(sid)
                if species is None:
                    continue
                clean = (
                    sid[len(_COBRA_MET_PREFIX):]
                    if sid.startswith(_COBRA_MET_PREFIX) else sid
                )
                if not clean.startswith(PROT_PREFIX):
                    continue
                uniprot = clean[len(PROT_PREFIX):]
                if uniprot == POOL_ID[len(PROT_PREFIX):]:
                    continue
                enzyme_ids.append(uniprot)
                # Prefer cobra's parsed notes (reliable); fall back to
                # raw libsbml notes (works for hand-written SBMLs that
                # don't follow cobra's `<p>key: value</p>` convention).
                mw = None
                if clean in cobra_met_ids:
                    notes = cobra_model.metabolites.get_by_id(clean).notes
                    if "mw" in notes:
                        try:
                            mw = float(notes["mw"])
                        except (TypeError, ValueError):
                            mw = None
                if mw is None:
                    mw = _parse_mw_from_notes(species.getNotesString())
                if mw is not None:
                    enzyme_mw[uniprot] = mw
                amt = species.getInitialAmount()
                enzyme_conc[uniprot] = float(amt) if amt > 0 else float("nan")

    # 3. Wrap and populate.
    ec_model = EcModel.from_cobra(cobra_model, adapter)
    if not enzyme_ids:
        ec_model.ec = EcData()
        return ec_model

    _populate_ec_from_sbml(
        ec_model, enzyme_ids, enzyme_mw, enzyme_conc,
    )
    return ec_model


_MW_NOTE_RE = re.compile(r"mw:\s*([0-9.eE+\-]+)")


def _parse_mw_from_notes(notes_string: str) -> Optional[float]:
    if not notes_string:
        return None
    m = _MW_NOTE_RE.search(notes_string)
    return float(m.group(1)) if m else None


def _populate_ec_from_sbml(
    ec_model: "EcModel",
    enzyme_ids: list[str],
    enzyme_mw: dict[str, float],
    enzyme_conc: dict[str, float],
) -> None:
    """Reconstruct ``model.ec`` from per-enzyme metadata + reaction
    stoichiometry.

    For each metabolic reaction that has a ``prot_<id>`` in its
    stoichiometry, derive kcat from the coefficient using
    ``kcat = -mw / (coef * 3600)``. Subunits are assumed 1; the
    coefficient already encodes the per-reaction stoichiometry.

    Ported in concept from the legacy geckopy package
    (Carrasco et al., 2023, https://doi.org/10.1128/spectrum.01705-23),
    geckopy/io/sbml.py:540-787.
    """
    from ..ec_model.ec_data import EcData

    # Per-enzyme arrays (order preserved from group iteration).
    n_e = len(enzyme_ids)
    mw_arr = np.array(
        [enzyme_mw.get(u, np.nan) for u in enzyme_ids], dtype=float,
    )
    concs_arr = np.array(
        [enzyme_conc.get(u, np.nan) for u in enzyme_ids], dtype=float,
    )
    enzyme_index = {u: i for i, u in enumerate(enzyme_ids)}

    # Walk metabolic reactions; collect (rxn_id, enzyme, coef) tuples.
    catalysed_rxns: list[str] = []
    rxn_index: dict[str, int] = {}
    triples: list[tuple[int, int, float]] = []  # (rxn_idx, enz_idx, kcat)

    for rxn in ec_model.reactions:
        if rxn.id == POOL_EXCHANGE_ID or rxn.id.startswith(USAGE_PREFIX):
            continue
        for met, coef in rxn.metabolites.items():
            if met.id == POOL_ID or not met.id.startswith(PROT_PREFIX):
                continue
            enzyme = met.id[len(PROT_PREFIX):]
            j = enzyme_index.get(enzyme)
            if j is None:
                continue
            mw = mw_arr[j]
            if not np.isfinite(mw) or coef == 0:
                continue
            # coef = -mw / (3600 * kcat)  =>  kcat = -mw / (3600 * coef)
            kcat = float(-mw / (3600.0 * coef))
            if rxn.id not in rxn_index:
                rxn_index[rxn.id] = len(catalysed_rxns)
                catalysed_rxns.append(rxn.id)
            triples.append((rxn_index[rxn.id], j, kcat))

    n_r = len(catalysed_rxns)
    # 0 marks "no kcat assigned".
    kcat_arr = np.zeros(n_r, dtype=float)
    mat = sparse.lil_matrix((n_r, n_e), dtype=float)
    for rxn_idx, enz_idx, kcat in triples:
        # Multiple (rxn, enz) entries collapse to one (kcat per cell).
        kcat_arr[rxn_idx] = kcat
        mat[rxn_idx, enz_idx] = 1.0

    ec_model.ec = EcData(
        rxns=catalysed_rxns,
        kcat=kcat_arr,
        source=[""] * n_r,
        notes=[""] * n_r,
        eccodes=[""] * n_r,
        genes=list(enzyme_ids),  # gene == uniprot when unknown
        enzymes=list(enzyme_ids),
        mw=mw_arr,
        sequence=[""] * n_e,
        concs=concs_arr,
        rxn_enz_mat=mat.tocsr(),
    )
