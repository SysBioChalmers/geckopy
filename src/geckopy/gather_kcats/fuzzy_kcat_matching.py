"""Fuzzy-match each ec reaction's EC + substrates against BRENDA kcat data.

Ported from GECKO MATLAB:
src/geckomat/gather_kcats/fuzzyKcatMatching.m.
"""
from __future__ import annotations

import logging
from typing import TYPE_CHECKING, Iterable, Optional

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    from ..databases import BrendaData, PhylDist
    from ..ec_model.ec_model import EcModel

logger = logging.getLogger(__name__)


# Diffusion-limited rate cap (Bar-Even et al. 2011), 1/s. Any kcat
# above this is clipped, matching MATLAB.
_DIFFUSION_LIMIT = 1e7

# Sentinel "no match" wildcard count; large so it never wins min-comparison.
_NO_MATCH_WC = 1000

# Output DataFrame schema. Used both in the main return and in the
# empty-result helper.
_RESULT_COLUMNS = [
    "rxn_id",
    "source",
    "eccode",
    "substrates",
    "genes",
    "kcat",
    "wildcard_level",
    "origin",
]


# Search order. Each entry is
#   (filter_organism, filter_substrate, use_sa_table, output_origin)
# The OUTPUT origin numbering follows the MATLAB docstring:
#   1: org+subs+kcat   2: any+subs+kcat   3: org+!subs+kcat
#   4: any+!subs+kcat  5: org+SA          6: any+SA
# but the SEARCH order tries org-SA (output 5) BEFORE any-no-subs-kcat
# (output 4), matching MATLAB exactly. See `D5` MATLAB-COMPAT note in
# docs/future_improvements.md.
_SEARCH_LEVELS: list[tuple[bool, bool, bool, int]] = [
    (True,  True,  False, 1),
    (False, True,  False, 2),
    (True,  False, False, 3),
    (True,  False, True,  5),
    (False, False, False, 4),
    (False, False, True,  6),
]


def fuzzy_kcat_matching(
    model: "EcModel",
    brenda: "BrendaData",
    phyl_dist: "PhylDist",
    *,
    ec_rxns: Optional[Iterable[str]] = None,
    force_wildcard_level: int = 0,
) -> pd.DataFrame:
    """Match each ec reaction to BRENDA kcat data with progressive
    relaxation of EC, substrate, and organism specificity.

    Ported from GECKO MATLAB:
    src/geckomat/gather_kcats/fuzzyKcatMatching.m.

    For each selected reaction, the EC string in
    ``model.ec.eccodes[i]`` is split on ``;`` into one or more EC
    tokens. Each token goes through:

    1. **Wildcard cascade** (per token): try the token as-is. If no
       match, replace its rightmost numeric level with ``-`` and try
       again. Repeat until matched or all four levels are wildcards.
    2. **Six-level matching** (per attempt): try increasingly
       relaxed criteria:

       * org + correct substrate + KCAT (output origin 1)
       * any organism + correct substrate + KCAT (output origin 2)
       * org + any substrate + KCAT (output origin 3)
       * org + specific activity (output origin 5)
       * any organism + any substrate + KCAT (output origin 4)
       * any organism + SA (output origin 6)

       When the EC has any wildcard, substrate-matching levels are
       skipped (the substrate field of an EC carries the catalytic
       specificity that wildcards already discard).

    Across the multiple EC tokens of one reaction, the best match is
    chosen by minimum wildcard count, then minimum origin (best output
    rank), then maximum kcat.

    The model's organism is read from ``model.adapter.params.org_name``.
    When BRENDA has no exact organism match, ``phyl_dist`` is used to
    pick the BRENDA row(s) for the phylogenetically closest
    organism(s) (with genus fallback for organisms missing from KEGG).

    MATLAB-COMPAT: The MATLAB SEARCH order tries org-SA BEFORE
    any-no-subs-kcat, but its OUTPUT ranking is the reverse. The
    consequence is that when both would match, org-SA wins (with output
    origin 5) even though any-no-subs-kcat (output origin 4) would have
    been ranked better by the docstring. geckopy replicates this
    behavior; tracked as a MATLAB-side bug in
    ``docs/future_improvements.md``.

    MATLAB-COMPAT: GECKO MATLAB takes a ``modelAdapter`` arg and reads
    the organism via ``adapter.params.org_name``. geckopy reads from
    ``model.adapter.params.org_name`` directly.

    MATLAB-COMPAT: GECKO MATLAB returns a struct of parallel arrays;
    geckopy returns a ``pandas.DataFrame`` for compatibility with the
    downstream ``apply_kcat_list`` and ``merge_dlkcat_and_fuzzy_kcats``
    functions, both of which are inherently relational operations.

    MATLAB-COMPAT: GECKO MATLAB tracks per-(origin, wildcard) match
    counts in an internal ``stats.matrix`` that is never returned.
    geckopy emits an aggregated ``logger.info`` summary instead.

    MATLAB-COMPAT: GECKO MATLAB's iterative EC escalation can produce
    invalid EC strings ("1.-.-.-.-") and crash with an index error if
    a token never matches at any wildcard level. geckopy caps escalation
    at 4 wildcards (full ``-.-.-.-``) and returns "no match" cleanly.

    MATLAB-COMPAT: GECKO MATLAB has dead code at
    ``if forceWClvl == 1`` (it is checked AFTER ``forceWClvl`` was
    decremented to 0 by the preceding ``while`` loop). geckopy
    interprets ``force_wildcard_level=N`` as "escalate every EC token
    by N wildcards from the right", which is the consistent extension
    of the loop's intent.

    Parameters
    ----------
    model
        EcModel with ``model.ec.rxns`` and ``model.ec.eccodes``
        already populated. Reactions whose ``ec.eccodes`` entry is
        empty contribute a row to the output with ``kcat=0`` and
        ``origin=NA``.
    brenda
        BRENDA tables loaded by ``load_brenda_data``.
    phyl_dist
        KEGG phylogenetic distance loaded by ``load_phyl_dist``.
    ec_rxns
        Optional iterable of reaction IDs (must appear in
        ``model.ec.rxns``) restricting which entries are matched.
        ``None`` means match all.
    force_wildcard_level
        Force at least this many wildcards on every EC token before
        the iterative escalation begins. Default 0 (no forcing).

    Returns
    -------
    pandas.DataFrame
        Columns:

        ============ ============================================================
        rxn_id       (str) reaction ID, matching ``model.ec.rxns``
        source       (str) ``"brenda"`` for every row
        eccode       (str) the ``;``-joined EC string used (after
                     ``force_wildcard_level`` but before iterative escalation)
        substrates   (list[str]) per-row substrate names
        genes        (list[str]) always empty for fuzzy matching; kept for
                     schema compatibility with future DLKcat output
        kcat         (float) selected kcat in 1/s, or 0.0 if no match
        wildcard_level (Int64) actual wildcard count at match time, or NA
        origin       (Int64) output origin 1..6, or NA if no match
        ============ ============================================================
    """
    if model.adapter is None:
        raise ValueError(
            "EcModel.adapter is None; fuzzy_kcat_matching needs an "
            "adapter for params.org_name."
        )
    if force_wildcard_level < 0:
        raise ValueError(
            f"force_wildcard_level must be >= 0, got {force_wildcard_level}"
        )

    organism_name = (model.adapter.params.org_name or "").strip()

    n = model.ec.n_rxns
    if n == 0:
        return _empty_result_df()

    if ec_rxns is None:
        positions: list[int] = list(range(n))
    else:
        ec_rxns_list = list(ec_rxns)
        index_by_id = {rid: i for i, rid in enumerate(model.ec.rxns)}
        unknown = [rid for rid in ec_rxns_list if rid not in index_by_id]
        if unknown:
            preview = unknown[:5]
            raise ValueError(
                f"{len(unknown)} reaction ID(s) in ec_rxns are not "
                f"present in model.ec.rxns (examples: {preview})"
            )
        positions = [index_by_id[rid] for rid in ec_rxns_list]

    if not positions:
        return _empty_result_df()

    kcat_ec_indices = _build_ec_indices(brenda.kcat)
    sa_ec_indices = _build_ec_indices(brenda.sa)
    org_index = _resolve_organism_index(organism_name, phyl_dist)

    out_rxn_ids: list[str] = []
    out_substrates: list[list[str]] = []
    out_kcats: list[float] = []
    out_eccodes: list[str] = []
    out_wcs: list[object] = []
    out_origins: list[object] = []

    for i in positions:
        rxn_id = model.ec.rxns[i]
        cobra_rxn_id = rxn_id[4:] if model.ec.gecko_light else rxn_id

        try:
            cobra_rxn = model.reactions.get_by_id(cobra_rxn_id)
            substrates, substrate_coeffs = _extract_substrates(cobra_rxn)
        except KeyError:
            substrates, substrate_coeffs = [], []

        out_rxn_ids.append(rxn_id)
        out_substrates.append(list(substrates))

        ec_string = model.ec.eccodes[i] or ""
        ec_tokens = [
            _apply_force_wildcards(t.strip(), force_wildcard_level)
            for t in ec_string.split(";") if t.strip()
        ]

        if not ec_tokens:
            out_kcats.append(0.0)
            out_eccodes.append("")
            out_wcs.append(pd.NA)
            out_origins.append(pd.NA)
            continue

        per_token: list[tuple[float, int, int]] = []
        for token in ec_tokens:
            kcat, origin, wc = _iterative_match_one_ec(
                token, substrates, substrate_coeffs,
                brenda.kcat, brenda.sa,
                kcat_ec_indices, sa_ec_indices,
                organism_name, phyl_dist, org_index,
            )
            per_token.append((kcat, origin, wc))

        matched = [r for r in per_token if r[1] > 0]
        out_eccodes.append(";".join(ec_tokens))
        if not matched:
            out_kcats.append(0.0)
            out_wcs.append(pd.NA)
            out_origins.append(pd.NA)
            continue

        min_wc = min(r[2] for r in matched)
        at_min_wc = [r for r in matched if r[2] == min_wc]
        min_origin = min(r[1] for r in at_min_wc)
        at_min_origin = [r for r in at_min_wc if r[1] == min_origin]
        best = max(at_min_origin, key=lambda r: r[0])

        out_kcats.append(best[0])
        out_wcs.append(best[2])
        out_origins.append(best[1])

    df = pd.DataFrame({
        "rxn_id": out_rxn_ids,
        "source": ["brenda"] * len(out_rxn_ids),
        "eccode": out_eccodes,
        "substrates": out_substrates,
        "genes": [[] for _ in out_rxn_ids],
        "kcat": out_kcats,
        "wildcard_level": pd.array(out_wcs, dtype="Int64"),
        "origin": pd.array(out_origins, dtype="Int64"),
    })

    matched_count = int(df["origin"].notna().sum())
    if matched_count > 0:
        origin_counts = df["origin"].value_counts(dropna=True).to_dict()
        wc_counts = df["wildcard_level"].value_counts(dropna=True).to_dict()
        logger.info(
            "fuzzy_kcat_matching: matched %d of %d reactions. "
            "Origin counts: %s. Wildcard-level counts: %s.",
            matched_count, len(df), dict(origin_counts), dict(wc_counts),
        )
    else:
        logger.info(
            "fuzzy_kcat_matching: no matches for any of %d reactions.",
            len(df),
        )

    return df


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #


def _empty_result_df() -> pd.DataFrame:
    """Return an empty DataFrame with the canonical schema."""
    return pd.DataFrame({
        "rxn_id": pd.Series([], dtype=str),
        "source": pd.Series([], dtype=str),
        "eccode": pd.Series([], dtype=str),
        "substrates": pd.Series([], dtype=object),
        "genes": pd.Series([], dtype=object),
        "kcat": pd.Series([], dtype=float),
        "wildcard_level": pd.array([], dtype="Int64"),
        "origin": pd.array([], dtype="Int64"),
    })


def _extract_substrates(rxn) -> tuple[list[str], list[float]]:
    """Return (names, abs_coeffs) for negative-coefficient metabolites.

    Uses ``metabolite.name`` (cobra equivalent of MATLAB ``metNames``).
    """
    names: list[str] = []
    coeffs: list[float] = []
    for met, coeff in rxn.metabolites.items():
        if coeff < 0:
            names.append(met.name)
            coeffs.append(-coeff)
    return names, coeffs


def _resolve_organism_index(
    org_name: str, phyl_dist: "PhylDist",
) -> Optional[int]:
    """Find the model organism's row in PhylDist. Returns None if no
    direct or genus-fallback match is found."""
    if not org_name:
        return None
    org_lower = org_name.lower()
    direct = phyl_dist.name_to_index.get(org_lower)
    if direct is not None:
        return direct
    parts = org_lower.split(None, 1)
    if parts:
        genus_indices = phyl_dist.genus_to_indices.get(parts[0], [])
        if genus_indices:
            return genus_indices[0]
    return None


def _escalate_wildcard(ec_token: str) -> Optional[str]:
    """Replace the rightmost numeric level with ``-``. Return None if
    every level is already ``-``."""
    parts = ec_token.split(".")
    if len(parts) != 4:
        return None
    for i in range(3, -1, -1):
        if parts[i] != "-":
            parts[i] = "-"
            return ".".join(parts)
    return None


def _apply_force_wildcards(ec_token: str, force_level: int) -> str:
    """Escalate ``ec_token`` by ``force_level`` wildcards from the right."""
    for _ in range(force_level):
        nxt = _escalate_wildcard(ec_token)
        if nxt is None:
            return ec_token
        ec_token = nxt
    return ec_token


def _build_ec_indices(table: pd.DataFrame) -> dict[str, np.ndarray]:
    """Group BRENDA table row indices by lowercased ``ec_code``."""
    if table.empty:
        return {}
    return {
        ec_lower: np.asarray(idxs, dtype=int)
        for ec_lower, idxs in (
            table.groupby(table["ec_code"].str.lower()).indices.items()
        )
    }


def _find_ec_rows(
    ec_token: str,
    ec_indices: dict[str, np.ndarray],
) -> np.ndarray:
    """Return BRENDA row indices matching ``ec_token``.

    Wildcard handling: the prefix up to the first ``-`` is matched
    case-insensitively against the start of each EC key. A fully-
    wildcarded ``-.-.-.-`` matches every row.
    """
    if not ec_indices:
        return np.array([], dtype=int)
    ec_lower = ec_token.lower()
    if "-" not in ec_token:
        return ec_indices.get(ec_lower, np.array([], dtype=int))
    prefix = ec_lower.split("-", 1)[0]
    if not prefix:
        all_idxs = list(ec_indices.values())
        return np.concatenate(all_idxs) if all_idxs else np.array([], dtype=int)
    matched = [arr for key, arr in ec_indices.items() if key.startswith(prefix)]
    return np.concatenate(matched) if matched else np.array([], dtype=int)


def _filter_by_substrate(
    table: pd.DataFrame,
    rows: np.ndarray,
    substrates: list[str],
) -> np.ndarray:
    if len(rows) == 0:
        return rows
    subs_lower = {s.lower() for s in substrates if s}
    if not subs_lower:
        return np.array([], dtype=int)
    sub_col = table["substrate"].iloc[rows].str.lower()
    keep = sub_col.isin(subs_lower).values
    return rows[keep]


def _filter_by_organism(
    table: pd.DataFrame,
    rows: np.ndarray,
    organism: Optional[str],
    phyl_dist: "PhylDist",
    org_index: Optional[int],
) -> np.ndarray:
    """Filter ``rows`` by organism.

    If ``organism`` is provided, exact case-insensitive match. If
    ``None``, use phylogenetic distance to keep only the rows whose
    organism is closest (in KEGG distance) to ``org_index``. If
    ``org_index`` is also None, no filtering is applied (matches MATLAB
    behaviour when the model organism has no KEGG entry).
    """
    if len(rows) == 0:
        return rows

    if organism is not None and organism != "":
        org_col = table["organism"].iloc[rows].str.lower()
        keep = org_col.values == organism.lower()
        return rows[keep]

    if org_index is None:
        return rows

    organisms = table["organism"].iloc[rows].values
    valid_kegg: list[int] = []
    valid_rows: list[int] = []
    for r, org in zip(rows, organisms):
        ol = str(org).lower()
        kegg_idx = phyl_dist.name_to_index.get(ol)
        if kegg_idx is None:
            parts = ol.split(None, 1)
            if parts:
                gen_indices = phyl_dist.genus_to_indices.get(parts[0], [])
                if gen_indices:
                    kegg_idx = gen_indices[0]
        if kegg_idx is not None:
            valid_kegg.append(kegg_idx)
            valid_rows.append(int(r))

    if not valid_rows:
        return np.array([], dtype=int)

    distances = phyl_dist.dist_matrix[org_index, valid_kegg]
    min_dist = float(distances.min())
    keep = [valid_rows[i] for i in range(len(valid_rows))
            if distances[i] == min_dist]
    return np.asarray(keep, dtype=int)


def _match_kcat(
    ec_token: str,
    substrates: list[str],
    substrate_coeffs: list[float],
    table: pd.DataFrame,
    ec_indices: dict[str, np.ndarray],
    organism: Optional[str],
    *,
    substrate_match: bool,
    sa: bool,
    phyl_dist: "PhylDist",
    org_index: Optional[int],
) -> tuple[float, int]:
    """One match attempt at one of the 6 search levels.

    Returns ``(max_kcat_in_1_per_s, num_matched_rows)``. ``num_matched
    == 0`` means no rows survived all filters.
    """
    rows = _find_ec_rows(ec_token, ec_indices)
    if len(rows) == 0:
        return 0.0, 0

    if substrate_match and not sa:
        rows = _filter_by_substrate(table, rows, substrates)
        if len(rows) == 0:
            return 0.0, 0

    rows = _filter_by_organism(table, rows, organism, phyl_dist, org_index)
    if len(rows) == 0:
        return 0.0, 0

    kcats = table["kcat"].iloc[rows].values

    if substrate_match and not sa:
        # MATLAB only counts kcat > 0 in this branch and divides by
        # min(substrCoeff). We match.
        positive = kcats > 0
        kcats = kcats[positive]
        if len(kcats) == 0:
            return 0.0, 0
        coeff = min(substrate_coeffs) if substrate_coeffs else 1.0
        if coeff > 0:
            kcats = kcats / coeff

    if len(kcats) == 0:
        return 0.0, 0

    n = len(kcats)
    kcat = float(np.max(kcats))
    if kcat > _DIFFUSION_LIMIT:
        kcat = _DIFFUSION_LIMIT
    return kcat, n


def _main_match(
    ec_token: str,
    substrates: list[str],
    substrate_coeffs: list[float],
    brenda_kcat: pd.DataFrame,
    brenda_sa: pd.DataFrame,
    kcat_ec_indices: dict[str, np.ndarray],
    sa_ec_indices: dict[str, np.ndarray],
    organism_name: str,
    phyl_dist: "PhylDist",
    org_index: Optional[int],
) -> tuple[float, int, int]:
    """Try each of the six search levels in MATLAB's order.

    Returns ``(kcat, output_origin, num_matches)`` from the first
    level that matches. ``output_origin == 0`` means no match.
    """
    has_wildcard = "-" in ec_token
    for filter_org, filter_subs, use_sa, output_origin in _SEARCH_LEVELS:
        if has_wildcard and filter_subs:
            continue
        table = brenda_sa if use_sa else brenda_kcat
        ec_idx = sa_ec_indices if use_sa else kcat_ec_indices
        org = organism_name if filter_org else None
        kcat, n = _match_kcat(
            ec_token, substrates, substrate_coeffs,
            table, ec_idx, org,
            substrate_match=filter_subs,
            sa=use_sa,
            phyl_dist=phyl_dist,
            org_index=org_index,
        )
        if n > 0:
            return kcat, output_origin, n
    return 0.0, 0, 0


def _iterative_match_one_ec(
    ec_token: str,
    substrates: list[str],
    substrate_coeffs: list[float],
    brenda_kcat: pd.DataFrame,
    brenda_sa: pd.DataFrame,
    kcat_ec_indices: dict[str, np.ndarray],
    sa_ec_indices: dict[str, np.ndarray],
    organism_name: str,
    phyl_dist: "PhylDist",
    org_index: Optional[int],
) -> tuple[float, int, int]:
    """Escalate wildcards in ``ec_token`` until a match is found or
    every level is wildcard.

    Returns ``(kcat, output_origin, wildcard_level)``. On no match,
    returns ``(0.0, 0, _NO_MATCH_WC)``.
    """
    current_ec = ec_token
    while True:
        wc = current_ec.count("-")
        kcat, origin, n = _main_match(
            current_ec, substrates, substrate_coeffs,
            brenda_kcat, brenda_sa,
            kcat_ec_indices, sa_ec_indices,
            organism_name, phyl_dist, org_index,
        )
        if origin > 0:
            return kcat, origin, wc
        nxt = _escalate_wildcard(current_ec)
        if nxt is None:
            return 0.0, 0, _NO_MATCH_WC
        current_ec = nxt
