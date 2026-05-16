"""Tests for find_ec_in_db."""
import logging

import numpy as np

from geckopy.get_enzyme_data import find_ec_in_db
from geckopy.get_enzyme_data.find_ec_in_db import (
    _subsumes,
    _wildcard_dedupe,
    _wildcard_intersection,
)


# --------------------------------------------------------------------------- #
# _subsumes
# --------------------------------------------------------------------------- #

def test_subsumes_equal_tokens_both_directions():
    assert _subsumes("1.1.1.1", "1.1.1.1")


def test_subsumes_full_wildcard_takes_anything():
    assert _subsumes("-.-.-.-", "1.1.1.1")
    assert _subsumes("-.-.-.-", "9.9.9.-")


def test_subsumes_partial_wildcard_takes_more_specific():
    assert _subsumes("1.1.1.-", "1.1.1.1")
    assert _subsumes("1.1.-.-", "1.1.5.5")
    assert _subsumes("1.-.-.-", "1.2.3.4")


def test_subsumes_specific_does_not_take_wildcard():
    assert not _subsumes("1.1.1.1", "1.1.1.-")


def test_subsumes_disjoint_returns_false():
    assert not _subsumes("1.1.1.1", "2.2.2.2")
    assert not _subsumes("1.1.1.1", "1.2.1.1")


def test_subsumes_multi_digit_levels_distinguished():
    """Crucial: `1.1.1.1` must NOT be claimed to subsume `1.1.1.10`."""
    assert not _subsumes("1.1.1.1", "1.1.1.10")
    assert not _subsumes("1.1.1.10", "1.1.1.1")


def test_subsumes_dash_at_intermediate_level():
    """A `-` at level N implies levels N..end are wildcard,
    regardless of what follows. So `1.-.5.5` and `1.-.-.-` are equivalent."""
    assert _subsumes("1.-.5.5", "1.7.5.5")
    assert _subsumes("1.-.5.5", "1.7.99.99")


def test_subsumes_rejects_wrong_arity():
    assert not _subsumes("1.1.1", "1.1.1.1")
    assert not _subsumes("1.1.1.1.1", "1.1.1.1")


# --------------------------------------------------------------------------- #
# _wildcard_dedupe
# --------------------------------------------------------------------------- #

def test_dedupe_empty():
    assert _wildcard_dedupe([]) == []


def test_dedupe_single():
    assert _wildcard_dedupe(["1.1.1.1"]) == ["1.1.1.1"]


def test_dedupe_unrelated_kept():
    assert _wildcard_dedupe(["1.1.1.1", "2.2.2.2"]) == ["1.1.1.1", "2.2.2.2"]


def test_dedupe_specific_subsumes_wildcard():
    assert _wildcard_dedupe(["1.1.1.1", "1.1.1.-"]) == ["1.1.1.1"]
    # Order independent
    assert _wildcard_dedupe(["1.1.1.-", "1.1.1.1"]) == ["1.1.1.1"]


def test_dedupe_equal_tokens_collapse_to_first():
    assert _wildcard_dedupe(["1.1.1.1", "1.1.1.1"]) == ["1.1.1.1"]


def test_dedupe_chain_of_wildcards():
    """`1.1.1.1` subsumes both `1.1.1.-` and `1.1.-.-`; only `1.1.1.1` remains."""
    assert _wildcard_dedupe(["1.1.-.-", "1.1.1.-", "1.1.1.1"]) == ["1.1.1.1"]


def test_dedupe_multi_digit_not_collapsed():
    assert _wildcard_dedupe(["1.1.1.1", "1.1.1.10"]) == ["1.1.1.1", "1.1.1.10"]


def test_dedupe_preserves_insertion_order_for_unrelated():
    assert _wildcard_dedupe(
        ["3.3.3.3", "1.1.1.1", "2.2.2.2"]
    ) == ["3.3.3.3", "1.1.1.1", "2.2.2.2"]


# --------------------------------------------------------------------------- #
# _wildcard_intersection
# --------------------------------------------------------------------------- #

def test_intersection_disjoint_yields_empty():
    assert _wildcard_intersection(["1.1.1.1"], ["2.2.2.2"]) == []


def test_intersection_equal_yields_one():
    assert _wildcard_intersection(["1.1.1.1"], ["1.1.1.1"]) == ["1.1.1.1"]


def test_intersection_picks_more_specific_of_a_pair():
    assert _wildcard_intersection(["1.1.1.-"], ["1.1.1.1"]) == ["1.1.1.1"]
    assert _wildcard_intersection(["1.1.1.1"], ["1.1.1.-"]) == ["1.1.1.1"]


def test_intersection_combines_multiple_pairs():
    """All four pairs are subsumption-pairs:
    (1.1.1.1, 1.1.1.-), (1.1.1.1, 2.2.2.2), (3.3.3.3, 1.1.1.-), (3.3.3.3, 2.2.2.2).
    The disjoint pairs contribute nothing; the subsumption pairs contribute
    the more specific one."""
    result = _wildcard_intersection(
        ["1.1.1.1", "3.3.3.3"], ["1.1.1.-", "2.2.2.2"]
    )
    assert result == ["1.1.1.1"]


# --------------------------------------------------------------------------- #
# Test fixture builder for find_ec_in_db
# --------------------------------------------------------------------------- #

def _make_db(rows: list[tuple[str, str, float]]) -> tuple[
    list[str], np.ndarray, dict[str, list[int]]
]:
    """Build (db_eccodes, db_mw, gene_to_protein_indices) from
    ``[(gene, ec_string, mw), ...]`` rows. Each row represents one
    protein in the DB; rows sharing a gene are grouped under that gene's
    indices in the order given.
    """
    db_eccodes = [ec for _, ec, _ in rows]
    db_mw = np.array([mw for _, _, mw in rows], dtype=float)
    gene_to_idx: dict[str, list[int]] = {}
    for i, (gene, _, _) in enumerate(rows):
        gene_to_idx.setdefault(gene, []).append(i)
    return db_eccodes, db_mw, gene_to_idx


# --------------------------------------------------------------------------- #
# Trivial cases
# --------------------------------------------------------------------------- #

def test_empty_gene_set_returns_empty():
    db_eccodes, db_mw, idx = _make_db([("g1", "1.1.1.1", 100.0)])
    assert find_ec_in_db([], db_eccodes, db_mw, idx) == ""


def test_single_gene_single_protein_single_ec():
    db_eccodes, db_mw, idx = _make_db([("g1", "1.1.1.1", 100.0)])
    assert find_ec_in_db(["g1"], db_eccodes, db_mw, idx) == "1.1.1.1"


def test_gene_not_in_db_silently_skipped():
    db_eccodes, db_mw, idx = _make_db([("g1", "1.1.1.1", 100.0)])
    # Only g_unknown requested; not in idx -> no EC, return empty.
    assert find_ec_in_db(["g_unknown"], db_eccodes, db_mw, idx) == ""


def test_gene_in_db_but_protein_has_empty_ec_returns_empty():
    db_eccodes, db_mw, idx = _make_db([("g1", "", 100.0)])
    assert find_ec_in_db(["g1"], db_eccodes, db_mw, idx) == ""


def test_mix_of_known_and_unknown_genes():
    db_eccodes, db_mw, idx = _make_db([("g1", "1.1.1.1", 100.0)])
    # g_unknown contributes nothing; g1 contributes its EC.
    assert find_ec_in_db(
        ["g_unknown", "g1"], db_eccodes, db_mw, idx
    ) == "1.1.1.1"


# --------------------------------------------------------------------------- #
# MW tiebreak (single distinct EC across multiple proteins)
# --------------------------------------------------------------------------- #

def test_mw_tiebreak_picks_lightest_protein():
    """When all proteins for a gene share an EC, the lightest-MW one is
    chosen. Since the EC string is the same for all, this only matters
    structurally; the visible result is the EC string."""
    db_eccodes, db_mw, idx = _make_db([
        ("g1", "1.1.1.1", 200.0),
        ("g1", "1.1.1.1", 100.0),  # lightest
        ("g1", "1.1.1.1", 300.0),
    ])
    assert find_ec_in_db(["g1"], db_eccodes, db_mw, idx) == "1.1.1.1"


def test_mw_tiebreak_handles_nan():
    """NaN MW must not crash argmin; it is treated as +inf."""
    db_eccodes, db_mw, idx = _make_db([
        ("g1", "1.1.1.1", float("nan")),
        ("g1", "1.1.1.1", 100.0),  # lightest finite
    ])
    assert find_ec_in_db(["g1"], db_eccodes, db_mw, idx) == "1.1.1.1"


def test_mw_tiebreak_all_nan_falls_back_to_first():
    db_eccodes, db_mw, idx = _make_db([
        ("g1", "1.1.1.1", float("nan")),
        ("g1", "1.1.1.1", float("nan")),
    ])
    assert find_ec_in_db(["g1"], db_eccodes, db_mw, idx) == "1.1.1.1"


# --------------------------------------------------------------------------- #
# Conflict (multiple distinct ECs for one gene)
# --------------------------------------------------------------------------- #

def test_conflict_logs_and_picks_first(caplog):
    db_eccodes, db_mw, idx = _make_db([
        ("g1", "1.1.1.1", 100.0),
        ("g1", "2.2.2.2", 50.0),  # lighter, but a different EC
    ])
    with caplog.at_level(logging.WARNING):
        result = find_ec_in_db(["g1"], db_eccodes, db_mw, idx)
    assert result == "1.1.1.1"
    assert "g1" in caplog.text
    assert "1.1.1.1" in caplog.text
    assert "2.2.2.2" in caplog.text


def test_conflict_does_not_apply_mw_tiebreak():
    """When ECs disagree, MW is irrelevant; first DB-order EC wins."""
    db_eccodes, db_mw, idx = _make_db([
        ("g1", "1.1.1.1", 999.0),
        ("g1", "2.2.2.2", 1.0),  # much lighter, but second
    ])
    assert find_ec_in_db(["g1"], db_eccodes, db_mw, idx) == "1.1.1.1"


def test_no_conflict_warning_when_all_proteins_agree(caplog):
    db_eccodes, db_mw, idx = _make_db([
        ("g1", "1.1.1.1", 100.0),
        ("g1", "1.1.1.1", 200.0),
    ])
    with caplog.at_level(logging.WARNING):
        find_ec_in_db(["g1"], db_eccodes, db_mw, idx)
    assert "distinct" not in caplog.text


# --------------------------------------------------------------------------- #
# Multi-EC string per protein
# --------------------------------------------------------------------------- #

def test_multi_ec_protein_split_into_tokens():
    db_eccodes, db_mw, idx = _make_db([
        ("g1", "1.1.1.1;2.2.2.2", 100.0),
    ])
    assert find_ec_in_db(["g1"], db_eccodes, db_mw, idx) == "1.1.1.1;2.2.2.2"


def test_multi_ec_protein_with_wildcard_dedupe():
    """Within a single protein's EC string, wildcard dedupe applies."""
    db_eccodes, db_mw, idx = _make_db([
        ("g1", "1.1.1.1;1.1.1.-", 100.0),
    ])
    assert find_ec_in_db(["g1"], db_eccodes, db_mw, idx) == "1.1.1.1"


# --------------------------------------------------------------------------- #
# Cross-gene reconciliation
# --------------------------------------------------------------------------- #

def test_two_genes_same_ec_intersection_kept():
    db_eccodes, db_mw, idx = _make_db([
        ("g1", "1.1.1.1", 100.0),
        ("g2", "1.1.1.1", 100.0),
    ])
    assert find_ec_in_db(["g1", "g2"], db_eccodes, db_mw, idx) == "1.1.1.1"


def test_two_genes_disjoint_ecs_falls_back_to_union():
    db_eccodes, db_mw, idx = _make_db([
        ("g1", "1.1.1.1", 100.0),
        ("g2", "2.2.2.2", 100.0),
    ])
    # No intersection -> union (both kept, in input order).
    result = find_ec_in_db(["g1", "g2"], db_eccodes, db_mw, idx)
    assert set(result.split(";")) == {"1.1.1.1", "2.2.2.2"}


def test_two_genes_one_wildcard_intersection_picks_specific():
    """Intersection with wildcard subsumption: `1.1.1.1` and `1.1.1.-`
    intersect at `1.1.1.1`."""
    db_eccodes, db_mw, idx = _make_db([
        ("g1", "1.1.1.1", 100.0),
        ("g2", "1.1.1.-", 100.0),
    ])
    assert find_ec_in_db(["g1", "g2"], db_eccodes, db_mw, idx) == "1.1.1.1"


def test_three_genes_partial_intersection():
    """g1, g2 share `1.1.1.1`; g3 has only `2.2.2.2`. No 3-way
    intersection -> falls back to union (deduped)."""
    db_eccodes, db_mw, idx = _make_db([
        ("g1", "1.1.1.1", 100.0),
        ("g2", "1.1.1.1", 100.0),
        ("g3", "2.2.2.2", 100.0),
    ])
    result = find_ec_in_db(["g1", "g2", "g3"], db_eccodes, db_mw, idx)
    # Union is the survivor; 1.1.1.1 (deduped from g1 and g2) plus 2.2.2.2.
    assert set(result.split(";")) == {"1.1.1.1", "2.2.2.2"}


def test_one_gene_missing_does_not_break_intersection():
    """If a gene contributes nothing, the cross-gene loop just skips it."""
    db_eccodes, db_mw, idx = _make_db([
        ("g1", "1.1.1.1", 100.0),
        ("g3", "1.1.1.1", 100.0),
    ])
    # g2 not in DB; result should still be 1.1.1.1 from g1 + g3.
    assert find_ec_in_db(
        ["g1", "g2", "g3"], db_eccodes, db_mw, idx
    ) == "1.1.1.1"


def test_all_genes_missing_returns_empty():
    db_eccodes, db_mw, idx = _make_db([("g1", "1.1.1.1", 100.0)])
    assert find_ec_in_db(["x", "y", "z"], db_eccodes, db_mw, idx) == ""


# --------------------------------------------------------------------------- #
# Final dedupe (geckopy divergence from MATLAB)
# --------------------------------------------------------------------------- #

def test_intersection_result_is_deduped():
    """Two genes both with EC `1.1.1.1` should yield a SINGLE
    `1.1.1.1` in the result, not `1.1.1.1;1.1.1.1`. MATLAB would
    produce duplicates here; geckopy dedupes."""
    db_eccodes, db_mw, idx = _make_db([
        ("g1", "1.1.1.1", 100.0),
        ("g2", "1.1.1.1", 100.0),
    ])
    result = find_ec_in_db(["g1", "g2"], db_eccodes, db_mw, idx)
    assert result == "1.1.1.1"
    # Sanity: no duplicate token.
    assert result.count("1.1.1.1") == 1


# --------------------------------------------------------------------------- #
# Optional `conflicts` accumulator
# --------------------------------------------------------------------------- #

def test_conflicts_accumulator_is_optional():
    """Default behaviour (no accumulator) is unchanged."""
    db_eccodes, db_mw, idx = _make_db([
        ("g1", "1.1.1.1", 100.0),
        ("g1", "2.2.2.2", 200.0),
    ])
    result = find_ec_in_db(["g1"], db_eccodes, db_mw, idx)
    assert result == "1.1.1.1"


def test_conflicts_accumulator_stays_empty_when_no_conflict():
    db_eccodes, db_mw, idx = _make_db([
        ("g1", "1.1.1.1", 100.0),
        ("g1", "1.1.1.1", 200.0),  # same EC, no conflict
    ])
    conflicts: list = []
    find_ec_in_db(["g1"], db_eccodes, db_mw, idx, conflicts=conflicts)
    assert conflicts == []


def test_conflicts_accumulator_collects_one_per_conflicting_gene():
    db_eccodes, db_mw, idx = _make_db([
        ("g1", "1.1.1.1", 100.0),
        ("g1", "2.2.2.2", 50.0),
        ("g2", "3.3.3.3", 100.0),  # no conflict
        ("g3", "4.4.4.4", 100.0),
        ("g3", "5.5.5.5", 100.0),
    ])
    conflicts: list = []
    find_ec_in_db(
        ["g1", "g2", "g3"], db_eccodes, db_mw, idx, conflicts=conflicts,
    )
    assert len(conflicts) == 2
    genes = [c[0] for c in conflicts]
    assert genes == ["g1", "g3"]


def test_conflicts_accumulator_protein_indices_are_first_per_distinct_ec():
    """Protein indices appended should be the first DB row seen for
    each distinct EC, mirroring MATLAB ``unique('stable')``."""
    db_eccodes, db_mw, idx = _make_db([
        ("g1", "1.1.1.1", 100.0),  # idx 0, EC A first match
        ("g1", "1.1.1.1", 200.0),  # idx 1, EC A duplicate
        ("g1", "2.2.2.2", 300.0),  # idx 2, EC B first match
        ("g1", "2.2.2.2", 400.0),  # idx 3, EC B duplicate
    ])
    conflicts: list = []
    find_ec_in_db(["g1"], db_eccodes, db_mw, idx, conflicts=conflicts)
    assert len(conflicts) == 1
    gene, protein_indices = conflicts[0]
    assert gene == "g1"
    assert protein_indices == [0, 2]


def test_conflicts_accumulator_appends_across_calls():
    db_eccodes, db_mw, idx = _make_db([
        ("g1", "1.1.1.1", 100.0),
        ("g1", "2.2.2.2", 100.0),
        ("g2", "3.3.3.3", 100.0),
        ("g2", "4.4.4.4", 100.0),
    ])
    conflicts: list = []
    find_ec_in_db(["g1"], db_eccodes, db_mw, idx, conflicts=conflicts)
    find_ec_in_db(["g2"], db_eccodes, db_mw, idx, conflicts=conflicts)
    assert [c[0] for c in conflicts] == ["g1", "g2"]


# --------------------------------------------------------------------------- #
# Complex realistic case
# --------------------------------------------------------------------------- #

def test_complex_three_subunit_with_isozyme_and_wildcard():
    """Three-subunit complex: g1 has two proteins sharing one EC,
    g2 has one protein with a more specific version of g1's EC,
    g3 has only a wildcard EC matching the family."""
    db_eccodes, db_mw, idx = _make_db([
        ("g1", "1.1.1.1", 100.0),
        ("g1", "1.1.1.1", 200.0),
        ("g2", "1.1.1.1", 50.0),
        ("g3", "1.1.-.-", 80.0),
    ])
    # Intersection finds 1.1.1.1 across g1, g2 (both have it) and g3
    # (1.1.-.- subsumes 1.1.1.1 -> intersection picks 1.1.1.1).
    assert find_ec_in_db(
        ["g1", "g2", "g3"], db_eccodes, db_mw, idx
    ) == "1.1.1.1"
