"""Unit tests for the BRENDA bulk-JSON parser."""
from __future__ import annotations

from pathlib import Path

import pytest

from geckopy.databases.brenda import Row, parse_brenda_json

FIXTURE = Path(__file__).parent / "data" / "brenda_minimal.json"


@pytest.fixture(scope="module")
def rows() -> list[Row]:
    return list(parse_brenda_json(FIXTURE))


def test_emits_rows(rows):
    assert rows, "parser yielded no rows from the minimal fixture"


def test_skips_spontaneous(rows):
    assert all(r.ec != "spontaneous" for r in rows)


def test_kcat_normal_row(rows):
    matches = [
        r for r in rows
        if r.kind == "kcat"
        and r.ec == "1.1.1.1"
        and r.substrate == "ethanol"
        and r.organism == "saccharomyces cerevisiae"
        and r.value == 23.5
    ]
    assert len(matches) == 1
    assert matches[0].references == ("PMID:11111", "PMID:22222")


def test_kcat_fans_out_across_proteins(rows):
    organisms = {
        r.organism for r in rows
        if r.kind == "kcat" and r.ec == "1.1.1.1" and r.value == 23.5
    }
    assert organisms == {"saccharomyces cerevisiae", "escherichia coli"}


def test_kcat_range_takes_upper_bound(rows):
    matches = [
        r for r in rows
        if r.kind == "kcat" and r.ec == "1.1.1.1" and r.substrate == "nadh"
    ]
    assert len(matches) == 1
    assert matches[0].value == 2.5


def test_kcat_above_physical_limit_filtered(rows):
    assert not any(r.kind == "kcat" and r.value >= 5e7 for r in rows)


def test_missing_value_sentinel_filtered(rows):
    assert not any(r.value == -999.0 for r in rows)


def test_mutant_rows_filtered(rows):
    for r in rows:
        assert r.value != 99.0, "site-directed mutant should be filtered"
        assert r.value != 88.0, "'mutated' comment should be filtered"


def test_unknown_protein_id_skipped(rows):
    assert not any(
        r.kind == "kcat" and r.ec == "1.1.1.1" and r.value == 11.0
        for r in rows
    ), "protein id 99 not in protein map; row must be skipped"


def test_empty_organism_skipped(rows):
    assert not any(
        r.kind == "kcat" and r.ec == "1.1.1.1" and r.value == 7.0
        for r in rows
    ), "protein with empty organism string must be skipped"


def test_sa_emits_with_star_substrate(rows):
    matches = [r for r in rows if r.kind == "sa" and r.ec == "1.1.1.1"]
    assert len(matches) == 1
    assert matches[0].substrate == "*"
    assert matches[0].value == 598.0
    assert matches[0].organism == "saccharomyces cerevisiae"
    assert matches[0].references == ("PMID:22222",)


def test_mw_emits_with_star_substrate(rows):
    matches = [r for r in rows if r.kind == "mw" and r.ec == "1.1.1.1"]
    assert len(matches) == 1
    assert matches[0].substrate == "*"
    assert matches[0].value == 130000.0


def test_empty_proteins_list_skipped(rows):
    assert not any(
        r.kind == "kcat" and r.ec == "1.2.3.4" and r.substrate == "methanol"
        for r in rows
    )


def test_garbage_value_skipped(rows):
    assert not any(r.kind == "kcat" and r.ec == "1.2.3.4" and r.value == 0.0 for r in rows)


def test_row_with_no_references_emits_empty_tuple(rows):
    matches = [
        r for r in rows
        if r.kind == "kcat" and r.ec == "1.2.3.4" and r.value == 50.0
    ]
    assert len(matches) == 1
    assert matches[0].references == ()


def test_references_without_pmid_dropped(rows):
    matches = [
        r for r in rows
        if r.kind == "kcat" and r.ec == "1.2.3.4" and r.value == 60.0
    ]
    assert len(matches) == 1
    assert matches[0].references == ()
