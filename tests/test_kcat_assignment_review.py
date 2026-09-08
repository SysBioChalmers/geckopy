"""Tests for flagging kcat assignments that deserve a second look."""
import numpy as np
import pytest

from geckopy.gather_kcats.review import (
    EcStats, Finding, findings_tsv, review_assignment,
)


def _stats(**kw):
    base = dict(n_kcat=5, kcat_max=100.0, kcat_median=50.0,
                values=frozenset({100.0, 50.0}), n_sa=0)
    base.update(kw)
    return EcStats(**base)


def test_no_database_evidence_for_the_ec_is_flagged():
    """The highest-yield check: the value did not come from a kcat."""
    rows = review_assignment(
        ["r_1"], np.array([0.0019]), ["brenda"], ["5.4.99.7"],
        ec_stats={"5.4.99.7": _stats(n_kcat=0, kcat_max=float("nan"),
                                     kcat_median=float("nan"),
                                     values=frozenset(), n_sa=4)},
    )
    assert [c for c in rows[0].checks if c == "no-ec-evidence"]
    assert "specific-activity" in rows[0].detail


def test_a_value_reused_across_reactions_is_a_fallback():
    ids = [f"r_{i}" for i in range(6)]
    kcats = np.full(6, 0.0918)
    rows = review_assignment(ids, kcats, ["brenda"] * 6, ["1.1.1.1"] * 6,
                             ec_stats={"1.1.1.1": _stats()}, repeat_min=5)
    assert all("repeated-value" in r.checks for r in rows)
    assert "6 reactions" in rows[0].detail

    # Below the threshold it is not remarkable.
    quiet = review_assignment(ids[:3], kcats[:3], ["brenda"] * 3,
                              ["1.1.1.1"] * 3,
                              ec_stats={"1.1.1.1": _stats()}, repeat_min=5)
    assert all("repeated-value" not in r.checks for r in quiet)


def test_taking_the_ec_maximum_is_flagged_only_when_the_spread_is_wide():
    st_wide = {"1.1.1.1": _stats(kcat_max=230.0, kcat_median=37.0, n_kcat=7)}
    wide = review_assignment(["r_1"], np.array([230.0]), ["brenda"],
                             ["1.1.1.1"], ec_stats=st_wide)
    assert "ec-maximum" in wide[0].checks
    assert "6x its median" in wide[0].detail

    # A maximum close to the median says nothing.
    st_tight = {"1.1.1.1": _stats(kcat_max=60.0, kcat_median=50.0)}
    tight = review_assignment(["r_1"], np.array([60.0]), ["brenda"],
                              ["1.1.1.1"], ec_stats=st_tight)
    assert tight == []


def test_custom_duplicating_its_own_ec_value_is_flagged():
    """Curation that copies the database adds no evidence but tightens
    the prior, since custom carries a narrower sigma than brenda."""
    st = {"4.2.1.11": _stats(kcat_max=230.0, kcat_median=37.0,
                             values=frozenset({230.0, 37.0}))}
    rows = review_assignment(["r_0366"], np.array([230.0]), ["custom"],
                             ["4.2.1.11"], ec_stats=st)
    assert "custom-duplicate" in rows[0].checks

    # A curated value the database does not hold is doing real work.
    rows = review_assignment(["r_0366"], np.array([12.0]), ["custom"],
                             ["4.2.1.11"], ec_stats=st)
    assert all("custom-duplicate" not in r.checks for r in rows)


def test_standard_fallbacks_and_tied_groups_are_quiet_by_default():
    """Both are deliberate, so listing every one is noise."""
    st = {"1.1.1.1": _stats()}
    assert review_assignment(["r_1"], np.array([13.0]), ["standard"],
                             ["1.1.1.1"], ec_stats=st) == []

    rows = review_assignment(["r_1"], np.array([13.0]), ["standard"],
                             ["1.1.1.1"], ec_stats=st,
                             include=("standard-fallback",))
    assert "standard-fallback" in rows[0].checks


def test_isozyme_copies_are_reported_once():
    ids = ["r_1_EXP_1", "r_1_EXP_2", "r_1_EXP_3"]
    kcats = np.full(3, 0.001)
    rows = review_assignment(ids, kcats, ["okp"] * 3, ["1.1.1.1"] * 3,
                             ec_stats={"1.1.1.1": _stats()})
    assert len(rows) == 1
    assert rows[0].n_isozymes == 3
    assert "3 isozyme copies" in rows[0].detail


def test_findings_rank_by_leverage_and_truncate_by_coverage():
    ids = [f"r_{i}" for i in range(4)]
    kcats = np.array([0.001, 0.002, 0.003, 0.004])   # all flagged: magnitude
    lev = np.array([0.1, 5.0, 3.0, 1.9])
    st = {"1.1.1.1": _stats()}

    rows = review_assignment(ids, kcats, ["brenda"] * 4, ["1.1.1.1"] * 4,
                             ec_stats=st, leverage=lev)
    assert [r.rxn_id for r in rows] == ["r_1", "r_2", "r_3", "r_0"]

    # 80% of 10.0 total is reached after 5.0 + 3.0.
    short = review_assignment(ids, kcats, ["brenda"] * 4, ["1.1.1.1"] * 4,
                              ec_stats=st, leverage=lev, coverage=0.8)
    assert [r.rxn_id for r in short] == ["r_1", "r_2"]

    assert len(review_assignment(ids, kcats, ["brenda"] * 4, ["1.1.1.1"] * 4,
                                 ec_stats=st, leverage=lev, top=1)) == 1


def test_shared_enzymes_are_counted_not_filtered():
    """A protein used elsewhere makes a correction a wider question."""
    ids = ["r_1", "r_2", "r_3"]
    rows = review_assignment(
        ids, np.array([0.001, 50.0, 50.0]), ["brenda"] * 3, ["1.1.1.1"] * 3,
        ec_stats={"1.1.1.1": _stats()},
        enzymes=[["P1", "P2"], ["P1"], ["P2"]],
    )
    assert len(rows) == 1                       # only r_1 trips a check
    assert rows[0].other_reactions == 2         # P1 in r_2, P2 in r_3


def test_checks_needing_a_database_are_skipped_without_one():
    rows = review_assignment(["r_1"], np.array([0.001]), ["brenda"],
                             ["1.1.1.1"])
    assert rows[0].checks == ("magnitude",)


def test_review_rejects_mismatched_inputs():
    with pytest.raises(ValueError, match="sources has"):
        review_assignment(["a", "b"], np.ones(2), ["brenda"], ["1.1.1.1"] * 2)
    with pytest.raises(ValueError, match="leverage has shape"):
        review_assignment(["a"], np.ones(1), ["brenda"], ["1.1.1.1"],
                          leverage=np.ones(2))


def test_tsv_carries_every_field():
    rows = [Finding("r_0698", "lanosterol synthase", "5.4.99.7", "brenda",
                    0.0019, 1.45, 1, 3, ("no-ec-evidence", "magnitude"),
                    "no kcat rows for EC 5.4.99.7")]
    text = findings_tsv(rows)
    header, body = text.splitlines()
    assert header.split("\t")[0] == "rxn_id"
    assert len(body.split("\t")) == len(header.split("\t"))
    assert "no-ec-evidence,magnitude" in body
