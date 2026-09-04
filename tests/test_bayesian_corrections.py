"""Tests for reporting a tuning result as its corrections."""
import numpy as np
import pytest

from geckopy.kcat_sensitivity_analysis.bayesian.corrections import (
    Correction, annotate_from_model, corrections, corrections_tsv,
)


def test_corrections_lists_only_changed_kcats():
    kcat0 = np.array([1.0, 1.0, 1.0, 1.0])
    kcat = np.array([1.0, 1.01, 4.0, 0.5])   # unchanged, within tol, up, down
    rows = corrections(kcat, kcat0, ["a", "b", "c", "d"])

    assert [r.rxn_id for r in rows] == ["c", "d"]
    # Fold change is direction-free: quartered and quadrupled both read 4.
    assert rows[0].fold_change == pytest.approx(4.0)
    assert rows[1].fold_change == pytest.approx(2.0)
    assert rows[0].kcat_prior == 1.0 and rows[0].kcat_tuned == 4.0


def test_corrections_rank_by_leverage_not_by_movement():
    """Sorting by how far a kcat moved puts the unidentifiable first."""
    kcat0 = np.ones(3)
    kcat = np.array([100.0, 3.0, 2.0])
    leverage = np.array([0.001, 1.0, 0.5])    # the big mover barely matters

    by_leverage = corrections(kcat, kcat0, ["drifter", "big", "mid"],
                              leverage=leverage)
    assert [r.rxn_id for r in by_leverage] == ["big", "mid", "drifter"]

    # Without leverage there is nothing to rank by but movement.
    by_movement = corrections(kcat, kcat0, ["drifter", "big", "mid"])
    assert [r.rxn_id for r in by_movement] == ["drifter", "big", "mid"]
    assert all(r.leverage_share == 0.0 for r in by_movement)


def test_cumulative_share_shows_where_the_list_stops_mattering():
    kcat0 = np.ones(4)
    kcat = np.array([2.0, 2.0, 2.0, 1.0])     # three changed, one not
    leverage = np.array([0.6, 0.3, 0.1, 5.0])  # the unchanged one dominates

    rows = corrections(kcat, kcat0, ["a", "b", "c", "untouched"],
                       leverage=leverage)
    # Shares are of *total* leverage, so changing everything that matters
    # is distinguishable from changing everything.
    assert rows[0].leverage_share == pytest.approx(0.1)
    assert rows[-1].cumulative_share == pytest.approx(0.16666, abs=1e-4)
    assert [r.cumulative_share for r in rows] == sorted(
        r.cumulative_share for r in rows)


def test_missing_provenance_is_reported_not_dropped():
    kcat0 = np.ones(2)
    kcat = np.array([2.0, 3.0])
    rows = corrections(kcat, kcat0, ["a", "b"], sources=["brenda", "custom"])

    assert [r.source for r in rows] == ["custom", "brenda"]  # ordered by fold
    assert all(r.name == "" and r.ec_code == "" for r in rows)


def test_corrections_reject_mismatched_inputs():
    with pytest.raises(ValueError, match="rxn_ids has"):
        corrections(np.ones(3), np.ones(3), ["a"])
    with pytest.raises(ValueError, match="kcat0 has"):
        corrections(np.ones(3), np.ones(2), ["a", "b", "c"])
    with pytest.raises(ValueError, match="leverage has shape"):
        corrections(np.ones(3), np.ones(3), ["a", "b", "c"],
                    leverage=np.ones(2))


def test_tsv_round_trips_every_field():
    rows = corrections(np.array([4.0]), np.array([1.0]), ["r_0698"],
                       sources=["brenda"], names=["lanosterol synthase"],
                       ec_codes=["5.4.99.7"], leverage=np.array([2.0]))
    text = corrections_tsv(rows)
    header, body = text.splitlines()

    assert header.split("\t")[0] == "rxn_id"
    assert len(body.split("\t")) == len(header.split("\t"))
    assert "lanosterol synthase" in body and "5.4.99.7" in body
    assert text.endswith("\n")


def test_tsv_of_no_corrections_is_a_header_alone():
    rows = corrections(np.ones(2), np.ones(2), ["a", "b"])
    assert rows == []
    assert corrections_tsv(rows).strip().startswith("rxn_id")


class _FakeEc:
    rxns = ["r_1", "r_2_EXP_1", "r_2_EXP_2"]
    eccodes = ["1.1.1.1", "2.2.2.2", "2.2.2.2"]
    source = ["brenda", "custom", "custom"]


class _FakeReaction:
    def __init__(self, name):
        self.name = name


class _FakeReactions:
    def __init__(self, mapping):
        self._m = mapping

    def get_by_id(self, rxn_id):
        return _FakeReaction(self._m[rxn_id])


class _FakeModel:
    ec = _FakeEc()
    reactions = _FakeReactions({"r_1": "alcohol dehydrogenase",
                                "r_2": "lanosterol synthase"})


def test_annotations_strip_the_isozyme_suffix_to_find_the_name():
    got = annotate_from_model(_FakeModel(), ["r_2_EXP_2", "r_1"])

    # The name lives on the reaction, which has no _EXP_ suffix.
    assert got["names"] == ["lanosterol synthase", "alcohol dehydrogenase"]
    assert got["ec_codes"] == ["2.2.2.2", "1.1.1.1"]
    assert got["sources"] == ["custom", "brenda"]


def test_annotations_report_unknown_ids_as_blank():
    got = annotate_from_model(_FakeModel(), ["r_absent"])
    assert got == {"names": [""], "ec_codes": [""], "sources": [""]}
