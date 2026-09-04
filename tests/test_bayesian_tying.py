"""Tests for tying isozyme copies the kcat assignment could not tell apart."""
import numpy as np
import pytest

from geckopy.kcat_sensitivity_analysis.bayesian.tying import (
    apply_ties, base_reaction, isozyme_tie_map, n_free,
)


def test_base_reaction_strips_the_isozyme_suffix():
    assert base_reaction("r_0438_EXP_3") == "r_0438"
    assert base_reaction("r_0698") == "r_0698"


def test_copies_sharing_prior_and_source_are_tied():
    rxns = ["r_1_EXP_1", "r_1_EXP_2", "r_1_EXP_3", "r_2"]
    kcat0 = np.array([54.33, 54.33, 54.33, 7.0])
    sources = ["okp", "okp", "okp", "brenda"]

    tie = isozyme_tie_map(rxns, kcat0, sources)

    assert list(tie) == [0, 0, 0, 3]
    assert n_free(tie) == 2          # one for r_1's three copies, one for r_2


def test_copies_with_different_priors_stay_free():
    """A distinction the assignment made is a distinction worth keeping."""
    rxns = ["r_1_EXP_1", "r_1_EXP_2"]
    kcat0 = np.array([54.33, 12.0])

    tie = isozyme_tie_map(rxns, kcat0, ["okp", "okp"])

    assert list(tie) == [0, 1]
    assert n_free(tie) == 2


def test_a_shared_value_from_different_sources_is_a_coincidence():
    """Equal priors reached by different routes are not indistinguishable."""
    rxns = ["r_1_EXP_1", "r_1_EXP_2"]
    kcat0 = np.array([230.0, 230.0])

    assert list(isozyme_tie_map(rxns, kcat0, ["brenda", "custom"])) == [0, 1]
    # Without source information there is nothing to separate them by.
    assert list(isozyme_tie_map(rxns, kcat0)) == [0, 0]


def test_copies_of_different_reactions_are_never_tied():
    rxns = ["r_1_EXP_1", "r_2_EXP_1"]
    kcat0 = np.array([5.0, 5.0])

    assert list(isozyme_tie_map(rxns, kcat0, ["okp", "okp"])) == [0, 1]


def test_the_map_is_idempotent():
    """Representatives follow themselves, so applying twice changes nothing."""
    rxns = ["r_1_EXP_1", "r_1_EXP_2", "r_1_EXP_3"]
    tie = isozyme_tie_map(rxns, np.full(3, 2.0), ["okp"] * 3)

    assert list(tie[tie]) == list(tie)


def test_apply_ties_gives_copies_their_representatives_value():
    tie = np.array([0, 0, 0, 3])

    one = apply_ties(np.array([1.0, 9.0, 99.0, 5.0]), tie)
    assert list(one) == [1.0, 1.0, 1.0, 5.0]

    # A batch of particles is projected column by column.
    batch = np.array([[1.0, 2.0], [9.0, 8.0], [99.0, 88.0], [5.0, 6.0]])
    out = apply_ties(batch, tie)
    assert out.shape == (4, 2)
    assert list(out[:, 0]) == [1.0, 1.0, 1.0, 5.0]
    assert list(out[:, 1]) == [2.0, 2.0, 2.0, 6.0]


def test_apply_ties_rejects_a_mismatched_map():
    with pytest.raises(ValueError, match="tie_map has"):
        apply_ties(np.ones((3, 2)), np.array([0, 0]))


def test_tie_map_rejects_mismatched_inputs():
    with pytest.raises(ValueError, match="rxn_ids has"):
        isozyme_tie_map(["a"], np.ones(2))


def test_tying_reduces_the_dimension_a_search_sees():
    """A screen or optimiser working on free parameters sees fewer."""
    from geckopy.kcat_sensitivity_analysis.bayesian import tying

    rxns = ["r_1_EXP_1", "r_1_EXP_2", "r_2_EXP_1", "r_2_EXP_2", "r_3"]
    kcat0 = np.array([10.0, 10.0, 4.0, 9.0, 1.0])
    sources = ["brenda"] * 5
    tie = tying.isozyme_tie_map(rxns, kcat0, sources)

    # r_1's copies share a prior and tie; r_2's do not; r_3 stands alone.
    assert n_free(tie) == 4

    rng = np.random.default_rng(0)
    particles = rng.lognormal(size=(5, 50))
    tying.apply_ties(particles, tie)

    assert np.array_equal(particles[0], particles[1])
    assert not np.array_equal(particles[2], particles[3])
