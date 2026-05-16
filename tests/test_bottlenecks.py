"""Tests for get_enzyme_bottlenecks."""
from pathlib import Path

import cobra
import pytest

from geckopy import ModelAdapter, make_ec_model
from geckopy.ec_model.pipeline import apply_kcat_constraints
from geckopy.utilities import get_enzyme_bottlenecks

EXAMPLE_DIR = Path(__file__).parents[1] / "examples" / "ecTestGEM"


def _ectestgem_solvable():
    adapter = ModelAdapter.from_folder(EXAMPLE_DIR)
    cobra_model = cobra.io.read_sbml_model(str(adapter.params.conv_gem))
    ec_model = make_ec_model(cobra_model, adapter)
    ec_model.ec.kcat[:] = 10.0
    apply_kcat_constraints(ec_model)
    ec_model.objective = "R3"
    return ec_model


def test_returns_top_n():
    model = _ectestgem_solvable()
    df = get_enzyme_bottlenecks(model, top=3)
    assert len(df) == 3


def test_returns_at_most_top_when_fewer_enzymes():
    model = _ectestgem_solvable()
    df = get_enzyme_bottlenecks(model, top=999)
    assert len(df) == len(model.ec.enzymes)


def test_sorted_by_absolute_shadow_price():
    model = _ectestgem_solvable()
    df = get_enzyme_bottlenecks(model, top=10)
    abs_sp = df["shadow_price"].abs().to_list()
    assert abs_sp == sorted(abs_sp, reverse=True)


def test_columns():
    model = _ectestgem_solvable()
    df = get_enzyme_bottlenecks(model, top=3)
    assert list(df.columns) == [
        "gene", "shadow_price", "flux", "cap_usage", "upper_bound",
    ]
    assert df.index.name == "uniprot"


def test_raises_on_infeasible():
    model = _ectestgem_solvable()
    # Make infeasible: require biomass higher than achievable.
    r3 = model.reactions.get_by_id("R3")
    r3.lower_bound = 1e9
    with pytest.raises(RuntimeError, match="Solver status"):
        get_enzyme_bottlenecks(model, top=3)
