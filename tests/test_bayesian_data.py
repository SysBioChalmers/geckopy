"""Tests for kcat_sensitivity_analysis.bayesian.data."""
from pathlib import Path

from geckopy import ModelAdapter
from geckopy.kcat_sensitivity_analysis.bayesian import (
    BayesianData,
    load_bayesian_data,
)


def _adapter(tmp_path: Path) -> ModelAdapter:
    (tmp_path / "model_adapter.toml").write_text(
        'conv_gem = "dummy.xml"\n'
        'org_name = "test"\n'
    )
    return ModelAdapter.from_folder(tmp_path)


_FLUX_DATA_TSV = (
    "Condition\tPtot\tgrRate\tglucose (r_1714)\tethanol (r_1761)\t"
    "bayesianRMSEweight\tsource\n"
    "glucose\tNaN\t0.36\t-13.33\t19.83\t1\tDLKcat\n"
    "ethanol\tNaN\t0.12\tNaN\t-8.5\t1\tbrenda\n"
)

_MAX_GROWTH_TSV = (
    "Condition\tPtot\tgrRate\tglucose (r_1714)\tfructose (r_1709)\t"
    "bayesianRMSEweight\tsource\n"
    "glucose\tNaN\t0.41\t-1000\tNaN\t1\tDLKcat\n"
    "fructose\tNaN\t0.338\tNaN\t-1000\t1\tDLKcat\n"
)

_ZERO_EXCH_TSV = "Rxns\nr_1542\nr_1545\nr_1546\n"


def test_load_bayesian_data_all_files_present(tmp_path):
    adapter = _adapter(tmp_path)
    data_dir = tmp_path / "data"
    data_dir.mkdir()
    (data_dir / "bayesianFluxData.tsv").write_text(_FLUX_DATA_TSV)
    (data_dir / "bayesianMaxGrowth.tsv").write_text(_MAX_GROWTH_TSV)
    (data_dir / "bayesianZeroExch.tsv").write_text(_ZERO_EXCH_TSV)

    bay_data = load_bayesian_data(adapter)

    assert isinstance(bay_data, BayesianData)

    assert bay_data.flux_data is not None
    assert bay_data.flux_data.conds == ["glucose", "ethanol"]
    assert bay_data.flux_data.exch_rxn_ids == ["r_1714", "r_1761"]
    assert bay_data.flux_data.source == ["DLKcat", "brenda"]
    assert list(bay_data.flux_data.bayesian_rmse_weight) == [1.0, 1.0]

    assert bay_data.max_grate is not None
    assert bay_data.max_grate.conds == ["glucose", "fructose"]
    assert bay_data.max_grate.exch_rxn_ids == ["r_1714", "r_1709"]

    assert bay_data.zero_flux == ["r_1542", "r_1545", "r_1546"]


def test_load_bayesian_data_all_files_missing(tmp_path):
    adapter = _adapter(tmp_path)

    bay_data = load_bayesian_data(adapter)

    assert bay_data.flux_data is None
    assert bay_data.max_grate is None
    assert bay_data.zero_flux == []


def test_load_bayesian_data_partial(tmp_path):
    """Only flux_data present: max_grate stays None (not a bogus
    partial struct), matching MATLAB's *intended* behaviour for a
    missing optional file -- unlike the asymmetric-guard bug flagged
    in the MATLAB REVIEW.md, both fields are treated identically."""
    adapter = _adapter(tmp_path)
    data_dir = tmp_path / "data"
    data_dir.mkdir()
    (data_dir / "bayesianFluxData.tsv").write_text(_FLUX_DATA_TSV)

    bay_data = load_bayesian_data(adapter)

    assert bay_data.flux_data is not None
    assert bay_data.max_grate is None
    assert bay_data.zero_flux == []
