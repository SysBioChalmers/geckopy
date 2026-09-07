"""Tests for kcat_tuning.evolutionary_tuning.data."""
from pathlib import Path

from geckopy import ModelAdapter
from geckopy.kcat_tuning.evolutionary_tuning import (
    TuningData,
    load_tuning_data,
)


def _adapter(tmp_path: Path) -> ModelAdapter:
    (tmp_path / "model_adapter.toml").write_text(
        'conv_gem = "dummy.xml"\n'
        'org_name = "test"\n'
    )
    return ModelAdapter.from_folder(tmp_path)


_FLUX_DATA_TSV = (
    "Condition\tPtot\tgrRate\tglucose (r_1714)\tethanol (r_1761)\t"
    "tuningRMSEweight\tsource\n"
    "glucose\tNaN\t0.36\t-13.33\t19.83\t1\tDLKcat\n"
    "ethanol\tNaN\t0.12\tNaN\t-8.5\t1\tbrenda\n"
)

_MAX_GROWTH_TSV = (
    "Condition\tPtot\tgrRate\tglucose (r_1714)\tfructose (r_1709)\t"
    "tuningRMSEweight\tsource\n"
    "glucose\tNaN\t0.41\t-1000\tNaN\t1\tDLKcat\n"
    "fructose\tNaN\t0.338\tNaN\t-1000\t1\tDLKcat\n"
)

_ZERO_EXCH_TSV = "Rxns\nr_1542\nr_1545\nr_1546\n"


def test_load_tuning_data_all_files_present(tmp_path):
    adapter = _adapter(tmp_path)
    data_dir = tmp_path / "data"
    data_dir.mkdir()
    (data_dir / "tuningFluxData.tsv").write_text(_FLUX_DATA_TSV)
    (data_dir / "tuningMaxGrowth.tsv").write_text(_MAX_GROWTH_TSV)
    (data_dir / "tuningZeroExch.tsv").write_text(_ZERO_EXCH_TSV)

    tuning_data = load_tuning_data(adapter)

    assert isinstance(tuning_data, TuningData)

    assert tuning_data.flux_data is not None
    assert tuning_data.flux_data.conds == ["glucose", "ethanol"]
    assert tuning_data.flux_data.exch_rxn_ids == ["r_1714", "r_1761"]
    assert tuning_data.flux_data.source == ["DLKcat", "brenda"]
    assert list(tuning_data.flux_data.tuning_rmse_weight) == [1.0, 1.0]

    assert tuning_data.max_grate is not None
    assert tuning_data.max_grate.conds == ["glucose", "fructose"]
    assert tuning_data.max_grate.exch_rxn_ids == ["r_1714", "r_1709"]

    assert tuning_data.zero_flux == ["r_1542", "r_1545", "r_1546"]


def test_load_tuning_data_all_files_missing(tmp_path):
    adapter = _adapter(tmp_path)

    tuning_data = load_tuning_data(adapter)

    assert tuning_data.flux_data is None
    assert tuning_data.max_grate is None
    assert tuning_data.zero_flux == []


def test_load_tuning_data_partial(tmp_path):
    """Only flux_data present: max_grate stays None (not a bogus
    partial struct), matching MATLAB's *intended* behaviour for a
    missing optional file -- unlike the asymmetric-guard bug flagged
    in the MATLAB REVIEW.md, both fields are treated identically."""
    adapter = _adapter(tmp_path)
    data_dir = tmp_path / "data"
    data_dir.mkdir()
    (data_dir / "tuningFluxData.tsv").write_text(_FLUX_DATA_TSV)

    tuning_data = load_tuning_data(adapter)

    assert tuning_data.flux_data is not None
    assert tuning_data.max_grate is None
    assert tuning_data.zero_flux == []
