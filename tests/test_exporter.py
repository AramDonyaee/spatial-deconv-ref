"""Tests for cell2location and RCTD export formatting."""
import anndata as ad
import pandas as pd
from screfsuggest.exporter import export_for_cell2location, export_for_rctd


def test_export_cell2location(mock_single_cell_adata, tmp_path):
    out_file = tmp_path / "c2l_ref.h5ad"
    export_for_cell2location(mock_single_cell_adata, out_file)

    assert out_file.exists()
    loaded = ad.read_h5ad(out_file)
    assert "cell_type" in loaded.obs
    assert loaded.n_obs == mock_single_cell_adata.n_obs


def test_export_rctd(mock_single_cell_adata, tmp_path):
    out_dir = tmp_path / "rctd_ref"
    res = export_for_rctd(mock_single_cell_adata, out_dir)

    assert res["counts"].exists()
    assert res["metadata"].exists()

    meta_df = pd.read_csv(res["metadata"], index_col=0)
    assert "cell_type" in meta_df.columns
    assert "nUMI" in meta_df.columns