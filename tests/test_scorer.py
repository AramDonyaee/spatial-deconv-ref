"""Tests for candidate scoring logic."""
from screfsuggest.profiler import extract_spatial_profile
from screfsuggest.scorer import score_candidate


def test_score_candidate(mock_spatial_adata, mock_single_cell_adata):
    profile = extract_spatial_profile(mock_spatial_adata, organism="homo_sapiens")
    scores = score_candidate(mock_single_cell_adata, profile, min_cells_per_type=10)

    assert "composite_score" in scores
    assert "concordance" in scores
    assert "gene_overlap" in scores
    assert 0.0 <= scores["composite_score"] <= 1.0
    assert 0.0 <= scores["gene_overlap"] <= 1.0
    assert scores["n_valid_cell_types"] == 3
    assert scores["n_shared_genes"] > 0