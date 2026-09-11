"""Tests for profile extraction and species auto-detection."""
from screfsuggest.profiler import extract_spatial_profile, detect_species


def test_detect_species_human(mock_spatial_adata):
    organism = detect_species(mock_spatial_adata)
    assert organism == "homo_sapiens"


def test_detect_species_mouse(mock_mouse_spatial_adata):
    organism = detect_species(mock_mouse_spatial_adata)
    assert organism == "mus_musculus"


def test_extract_spatial_profile_auto(mock_spatial_adata):
    profile = extract_spatial_profile(mock_spatial_adata, organism="auto")
    assert profile.organism == "homo_sapiens"
    assert profile.n_spots == 40
    assert len(profile.genes) == 8
    assert "CD3D" in profile.genes
    assert profile.mean_counts_per_spot > 0
    assert (profile.pseudobulk >= 0).all()