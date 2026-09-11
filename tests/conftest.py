"""Pytest fixtures with synthetic data (runs completely offline)."""
import pytest
import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata import AnnData


@pytest.fixture
def mock_spatial_adata():
    """Generates miniature synthetic spatial AnnData (Human)."""
    np.random.seed(42)
    n_spots = 40
    genes = ["CD3D", "CD4", "CD8A", "MS4A1", "PECAM1", "EPCAM", "MALAT1", "ACTB"]
    counts = np.random.poisson(lam=5.0, size=(n_spots, len(genes)))

    return AnnData(
        X=sp.csr_matrix(counts),
        obs=pd.DataFrame(index=[f"spot_{i}" for i in range(n_spots)]),
        var=pd.DataFrame(index=genes)
    )


@pytest.fixture
def mock_mouse_spatial_adata():
    """Generates miniature synthetic spatial AnnData (Mouse)."""
    np.random.seed(42)
    n_spots = 30
    genes = ["Cd3d", "Cd4", "Cd8a", "Ms4a1", "Pecam1", "Epcam", "Malat1", "Actb"]
    counts = np.random.poisson(lam=4.0, size=(n_spots, len(genes)))

    return AnnData(
        X=sp.csr_matrix(counts),
        obs=pd.DataFrame(index=[f"spot_{i}" for i in range(n_spots)]),
        var=pd.DataFrame(index=genes)
    )


@pytest.fixture
def mock_single_cell_adata():
    """Generates miniature synthetic single-cell reference AnnData."""
    np.random.seed(42)
    n_cells = 90
    genes = ["CD3D", "CD4", "CD8A", "MS4A1", "PECAM1", "GAPDH", "ACTB"]
    counts = np.random.poisson(lam=3.0, size=(n_cells, len(genes)))
    cell_types = ["T cell"] * 30 + ["B cell"] * 30 + ["Endothelial"] * 30

    return AnnData(
        X=sp.csr_matrix(counts),
        obs=pd.DataFrame({
            "cell_type": cell_types,
            "assay": ["10x 3' v3"] * n_cells,
            "disease": ["normal"] * n_cells
        }, index=[f"cell_{i}" for i in range(n_cells)]),
        var=pd.DataFrame(index=genes)
    )