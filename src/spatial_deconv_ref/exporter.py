"""Export single-cell references formatted for deconvolution engines."""
from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata import AnnData


def export_for_cell2location(
    adata: AnnData,
    output_path: str | Path,
    cell_type_col: str = "cell_type"
) -> Path:
    """
    Format reference for cell2location:
    - Retains unnormalized counts
    - Removes zero-count genes and cells
    - Casts cell_type to categorical
    - Clears index names to prevent AnnData HDF5 writer collisions
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    ref = adata.copy()

    # Clear index names to prevent ValueError: DataFrame.index.name ('feature_name')
    ref.var.index.name = None
    ref.obs.index.name = None

    if sp.issparse(ref.X):
        ref.X = ref.X.astype(np.float32)
    else:
        ref.X = np.asarray(ref.X, dtype=np.float32)

    # Filter empty cells and empty genes
    gene_counts = np.asarray((ref.X > 0).sum(axis=0)).ravel()
    cell_counts = np.asarray((ref.X > 0).sum(axis=1)).ravel()

    ref = ref[cell_counts > 0, gene_counts > 0].copy()
    ref.obs[cell_type_col] = ref.obs[cell_type_col].astype("category")

    ref.write_h5ad(output_path)
    return output_path


def export_for_rctd(
    adata: AnnData,
    output_dir: str | Path,
    cell_type_col: str = "cell_type"
) -> dict[str, Path]:
    """
    Format reference for RCTD:
    - Exports reference count matrix (ref_counts.h5ad)
    - Exports barcode-to-cell-type annotations (ref_annotations.csv)
    """
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    counts_file = out_dir / "ref_counts.h5ad"
    meta_file = out_dir / "ref_annotations.csv"

    if sp.issparse(adata.X):
        umi_counts = np.asarray(adata.X.sum(axis=1)).ravel()
    else:
        umi_counts = np.asarray(adata.X).sum(axis=1).ravel()

    meta_df = pd.DataFrame({
        "barcode": adata.obs_names,
        "cell_type": adata.obs[cell_type_col].values,
        "nUMI": umi_counts
    })
    meta_df.set_index("barcode", inplace=True)
    meta_df.to_csv(meta_file)

    ref_counts = adata.copy()
    ref_counts.var.index.name = None
    ref_counts.obs.index.name = None
    ref_counts.write_h5ad(counts_file)

    return {"counts": counts_file, "metadata": meta_file}