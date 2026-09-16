"""Scoring candidate single-cell references against spatial transcriptomics data."""
from __future__ import annotations

import warnings
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import spearmanr
import cellxgene_census
import tiledbsoma
from anndata import AnnData
from .profiler import SpatialProfile


def get_resilient_soma_context() -> tiledbsoma.SOMATileDBContext:
    """Create TileDB context optimized for high-latency or slow connections."""
    return tiledbsoma.SOMATileDBContext(
        tiledb_config={
            "vfs.s3.connect_timeout_ms": "60000",
            "vfs.s3.request_timeout_ms": "120000",
            "vfs.s3.max_parallel_ops": "4",
            "vfs.s3.connect_max_retries": "5",
        }
    )


def fetch_candidate_anndata(
    dataset_id: str,
    organism: str = "homo_sapiens",
    genes_subset: list[str] | None = None,
    max_cells: int | None = 2500,
    census_version: str = "stable"
) -> AnnData:
    """
    Fetch candidate AnnData from Census with targeted gene and cell downsampling.
    Uses indexed integer coordinates to avoid scanning the entire database.
    """
    context = get_resilient_soma_context()

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=UserWarning)
        warnings.filterwarnings("ignore", category=FutureWarning)

        with cellxgene_census.open_soma(census_version=census_version, context=context) as census:
            exp = census["census_data"][organism]

            # 1. Retrieve the integer soma_joinids for this candidate dataset
            obs_indexer = exp.obs.read(
                column_names=["soma_joinid"],
                value_filter=f"dataset_id == '{dataset_id}'"
            ).concat().to_pandas()

            if obs_indexer.empty:
                raise ValueError(f"No cells found in Census for dataset_id: {dataset_id}")

            join_ids = obs_indexer["soma_joinid"].values

            # 2. Subsample to max_cells locally
            if max_cells is not None and len(join_ids) > max_cells:
                np.random.seed(42)
                join_ids = np.random.choice(join_ids, size=max_cells, replace=False)

            join_ids.sort()
            obs_coords = [int(x) for x in join_ids]

            # 3. Filter down to spatial panel genes if provided (slice to 500 genes)
            var_filter = None
            if genes_subset:
                clean_genes = [g.replace("'", "") for g in genes_subset[:500] if isinstance(g, str)]
                if clean_genes:
                    gene_list_str = ", ".join(f"'{g}'" for g in clean_genes)
                    var_filter = f"feature_name in [{gene_list_str}]"

            # 4. Stream slice via direct indexed coordinates
            adata = cellxgene_census.get_anndata(
                census=census,
                organism=organism,
                obs_coords=obs_coords,
                var_value_filter=var_filter,
                obs_column_names=["cell_type", "assay", "tissue", "disease"]
            )

            adata.var_names = adata.var["feature_name"].astype(str)
            adata.var_names_make_unique()

            # Prevent index name collision in AnnData HDF5 writer
            adata.var.index.name = None
            adata.obs.index.name = None
            return adata


def score_candidate(
    candidate_adata: AnnData, 
    profile: SpatialProfile,
    min_cells_per_type: int = 15
) -> dict[str, float]:
    """Score candidate single-cell reference against spatial profile."""
    sc_genes = set(candidate_adata.var_names)
    spatial_genes = set(profile.genes)
    shared_genes = list(sc_genes.intersection(spatial_genes))
    jaccard_overlap = len(shared_genes) / max(len(spatial_genes), 1)

    ct_counts = candidate_adata.obs["cell_type"].value_counts()
    valid_ct = int((ct_counts >= min_cells_per_type).sum())
    diversity_score = float(min(np.log1p(valid_ct) / np.log1p(25), 1.0))

    concordance_score = 0.0
    if len(shared_genes) >= 20:
        sub_adata = candidate_adata[:, shared_genes]
        if sp.issparse(sub_adata.X):
            sc_bulk = np.asarray(sub_adata.X.sum(axis=0)).ravel()
        else:
            sc_bulk = np.asarray(sub_adata.X).sum(axis=0).ravel()

        sc_total = float(sc_bulk.sum())
        if sc_total > 0:
            sc_cpm = np.log1p((sc_bulk / sc_total) * 1e6)
            spatial_vec = profile.pseudobulk.loc[shared_genes].values
            corr, _ = spearmanr(sc_cpm, spatial_vec)
            concordance_score = float(np.nan_to_num(corr, nan=0.0))

    composite_score = (
        0.40 * max(concordance_score, 0.0) +
        0.35 * jaccard_overlap +
        0.25 * diversity_score
    )

    return {
        "composite_score": round(float(composite_score), 4),
        "concordance": round(float(concordance_score), 4),
        "gene_overlap": round(float(jaccard_overlap), 4),
        "cell_type_score": round(float(diversity_score), 4),
        "n_shared_genes": len(shared_genes),
        "n_valid_cell_types": valid_ct
    }