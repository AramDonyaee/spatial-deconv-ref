"""CZ CELLxGENE Census discovery module."""
from __future__ import annotations

import cellxgene_census
import pandas as pd


def find_candidate_datasets(
    tissue: str,
    organism: str = "homo_sapiens",
    min_cells: int = 2000,
    census_version: str = "stable"
) -> pd.DataFrame:
    """
    Search CELLxGENE Census for candidate single-cell datasets matching tissue.
    Reads minimal columns to ensure fast response over S3.
    """
    tissue_clean = tissue.strip().lower()
    print(f"      Connecting to Census S3 storage and scanning '{tissue_clean}' metadata...")

    with cellxgene_census.open_soma(census_version=census_version) as census:
        exp = census["census_data"][organism]

        # Primary search: broad tissue_general ontology
        obs_df = exp.obs.read(
            column_names=["dataset_id", "cell_type", "assay"],
            value_filter=f"tissue_general == '{tissue_clean}'"
        ).concat().to_pandas()

        # Fallback search: exact tissue name
        if obs_df.empty:
            obs_df = exp.obs.read(
                column_names=["dataset_id", "cell_type", "assay"],
                value_filter=f"tissue == '{tissue_clean}'"
            ).concat().to_pandas()

        if obs_df.empty:
            return pd.DataFrame()

        print(f"      Retrieved {len(obs_df):,} cell records. Aggregating candidate statistics...")

        summary_records = []
        for dataset_id, group in obs_df.groupby("dataset_id", observed=False):
            n_cells = len(group)
            if n_cells < min_cells:
                continue

            summary_records.append({
                "dataset_id": dataset_id,
                "n_cells": n_cells,
                "n_cell_types": int(group["cell_type"].nunique()),
                "dominant_assay": str(group["assay"].mode().iloc[0]) if not group["assay"].empty else "unknown"
            })

        candidate_df = pd.DataFrame(summary_records)
        if candidate_df.empty:
            return candidate_df

        # Fetch readable study titles and collection info
        datasets_info = census["census_info"]["datasets"].read().concat().to_pandas()
        candidate_df = candidate_df.merge(
            datasets_info[["dataset_id", "dataset_title", "collection_name"]],
            on="dataset_id",
            how="left"
        )
        return candidate_df.sort_values(by="n_cells", ascending=False).reset_index(drop=True)