"""Spatial profile extraction and species auto-detection."""
from __future__ import annotations

import re
from dataclasses import dataclass
import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata import AnnData


@dataclass
class SpatialProfile:
    """Summary profile of a spatial transcriptomics dataset."""
    organism: str
    genes: list[str]
    pseudobulk: pd.Series
    n_spots: int
    mean_counts_per_spot: float


def detect_species(adata: AnnData, sample_size: int = 500) -> str:
    """
    Auto-detect species (homo_sapiens vs. mus_musculus) using:
    1. Ensembl Gene ID prefix (adata.var['gene_ids'] or adata.var_names)
    2. Reference genome metadata (adata.var['genome'])
    3. HGNC (ALL UPPERCASE) vs MGI (Title Case) gene symbol naming conventions.
    """
    # 1. Ensembl ID check
    id_series = None
    if "gene_ids" in adata.var:
        id_series = adata.var["gene_ids"].dropna().astype(str)
    elif adata.var_names.str.startswith("ENS").any():
        id_series = adata.var_names.to_series()

    if id_series is not None and not id_series.empty:
        ensg_ratio = float(id_series.str.startswith("ENSG").mean())
        ensmusg_ratio = float(id_series.str.startswith("ENSMUSG").mean())
        if ensg_ratio > 0.5:
            return "homo_sapiens"
        if ensmusg_ratio > 0.5:
            return "mus_musculus"

    # 2. Reference Genome metadata check
    if "genome" in adata.var:
        genomes = adata.var["genome"].dropna().astype(str).str.lower().unique()
        for g in genomes:
            if any(h in g for h in ["grch", "hg38", "hg19"]):
                return "homo_sapiens"
            if any(m in g for m in ["mm10", "mm39", "grcm"]):
                return "mus_musculus"

    # 3. Gene symbol casing heuristics (HGNC uppercase vs MGI titlecase)
    genes = [str(g) for g in adata.var_names[:sample_size] if isinstance(g, str) and len(str(g)) > 1]
    clean_genes = [g for g in genes if re.match(r"^[A-Za-z0-9\-]+$", g)]

    if not clean_genes:
        raise ValueError("Could not auto-detect species: No standard alphanumeric gene symbols found.")

    n_upper = sum(1 for g in clean_genes if g.isupper())
    n_title = sum(1 for g in clean_genes if g[0].isupper() and g[1:].islower())

    ratio_upper = n_upper / len(clean_genes)
    ratio_title = n_title / len(clean_genes)

    if ratio_upper >= 0.65:
        return "homo_sapiens"
    elif ratio_title >= 0.65:
        return "mus_musculus"

    raise ValueError(
        f"Unable to auto-detect species (uppercase ratio: {ratio_upper:.2f}, "
        f"titlecase ratio: {ratio_title:.2f}). Please specify organism='homo_sapiens' "
        f"or organism='mus_musculus' explicitly."
    )


def extract_spatial_profile(
    adata: AnnData,
    organism: str = "auto",
    gene_column: str | None = None
) -> SpatialProfile:
    """
    Summarize spatial AnnData into an atlas query profile.
    Computes library-size normalized pseudobulk (log1p CPM).
    """
    if organism.lower() == "auto":
        resolved_organism = detect_species(adata)
        print(f"      [Species Auto-Detection] Identified organism: '{resolved_organism}'")
    else:
        resolved_organism = organism.lower().replace(" ", "_")

    if resolved_organism not in ["homo_sapiens", "mus_musculus"]:
        raise ValueError(
            f"Unsupported organism '{resolved_organism}'. "
            f"CELLxGENE Census currently supports 'homo_sapiens' and 'mus_musculus'."
        )

    if gene_column and gene_column in adata.var:
        genes = adata.var[gene_column].astype(str).tolist()
    else:
        genes = adata.var_names.astype(str).tolist()

    X = adata.X
    if sp.issparse(X):
        bulk_counts = np.asarray(X.sum(axis=0)).ravel()
        spot_sums = np.asarray(X.sum(axis=1)).ravel()
    else:
        bulk_counts = np.asarray(X).sum(axis=0).ravel()
        spot_sums = np.asarray(X).sum(axis=1).ravel()

    total_counts = float(bulk_counts.sum())
    if total_counts > 0:
        cpm = (bulk_counts / total_counts) * 1e6
        log_pseudobulk = np.log1p(cpm)
    else:
        log_pseudobulk = bulk_counts

    pseudobulk_series = pd.Series(log_pseudobulk, index=genes)
    pseudobulk_series = pseudobulk_series.groupby(pseudobulk_series.index).mean()

    return SpatialProfile(
        organism=resolved_organism,
        genes=list(pseudobulk_series.index),
        pseudobulk=pseudobulk_series,
        n_spots=int(adata.n_obs),
        mean_counts_per_spot=float(np.mean(spot_sums)) if len(spot_sums) > 0 else 0.0
    )