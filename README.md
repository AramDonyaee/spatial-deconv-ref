# spatial-deconv-ref

[![PyPI version](https://img.shields.io/badge/pypi-v0.1.0-blue.svg)](https://pypi.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Supported OS](https://img.shields.io/badge/OS-Linux%20%7C%20macOS%20%7C%20WSL2-success.svg)]()

> **Automated single-cell reference recommendation, scoring, and formatting for spatial transcriptomics deconvolution.**

---

## Overview

Spatial transcriptomics deconvolution tools (**cell2location**, **RCTD**, **Stereoscope**) depend heavily on the quality and biological match of the single-cell RNA-seq reference. An unmatched reference leads to biased cell-type proportions and missed rare cell states.

`spatial-deconv-ref` connects your spatial data directly to the **CZ CELLxGENE Census**, evaluates candidate references using a multi-factor score, and exports the top-ranked reference formatted for immediate downstream deconvolution.

---

## Key Features

- **Species Auto-Detection:** Automatically infers whether your spatial panel is Human (*Homo sapiens*) or Mouse (*Mus musculus*) using HGNC vs. MGI casing conventions, Ensembl IDs, and genome metadata.
- **Atlas-Scale Discovery:** Queries millions of cells across CZ CELLxGENE Census in seconds using indexed metadata streaming.
- **Multi-Factor Scoring Index:**
  1. **Pseudobulk Concordance:** Spearman rank correlation of log1p-CPM expression vectors between spatial spots and candidate cells.
  2. **Gene Panel Overlap:** Jaccard overlap fraction (critical for targeted platforms like 10x Xenium and NanoString CosMx).
  3. **Cell-Type Diversity:** Penalizes monocultures and ensures clusters have sufficient cell representation ($\ge 15$ cells).
- **Turnkey Exporters:**
  - **`cell2location`:** Exports `.h5ad` with unnormalized integer counts, non-empty genes/cells, and categorical cell types.
  - **`RCTD` / Seurat:** Exports `ref_counts.h5ad` and `ref_annotations.csv` with UMI counts.
- **Dual Interface:** Full Python API for Jupyter/Scanpy workflows + standalone CLI for bash/R pipelines.

---

## Installation

Because `cellxgene-census` relies on `tiledbsoma`, a **Linux**, **macOS**, or **Windows Subsystem for Linux (WSL2)** environment is required:

```bash
git clone https://github.com/your-username/spatial-deconv-ref.git
cd spatial-deconv-ref
pip install -e .
Usage: Python API
code
Python
import scanpy as sc
from spatial_deconv_ref import ReferenceSuggester

# 1. Load your spatial dataset (Visium, Visium HD, Xenium, etc.)
adata_spatial = sc.read_10x_h5("filtered_feature_bc_matrix.h5")
adata_spatial.var_names_make_unique()

# 2. Recommend candidate references (species is auto-detected)
suggester = ReferenceSuggester(organism="auto")
ranked_df = suggester.recommend(
    spatial_adata=adata_spatial,
    tissue="lymph node",
    top_n_eval=3,
    min_cells=1000
)

# 3. View the recommendation leaderboard
print(ranked_df[["dataset_title", "composite_score", "concordance", "n_valid_cell_types", "portal_url"]])

# 4. Export the winning reference
best_id = ranked_df.iloc[0]["dataset_id"]

# For cell2location:
suggester.fetch_and_export(
    dataset_id=best_id,
    export_format="cell2location",
    output_path="c2l_reference.h5ad",
    max_cells=3000
)

# For RCTD:
suggester.fetch_and_export(
    dataset_id=best_id,
    export_format="rctd",
    output_path="./rctd_reference/",
    max_cells=3000
)
Usage: Command-Line Interface (CLI)
Ideal for command-line users, R researchers running RCTD, and workflow managers (Nextflow / Snakemake):
code
Bash
# Basic candidate ranking
spatial-deconv-ref recommend \
    --spatial filtered_feature_bc_matrix.h5 \
    --tissue "lymph node"

# Rank and export directly for RCTD (R ingest)
spatial-deconv-ref recommend \
    --spatial filtered_feature_bc_matrix.h5 \
    --tissue "lymph node" \
    --export rctd \
    --out ./rctd_ref/ \
    --max-cells 3000

# Rank and export directly for cell2location
spatial-deconv-ref recommend \
    --spatial filtered_feature_bc_matrix.h5 \
    --tissue "colon" \
    --export cell2location \
    --out ./c2l_ref.h5ad \
    --max-cells 3000

## CLI Parameters

| Flag | Shorthand | Description | Default |
|------|-----------|-------------|---------|
| `--spatial` | `-s` | Path to spatial matrix file (`.h5`, `.h5ad`) | Required |
| `--tissue` | `-t` | Target tissue ontology term (e.g. `'colon'`, `'lymph node'`) | Required |
| `--organism` | `-o` | Organism: `'auto'`, `'homo_sapiens'`, `'mus_musculus'` | `'auto'` |
| `--top-n` |  | Number of candidate references to score | `3` |
| `--min-cells` |  | Minimum cells required for candidate dataset | `1000` |
| `--export` | `-e` | Export target: `'cell2location'` or `'rctd'` | `None` |
| `--out` |  | Output destination path or folder | `None` |
| `--max-cells` |  | Max cells to download and export for reference | `3000` |

## Scoring Formula

Candidates are scored using a composite index:

$$
\text{Composite Score} = 0.40 \cdot \max(r_s, 0) + 0.35 \cdot S_{\text{overlap}} + 0.25 \cdot S_{\text{diversity}}
$$

Where:

- $r_s$ (**Concordance**): Spearman rank correlation of normalized pseudobulk expression across shared features.
- $S_{\text{overlap}}$ (**Gene Overlap**): Fraction of spatial panel genes captured in the reference.
- $S_{\text{diversity}}$ (**Cluster Diversity**): Log-normalized count of distinct cell types with $\ge 15$ cells.