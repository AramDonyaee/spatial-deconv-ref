"""
test_my_data.py: Run screfsuggest on local Visium HD or Xenium datasets.
"""
import warnings
from pathlib import Path
import scanpy as sc
from screfsuggest import ReferenceSuggester, extract_spatial_profile

# ==============================================================================
# CONFIGURATION
# ==============================================================================
DATA_TYPE = "visium_hd"   # Options: "visium_hd" or "xenium"
TISSUE = "colon"          # e.g., 'colon', 'breast', 'brain', 'lung', 'lymph node'

if DATA_TYPE == "visium_hd":
    # Point to filtered_feature_bc_matrix.h5 inside square_016um or square_008um
    DATA_PATH = Path("/mnt/e/FDM Downloads/Visium HD datasets/binned_outputs/square_016um/filtered_feature_bc_matrix.h5")
elif DATA_TYPE == "xenium":
    # Point to cell_feature_matrix.h5 (or cell_feature_matrix.tar.gz)
    DATA_PATH = Path("/mnt/d/Xenium-Adenocarcinoma/cell_feature_matrix.h5")

OUTPUT_DIR = Path("exported_references")
OUTPUT_DIR.mkdir(exist_ok=True)
# ==============================================================================


def load_dataset(path: Path, data_type: str):
    print(f"\n[1/4] Loading {data_type.upper()} dataset from: {path}")
    if not path.exists():
        raise FileNotFoundError(f"File not found: {path}")

    # Suppress Scanpy's duplicate var_names warning during read
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=UserWarning, message=".*Variable names are not unique.*")
        adata = sc.read_10x_h5(path)

    # Immediately make them unique
    adata.var_names_make_unique()

    if data_type == "xenium":
        is_real_gene = ~adata.var_names.str.startswith(
            ("Blank-", "NegControlProbe-", "NegControlCodeword-", "BLANK_")
        )
        adata = adata[:, is_real_gene].copy()
        print(f"      Filtered out probe controls; remaining targeted genes: {adata.n_vars}")

    print(f"      Loaded: {adata.n_obs:,} spots/cells x {adata.n_vars:,} unique genes.")
    return adata


def main():
    print("=" * 65)
    print(f" screfsuggest: Local {DATA_TYPE.upper()} Pipeline")
    print("=" * 65)

    adata_spatial = load_dataset(DATA_PATH, DATA_TYPE)

    print("\n[2/4] Profiling spatial expression vector...")
    profile = extract_spatial_profile(adata_spatial, organism="auto")

    print(f"\n[3/4] Querying CELLxGENE Census for '{TISSUE}' candidates...")
    suggester = ReferenceSuggester(organism="auto")
    
    # min_cells=1000 allows multiple mouse colon datasets to qualify
    ranked_df = suggester.recommend(
        spatial_adata=adata_spatial,
        tissue=TISSUE,
        top_n_eval=3,
        min_cells=1000
    )

    print("\n--- Recommendation Leaderboard ---")
    cols = ["dataset_title", "composite_score", "concordance", "n_valid_cell_types", "n_cells", "portal_url"]
    print(ranked_df[[c for c in cols if c in ranked_df.columns]].to_string(index=False))

    best_dataset = ranked_df.iloc[0]
    best_id = best_dataset["dataset_id"]
    print(f"\n[4/4] Exporting winning dataset: {best_dataset.get('dataset_title', best_id)}")

    c2l_file = OUTPUT_DIR / f"{DATA_TYPE}_{TISSUE}_c2l_ref.h5ad"
    # Capped at 5,000 cells to download in seconds (~20-40 MB)
    suggester.fetch_and_export(best_id, "cell2location", c2l_file, max_cells=5000)
    print(f"      Exported cell2location file: {c2l_file}")

    rctd_dir = OUTPUT_DIR / f"{DATA_TYPE}_{TISSUE}_rctd_ref"
    suggester.fetch_and_export(best_id, "rctd", rctd_dir, max_cells=5000)
    print(f"      Exported RCTD folder: {rctd_dir}/")
    print("\nProcessing complete! All warnings silenced and files ready for deconvolution.")


if __name__ == "__main__":
    main()