"""
run_demo.py: Live end-to-end verification script.
"""
from pathlib import Path
import scanpy as sc
from spatial_deconv_ref import ReferenceSuggester, extract_spatial_profile


def main():
    print("=" * 65)
    print(" spatial_deconv_ref: End-to-End Verification Demo")
    print("=" * 65)

    print("\n[1/3] Loading 10x Visium Human Lymph Node sample...")
    adata_spatial = sc.datasets.visium_sge(sample_id="V1_Human_Lymph_Node")
    adata_spatial.var_names_make_unique()
    print(f"      Loaded: {adata_spatial.n_obs:,} spots x {adata_spatial.n_vars:,} genes.")

    print("\n[2/3] Querying CELLxGENE Census for 'lymph node' candidates...")
    suggester = ReferenceSuggester(organism="auto")
    ranked_df = suggester.recommend(
        spatial_adata=adata_spatial,
        tissue="lymph node",
        top_n_eval=3,
        min_cells=1000
    )

    print("\n--- Recommendation Leaderboard ---")
    cols = ["dataset_title", "composite_score", "concordance", "n_valid_cell_types", "n_cells"]
    print(ranked_df[[c for c in cols if c in ranked_df.columns]].to_string(index=False))

    best_id = ranked_df.iloc[0]["dataset_id"]
    out_dir = Path("exported_references")
    out_dir.mkdir(exist_ok=True)

    print(f"\n[3/3] Exporting winning dataset: {best_id}")
    c2l_path = out_dir / "lymph_node_c2l_ref.h5ad"
    suggester.fetch_and_export(best_id, "cell2location", c2l_path, max_cells=2000)
    print(f"      Saved cell2location reference: {c2l_path}")

    rctd_dir = out_dir / "lymph_node_rctd_ref"
    suggester.fetch_and_export(best_id, "rctd", rctd_dir, max_cells=2000)
    print(f"      Saved RCTD reference: {rctd_dir}/")
    print("\nVerification finished successfully with zero errors!")


if __name__ == "__main__":
    main()