"""
run_demo.py: Live demo using Scanpy's built-in 10x Visium Human Lymph Node dataset.
"""
from pathlib import Path
import scanpy as sc
from screfsuggest import ReferenceSuggester, extract_spatial_profile


def main():
    print("=" * 60)
    print(" screfsuggest: End-to-End Pipeline Demo")
    print("=" * 60)

    print("\n[1/4] Loading 10x Visium Human Lymph Node sample...")
    adata_spatial = sc.datasets.visium_sge(sample_id="V1_Human_Lymph_Node")
    adata_spatial.var_names_make_unique()
    print(f"      Loaded: {adata_spatial.n_obs} spots x {adata_spatial.n_vars} genes.")

    print("\n[2/4] Extracting spatial profile (auto-detecting species)...")
    profile = extract_spatial_profile(adata_spatial, organism="auto")
    print(f"      Profile: {profile.organism}, {len(profile.genes)} genes, {profile.mean_counts_per_spot:.1f} mean UMI/spot.")

    print("\n[3/4] Querying CELLxGENE Census for 'lymph node' candidates...")
    suggester = ReferenceSuggester(organism="auto")
    ranked_df = suggester.recommend(
        spatial_adata=adata_spatial,
        tissue="lymph node",
        top_n_eval=3,
        min_cells=2000
    )

    print("\n--- Recommendation Leaderboard ---")
    cols = ["dataset_title", "composite_score", "gene_overlap", "concordance", "n_valid_cell_types"]
    print(ranked_df[[c for c in cols if c in ranked_df.columns]].to_string(index=False))

    best_id = ranked_df.iloc[0]["dataset_id"]
    out_dir = Path("exported_references")
    out_dir.mkdir(exist_ok=True)

    print(f"\n[4/4] Exporting winning dataset: {best_id}")
    c2l_path = out_dir / "lymph_node_c2l_ref.h5ad"
    suggester.fetch_and_export(best_id, "cell2location", c2l_path)
    print(f"      Saved cell2location reference: {c2l_path}")

    rctd_dir = out_dir / "lymph_node_rctd_ref"
    suggester.fetch_and_export(best_id, "rctd", rctd_dir)
    print(f"      Saved RCTD reference: {rctd_dir}/")
    print("\nDemo finished successfully!")


if __name__ == "__main__":
    main()