"""Command-line interface for spatial_deconv_ref."""
from __future__ import annotations

import argparse
import sys
import warnings
from pathlib import Path
import scanpy as sc
from .pipeline import ReferenceSuggester


def main():
    parser = argparse.ArgumentParser(
        prog="spatial-deconv-ref",
        description="Recommend and export single-cell references for spatial deconvolution."
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    # Subcommand: recommend
    rec_parser = subparsers.add_parser("recommend", help="Find, rank, and optionally export references from Census.")
    rec_parser.add_argument("--spatial", "-s", required=True, help="Path to spatial input file (.h5, .h5ad)")
    rec_parser.add_argument("--tissue", "-t", required=True, help="Target tissue name (e.g. 'colon', 'lymph node')")
    rec_parser.add_argument("--organism", "-o", default="auto", help="Organism: 'auto', 'homo_sapiens', 'mus_musculus'")
    rec_parser.add_argument("--top-n", type=int, default=3, help="Number of candidate datasets to evaluate (default: 3)")
    rec_parser.add_argument("--min-cells", type=int, default=1000, help="Minimum cells required for candidate dataset (default: 1000)")
    rec_parser.add_argument("--export", "-e", choices=["cell2location", "rctd"], default=None, help="Export format for the winning reference")
    rec_parser.add_argument("--out", required=False, help="Output destination (file path for cell2location, directory for rctd)")
    rec_parser.add_argument("--max-cells", type=int, default=3000, help="Max cells to export (default: 3000)")

    args = parser.parse_args()

    if args.command == "recommend":
        spatial_path = Path(args.spatial)
        if not spatial_path.exists():
            print(f"Error: Spatial file not found at: {spatial_path}", file=sys.stderr)
            sys.exit(1)

        print(f"[1/3] Loading spatial data from: {spatial_path}")
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=UserWarning)
            if spatial_path.suffix == ".h5ad":
                adata = sc.read_h5ad(spatial_path)
            else:
                adata = sc.read_10x_h5(spatial_path)
        adata.var_names_make_unique()
        print(f"      Loaded {adata.n_obs:,} spots/cells x {adata.n_vars:,} genes.")

        print(f"\n[2/3] Querying CELLxGENE Census for '{args.tissue}' candidates...")
        suggester = ReferenceSuggester(organism=args.organism)
        ranked = suggester.recommend(
            spatial_adata=adata,
            tissue=args.tissue,
            top_n_eval=args.top_n,
            min_cells=args.min_cells
        )

        cols = ["dataset_title", "composite_score", "concordance", "n_valid_cell_types", "n_cells"]
        print("\n--- Recommendation Leaderboard ---")
        print(ranked[[c for c in cols if c in ranked.columns]].to_string(index=False))

        # Handle export if requested
        if args.export:
            if not args.out:
                print("\nError: --out is required when --export is specified.", file=sys.stderr)
                sys.exit(1)

            best_dataset = ranked.iloc[0]
            best_id = best_dataset["dataset_id"]
            print(f"\n[3/3] Exporting top candidate '{best_dataset.get('dataset_title', best_id)}'...")

            suggester.fetch_and_export(
                dataset_id=best_id,
                export_format=args.export,
                output_path=args.out,
                max_cells=args.max_cells
            )
            print(f"      Export complete: {args.out}")


if __name__ == "__main__":
    main()