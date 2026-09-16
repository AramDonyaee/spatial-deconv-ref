"""Command-line interface for spatial_deconv_ref."""
import argparse
import sys
from pathlib import Path
import scanpy as sc
from .pipeline import ReferenceSuggester


def main():
    parser = argparse.ArgumentParser(
        description="Recommend and export single-cell references for spatial deconvolution."
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    # Subcommand: recommend
    rec_parser = subparsers.add_parser("recommend", help="Find and rank references from Census.")
    rec_parser.add_argument("--spatial", "-s", required=True, help="Path to spatial data (.h5, .h5ad)")
    rec_parser.add_argument("--tissue", "-t", required=True, help="Tissue name (e.g. 'colon', 'lymph node')")
    rec_parser.add_argument("--organism", "-o", default="auto", help="Organism: 'auto', 'homo_sapiens', 'mus_musculus'")
    rec_parser.add_argument("--top-n", type=int, default=3, help="Number of candidates to score (default: 3)")
    rec_parser.add_argument("--min-cells", type=int, default=1000, help="Minimum cells in reference candidate")
    rec_parser.add_argument("--export", "-e", choices=["cell2location", "rctd"], default=None, help="Export format for winner")
    rec_parser.add_argument("--out", required=False, help="Output file path (c2l) or directory (rctd)")
    rec_parser.add_argument("--max-cells", type=int, default=3000, help="Max cells to export (default: 3000)")

    args = parser.parse_args()

    if args.command == "recommend":
        print(f"Loading spatial data from {args.spatial}...")
        spatial_path = Path(args.spatial)
        if spatial_path.suffix == ".h5ad":
            adata = sc.read_h5ad(spatial_path)
        else:
            adata = sc.read_10x_h5(spatial_path)
        adata.var_names_make_unique()

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
                print("Error: --out is required when --export is specified.", file=sys.stderr)
                sys.exit(1)
            best_id = ranked.iloc[0]["dataset_id"]
            print(f"\nExporting top candidate ({best_id}) to {args.out}...")
            suggester.fetch_and_export(
                dataset_id=best_id,
                export_format=args.export,
                output_path=args.out,
                max_cells=args.max_cells
            )
            print("Export complete!")


if __name__ == "__main__":
    main()