"""High-level recommendation pipeline orchestrator."""
from __future__ import annotations

import io
import sys
import contextlib
from pathlib import Path
import pandas as pd
from tqdm import tqdm
from anndata import AnnData

from .profiler import extract_spatial_profile
from .query import find_candidate_datasets
from .scorer import fetch_candidate_anndata, score_candidate
from .exporter import export_for_cell2location, export_for_rctd


class ReferenceSuggester:
    """Engine to discover, rank, and export scRNA-seq references for deconvolution."""

    def __init__(self, organism: str = "auto"):
        self.organism = organism

    def recommend(
        self,
        spatial_adata: AnnData,
        tissue: str,
        top_n_eval: int = 3,
        min_cells: int = 1000,
        max_cells_per_candidate: int = 2500
    ) -> pd.DataFrame:
        """Query Census and score top candidates with clean progress tracking."""
        profile = extract_spatial_profile(spatial_adata, organism=self.organism)
        self.organism = profile.organism

        # 1. Candidate discovery
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            candidates = find_candidate_datasets(
                tissue=tissue, 
                organism=self.organism, 
                min_cells=min_cells
            )

        if candidates.empty:
            raise ValueError(f"No candidate datasets found in Census for tissue '{tissue}' ({self.organism}).")

        top_candidates = candidates.head(top_n_eval).copy()

        print(f"\n      Top {len(top_candidates)} Candidate Datasets Identified:")
        for idx, (_, row) in enumerate(top_candidates.iterrows()):
            d_id = row["dataset_id"]
            title = str(row.get("dataset_title", "Unknown"))
            cells = row.get("n_cells", 0)
            url = f"https://cellxgene.cziscience.com/e/{d_id}.cxg/"
            print(f"        [{idx + 1}] {title}")
            print(f"            Cells: {cells:,} | URL: {url}")

        scored_records = []
        print(f"\n      Scoring candidate references against spatial profile...")

        # Single cleanly updating progress bar
        with tqdm(
            total=len(top_candidates),
            desc="      Overall Progress",
            unit="atlas",
            leave=True,
            dynamic_ncols=True,
            file=sys.stdout
        ) as pbar:
            for _, row in top_candidates.iterrows():
                d_id = row["dataset_id"]
                short_title = str(row.get("dataset_title", d_id))[:30]

                # Step 1: Downloading slice
                pbar.set_postfix_str(f"[{short_title}] 1/2 Fetching S3 slice...")
                try:
                    # Redirect stdout/stderr to prevent census logs from duplicating the bar
                    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
                        cand_adata = fetch_candidate_anndata(
                            dataset_id=d_id, 
                            organism=self.organism, 
                            genes_subset=profile.genes,
                            max_cells=max_cells_per_candidate
                        )

                    # Step 2: Scoring
                    pbar.set_postfix_str(f"[{short_title}] 2/2 Correlating pseudobulk...")
                    score_dict = score_candidate(cand_adata, profile)
                    score_dict["dataset_id"] = d_id
                    score_dict["dataset_title"] = row.get("dataset_title", "Unknown")
                    score_dict["portal_url"] = f"https://cellxgene.cziscience.com/e/{d_id}.cxg/"
                    score_dict["dominant_assay"] = row.get("dominant_assay", "Unknown")
                    score_dict["n_cells"] = row.get("n_cells", 0)
                    scored_records.append(score_dict)

                except Exception as exc:
                    pbar.write(f"\n      [Warning] Skipped candidate {d_id}: {exc}")

                pbar.update(1)
            pbar.set_postfix_str("Complete")

        if not scored_records:
            raise RuntimeError("Failed to score candidate references from Census.")

        ranked_df = pd.DataFrame(scored_records).sort_values(
            by="composite_score", ascending=False
        ).reset_index(drop=True)
        return ranked_df

    def fetch_and_export(
        self,
        dataset_id: str,
        export_format: str,
        output_path: str | Path,
        max_cells: int | None = 3000
    ) -> Path | dict[str, Path]:
        """Fetch winning reference and export it cleanly."""
        target_organism = "homo_sapiens" if self.organism == "auto" else self.organism
        cap_str = f"{max_cells:,} cells" if max_cells else "all cells"
        print(f"      Downloading reference data ({cap_str}) for export...")

        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            cand_adata = fetch_candidate_anndata(
                dataset_id=dataset_id,
                organism=target_organism,
                max_cells=max_cells
            )

        fmt = export_format.strip().lower()
        if fmt == "cell2location":
            return export_for_cell2location(cand_adata, output_path)
        elif fmt == "rctd":
            return export_for_rctd(cand_adata, output_path)
        else:
            raise ValueError(f"Unsupported format '{export_format}'. Choose 'cell2location' or 'rctd'.")