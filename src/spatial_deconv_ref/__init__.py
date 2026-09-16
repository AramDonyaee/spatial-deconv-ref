"""screfsuggest: Automated single-cell reference recommendation for spatial transcriptomics."""

from .profiler import SpatialProfile, extract_spatial_profile, detect_species
from .query import find_candidate_datasets
from .scorer import fetch_candidate_anndata, score_candidate
from .exporter import export_for_cell2location, export_for_rctd
from .pipeline import ReferenceSuggester

import logging
import warnings


# 1. Silence cellxgene_census logger INFO messages
logging.getLogger("cellxgene_census").setLevel(logging.WARNING)

# 2. Silence the UserWarning about specifying census_version
warnings.filterwarnings("ignore", message=r".*Specify 'census_version=.*")
warnings.filterwarnings("ignore", category=FutureWarning, module="cellxgene_census")

__version__ = "0.1.0"

__all__ = [
    "SpatialProfile",
    "extract_spatial_profile",
    "detect_species",
    "find_candidate_datasets",
    "fetch_candidate_anndata",
    "score_candidate",
    "export_for_cell2location",
    "export_for_rctd",
    "ReferenceSuggester",
]