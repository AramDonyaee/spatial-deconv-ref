import logging
import warnings

# Completely silence Census log notices and stable release warnings globally
logging.getLogger("cellxgene_census").setLevel(logging.ERROR)
logging.getLogger("tiledbsoma").setLevel(logging.ERROR)
warnings.filterwarnings("ignore", message=r".*Specify 'census_version=.*")
warnings.filterwarnings("ignore", message=r".*column_names.*is deprecated.*")
warnings.filterwarnings("ignore", category=FutureWarning, module="cellxgene_census")

from .profiler import SpatialProfile, extract_spatial_profile, detect_species
from .query import find_candidate_datasets
from .scorer import fetch_candidate_anndata, score_candidate
from .exporter import export_for_cell2location, export_for_rctd
from .pipeline import ReferenceSuggester

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