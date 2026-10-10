"""cargofinder: candidate cargo discovery for encapsulin nanocompartments."""

from cargofinder.cluster import ClusterResult, cluster_candidates
from cargofinder.config import RunConfig, ScoringWeights
from cargofinder.features import PeptideFeatures, compute_peptide_features
from cargofinder.genecalling import call_genes_from_genome_fasta, call_genes_pyrodigal
from cargofinder.io import Gene, GenomeData, load_protein_and_gff3
from cargofinder.neighborhood import NeighborCandidate, extract_neighborhood
from cargofinder.pipeline import RunResult, run_pipeline
from cargofinder.plotting import plot_neighborhood
from cargofinder.profiles import fetch_default_profiles, fetch_pfam_profile
from cargofinder.reporting import (
    export_candidates_json,
    export_candidates_tsv,
    generate_run_summary,
)
from cargofinder.scoring import ScoredCandidate, score_candidates
from cargofinder.shells import ShellHit, detect_shells

__version__ = "0.1.0"

__all__ = [
    "Gene",
    "GenomeData",
    "ShellHit",
    "NeighborCandidate",
    "ClusterResult",
    "ScoredCandidate",
    "PeptideFeatures",
    "RunConfig",
    "ScoringWeights",
    "RunResult",
    "load_protein_and_gff3",
    "call_genes_pyrodigal",
    "call_genes_from_genome_fasta",
    "detect_shells",
    "extract_neighborhood",
    "cluster_candidates",
    "compute_peptide_features",
    "score_candidates",
    "plot_neighborhood",
    "export_candidates_tsv",
    "export_candidates_json",
    "generate_run_summary",
    "fetch_pfam_profile",
    "fetch_default_profiles",
    "run_pipeline",
]
