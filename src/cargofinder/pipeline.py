"""High-level Python API for encapsulin cargo candidate discovery."""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional

from cargofinder.cluster import cluster_candidates
from cargofinder.config import DEFAULT_SHELL_PROFILES, RunConfig, ScoringWeights
from cargofinder.features import load_motif_patterns
from cargofinder.genecalling import call_genes_from_genome_fasta
from cargofinder.io import GenomeData, load_protein_and_gff3
from cargofinder.neighborhood import NeighborCandidate, extract_neighborhood
from cargofinder.plotting import plot_neighborhood
from cargofinder.profiles import fetch_default_profiles
from cargofinder.reporting import (
    export_candidates_json,
    export_candidates_tsv,
    generate_run_summary,
    save_run_summary,
)
from cargofinder.scoring import ScoredCandidate, score_candidates
from cargofinder.shells import ShellHit, detect_shells

logger = logging.getLogger(__name__)


@dataclass
class RunResult:
    """Consolidated results of an end-to-end cargofinder run."""

    genome: GenomeData
    shells: List[ShellHit]
    candidates: List[ScoredCandidate]
    output_files: Dict[str, Path]
    summary: Dict[str, object]


def run_pipeline(
    input_path: Path | str,
    gff_path: Optional[Path | str] = None,
    out_dir: Path | str = "results",
    profile_dir: Optional[Path | str] = None,
    motif_file: Optional[Path | str] = None,
    window_size: int = 5,
    tail_length: int = 20,
    max_operon_gap_bp: int = 150,
    evalue_cutoff: float = 1e-5,
    min_shell_length_aa: int = 150,
    weights: Optional[ScoringWeights] = None,
    genome_id: Optional[str] = None,
    mmseqs_bin: Optional[Path | str] = None,
    generate_plots: bool = True,
) -> RunResult:
    """Execute the full cargofinder discovery pipeline.

    Parameters
    ----------
    input_path : Path or str
        Path to nucleotide FASTA (if gff_path is None) or protein FASTA
        (if gff_path is provided).
    gff_path : Path or str, optional
        Path to GFF3 annotation file.
    out_dir : Path or str
        Output directory for reports, tables, and figures.
    profile_dir : Path or str, optional
        Directory containing encapsulin shell .hmm profiles. If None, checks default
        locations or downloads standard profiles into out_dir/profiles.
    motif_file : Path or str, optional
        JSON file with targeting peptide regex patterns.
    window_size : int
        Number of genes to search upstream and downstream of each shell (default 5).
    tail_length : int
        Residue count for terminal peptide feature extraction (default 20).
    max_operon_gap_bp : int
        Max intergenic gap for likely operon co-transcription (default 150).
    evalue_cutoff : float
        HMMER E-value threshold for encapsulin shell detection (default 1e-5).
    min_shell_length_aa : int
        Minimum length in amino acids for encapsulin shell monomers (default 150).
    weights : ScoringWeights, optional
        Scoring weights. Defaults to standard weights.
    genome_id : str, optional
        Explicit genome identifier.
    mmseqs_bin : Path or str, optional
        Explicit path to mmseqs binary.
    generate_plots : bool
        Whether to generate neighborhood diagrams (default True).

    Returns
    -------
    RunResult
        Container with all results, candidate rankings, and output file paths.
    """
    in_p = Path(input_path)
    if not in_p.is_file():
        raise FileNotFoundError(f"Input sequence file not found: {input_path}")

    out_p = Path(out_dir)
    out_p.mkdir(parents=True, exist_ok=True)

    config = RunConfig(
        window_size=window_size,
        tail_length=tail_length,
        max_operon_gap_bp=max_operon_gap_bp,
        evalue_cutoff=evalue_cutoff,
        min_shell_length_aa=min_shell_length_aa,
        weights=weights or ScoringWeights(),
    )

    # 1. Load Genome Data
    logger.info("Step 1: Ingesting genome and gene annotations...")
    if gff_path is not None:
        gff_p = Path(gff_path)
        if not gff_p.is_file():
            raise FileNotFoundError(f"GFF3 annotation file not found: {gff_path}")
        genome = load_protein_and_gff3(in_p, gff_p, genome_id=genome_id)
    else:
        genome = call_genes_from_genome_fasta(in_p, genome_id=genome_id)

    logger.info(
        "Genome %s loaded: %d protein-coding genes across %d contig(s)",
        genome.genome_id,
        len(genome.genes),
        len(genome.contig_ids),
    )

    # 2. Prepare Shell HMM Profiles
    logger.info("Step 2: Locating encapsulin shell HMM profiles...")
    candidate_profile_dirs: List[Path] = []
    if profile_dir is not None:
        candidate_profile_dirs.append(Path(profile_dir))
    candidate_profile_dirs.extend(
        [
            Path("profiles"),
            out_p / "profiles",
            Path(__file__).parent.parent.parent / "profiles",
            Path(__file__).parent / "data",
        ]
    )

    profile_files: List[Path] = []
    for pdir in candidate_profile_dirs:
        if pdir.is_dir():
            found = list(pdir.glob("*.hmm")) + list(pdir.glob("*.hmm.gz"))
            if found:
                profile_files = sorted(found)
                logger.info("Found %d HMM profile(s) in %s", len(profile_files), pdir)
                break

    if not profile_files:
        logger.info(
            "No local profiles found; downloading default Linocin_M18 profile..."
        )
        dl_dir = out_p / "profiles"
        profile_files = fetch_default_profiles(dl_dir, DEFAULT_SHELL_PROFILES)

    # 3. Detect Shell Genes
    logger.info("Step 3: Scanning genes for encapsulin shell profiles...")
    shells = detect_shells(
        genome=genome,
        profile_paths=profile_files,
        evalue_cutoff=config.evalue_cutoff,
        min_length_aa=config.min_shell_length_aa,
    )

    output_files: Dict[str, Path] = {}

    if not shells:
        logger.warning(
            "No encapsulin shell genes found in genome %s (E-value <= %e).",
            genome.genome_id,
            config.evalue_cutoff,
        )
        # Generate empty outputs and summary
        tsv_path = export_candidates_tsv([], out_p / "candidates.tsv")
        json_path = export_candidates_json([], out_p / "candidates.json")
        summary = generate_run_summary(genome, [], [], config)
        summary_path = save_run_summary(summary, out_p / "run_summary.json")
        output_files["tsv"] = tsv_path
        output_files["json"] = json_path
        output_files["summary"] = summary_path
        return RunResult(genome, [], [], output_files, summary)

    # 4. Extract Neighborhood
    logger.info(
        "Step 4: Extracting genomic neighborhoods (+/- %d genes)...", config.window_size
    )
    neighbor_cands: List[NeighborCandidate] = extract_neighborhood(
        genome=genome,
        shell_hits=shells,
        window_size=config.window_size,
        max_operon_gap_bp=config.max_operon_gap_bp,
    )

    # 5. Cluster Candidate Sequences
    logger.info("Step 5: Clustering neighbor sequences via MMseqs2...")
    cluster_res = cluster_candidates(
        candidates=neighbor_cands,
        mmseqs_bin=mmseqs_bin,
    )

    # 6. Score and Rank
    logger.info("Step 6: Computing terminal features and scoring candidates...")
    motifs: Optional[Dict[str, str]] = None
    if motif_file is not None:
        motifs = load_motif_patterns(motif_file)

    scored = score_candidates(
        candidates=neighbor_cands,
        cluster_result=cluster_res,
        weights=config.weights,
        window_size=config.window_size,
        tail_length=config.tail_length,
        motifs=motifs,
    )

    # 7. Export Reports
    logger.info("Step 7: Writing output reports and diagrams...")
    output_files["tsv"] = export_candidates_tsv(scored, out_p / "candidates.tsv")
    output_files["json"] = export_candidates_json(scored, out_p / "candidates.json")

    summary = generate_run_summary(genome, shells, scored, config)
    output_files["summary"] = save_run_summary(summary, out_p / "run_summary.json")

    # 8. Generate Visualizations
    if generate_plots:
        plots_dir = out_p / "plots"
        plots_dir.mkdir(parents=True, exist_ok=True)
        for s in shells:
            plot_file = plots_dir / f"neighborhood_{s.gene_id}.png"
            plot_neighborhood(
                shell=s,
                genome=genome,
                candidates=scored,
                out_path=plot_file,
                flank_genes=config.window_size + 1,
            )
            output_files[f"plot_{s.gene_id}"] = plot_file

    logger.info(
        "Run completed successfully: %d candidate(s) ranked across %d shell(s)",
        len(scored),
        len(shells),
    )
    return RunResult(genome, shells, scored, output_files, summary)
