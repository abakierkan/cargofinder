"""High-level Python API for encapsulin cargo candidate discovery."""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence

from cargofinder.cargo import (
    CargoClassHit,
    classify_cargo_candidates,
    find_cargo_class_profiles,
)
from cargofinder.cluster import (
    cluster_candidates,
    get_mmseqs_version,
    singleton_clusters,
)
from cargofinder.config import DEFAULT_SHELL_PROFILES, RunConfig, ScoringWeights
from cargofinder.features import load_motif_patterns
from cargofinder.genecalling import call_genes_from_genome_fasta
from cargofinder.io import Gene, GenomeData, load_protein_and_gff3, namespace_genome
from cargofinder.neighborhood import NeighborCandidate, extract_neighborhood
from cargofinder.plotting import plot_neighborhood
from cargofinder.profiles import (
    describe_profiles,
    fetch_default_profiles,
    find_profile_files,
)
from cargofinder.prophage import (
    ProphageEvidence,
    assess_prophage_context,
    collect_prophage_windows,
    find_phage_marker_profiles,
)
from cargofinder.reporting import (
    export_candidates_json,
    export_candidates_tsv,
    export_shells_tsv,
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
    genomes: Optional[List[GenomeData]] = None
    prophage: Dict[str, ProphageEvidence] = field(default_factory=dict)
    cargo_classes: Dict[str, CargoClassHit] = field(default_factory=dict)


def resolve_profile_files(
    profile_dir: Optional[Path | str], out_dir: Path
) -> List[Path]:
    """Find the shell profiles for a run.

    An explicit directory must contain profiles. Without one, ./profiles is used if
    it has any, and otherwise the default profile is downloaded into the output
    directory.
    """
    if profile_dir is not None:
        return find_profile_files(profile_dir)

    local = Path("profiles")
    if local.is_dir() and (list(local.glob("*.hmm")) + list(local.glob("*.hmm.gz"))):
        logger.info("Using HMM profiles found in %s", local)
        return find_profile_files(local)

    logger.info("No local profiles found; downloading default Linocin_M18 profile...")
    return fetch_default_profiles(out_dir / "profiles", DEFAULT_SHELL_PROFILES)


def _profile_root(profile_dir: Optional[Path | str], out_dir: Path) -> Path:
    """Return the directory that holds this run's profiles.

    Mirrors the search order of `resolve_profile_files` so that prophage markers
    are looked for beside the shell profiles actually used.
    """
    if profile_dir is not None:
        return Path(profile_dir)
    local = Path("profiles")
    if local.is_dir() and (list(local.glob("*.hmm")) + list(local.glob("*.hmm.gz"))):
        return local
    return out_dir / "profiles"


def _load_genomes(
    inputs: List[Path], gff_path: Optional[Path | str], genome_id: Optional[str]
) -> List[GenomeData]:
    """Load every input genome, namespacing IDs when there is more than one."""
    if gff_path is not None and len(inputs) > 1:
        raise ValueError(
            "--gff supports one protein FASTA per run. "
            "Pass several nucleotide FASTA files, or run annotated genomes separately."
        )
    if genome_id is not None and len(inputs) > 1:
        raise ValueError("genome_id can only be set for a single-genome run.")

    genomes: List[GenomeData] = []
    for in_p in inputs:
        if not in_p.is_file():
            raise FileNotFoundError(f"Input sequence file not found: {in_p}")
        if gff_path is not None:
            gff_p = Path(gff_path)
            if not gff_p.is_file():
                raise FileNotFoundError(f"GFF3 annotation file not found: {gff_path}")
            genomes.append(load_protein_and_gff3(in_p, gff_p, genome_id=genome_id))
        else:
            genomes.append(call_genes_from_genome_fasta(in_p, genome_id=genome_id))

    if len(genomes) > 1:
        ids = [g.genome_id for g in genomes]
        duplicates = sorted({i for i in ids if ids.count(i) > 1})
        if duplicates:
            raise ValueError(
                f"Input files share the genome name(s) {duplicates}; genome names "
                "come from file names, so rename the files to be unique."
            )
        genomes = [namespace_genome(g) for g in genomes]
    return genomes


def safe_filename(name: str) -> str:
    """Replace characters that are not portable in file names."""
    return re.sub(r"[^A-Za-z0-9._-]", "_", name)


def run_pipeline(
    input_path: Path | str | Sequence[Path | str],
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
    termini: str = "c",
    drop_prophage_loci: bool = False,
    prophage_window: int = 20,
) -> RunResult:
    """Execute the full cargofinder discovery pipeline.

    Parameters
    ----------
    input_path : Path or str, or a sequence of them
        Nucleotide FASTA file(s), or one protein FASTA if gff_path is provided.
        With several genomes, gene and contig IDs are prefixed with "<genome>|" and
        cluster recurrence is measured across them.
    gff_path : Path or str, optional
        Path to GFF3 annotation file.
    out_dir : Path or str
        Output directory for reports, tables, and figures.
    profile_dir : Path or str, optional
        Directory containing encapsulin shell .hmm profiles. Must contain profiles
        if given. If None, ./profiles is used if present, otherwise the default
        profile is downloaded into out_dir/profiles.
    motif_file : Path or str, optional
        JSON file with targeting peptide regex patterns.
    window_size : int
        Number of genes to search upstream and downstream of each shell (default 5).
    tail_length : int
        Residue count for terminal peptide feature extraction (default 20).
    max_operon_gap_bp : int
        Max intergenic gap for likely operon co-transcription (default 150).
    evalue_cutoff : float
        Full-sequence HMMER E-value threshold for shell detection (default 1e-5).
    min_shell_length_aa : int
        Minimum length in amino acids for encapsulin shell monomers (default 150).
    weights : ScoringWeights, optional
        Scoring weights. Defaults to standard weights.
    genome_id : str, optional
        Explicit genome identifier (single-genome runs only).
    mmseqs_bin : Path or str, optional
        Explicit path to mmseqs binary.
    generate_plots : bool
        Whether to generate neighborhood diagrams (default True).
    termini : str
        'c' scores C-terminal tails only (default); 'both' uses the better terminus.

    Returns
    -------
    RunResult
        Container with all results, candidate rankings, and output file paths.
    """
    if isinstance(input_path, (str, Path)):
        inputs = [Path(input_path)]
    else:
        inputs = [Path(p) for p in input_path]
    if not inputs:
        raise ValueError("No input files provided.")

    config = RunConfig(
        window_size=window_size,
        tail_length=tail_length,
        max_operon_gap_bp=max_operon_gap_bp,
        evalue_cutoff=evalue_cutoff,
        min_shell_length_aa=min_shell_length_aa,
        weights=weights or ScoringWeights(),
        termini=termini,
    )
    config.validate()

    for in_p in inputs:
        if not in_p.is_file():
            raise FileNotFoundError(f"Input sequence file not found: {in_p}")

    out_p = Path(out_dir)
    out_p.mkdir(parents=True, exist_ok=True)

    # 1. Load Genome Data
    logger.info("Step 1: Ingesting genome and gene annotations...")
    genomes = _load_genomes(inputs, gff_path, genome_id)
    for genome in genomes:
        logger.info(
            "Genome %s loaded: %d protein-coding genes across %d contig(s)",
            genome.genome_id,
            len(genome.genes),
            len(genome.contig_ids),
        )

    # 2. Prepare Shell HMM Profiles
    logger.info("Step 2: Locating encapsulin shell HMM profiles...")
    profile_files = resolve_profile_files(profile_dir, out_p)
    profile_info = describe_profiles(profile_files)

    # 3. Detect Shell Genes
    logger.info("Step 3: Scanning genes for encapsulin shell profiles...")
    shells_by_genome = [
        detect_shells(
            genome=genome,
            profile_paths=profile_files,
            evalue_cutoff=config.evalue_cutoff,
            min_length_aa=config.min_shell_length_aa,
        )
        for genome in genomes
    ]
    shells = [s for group in shells_by_genome for s in group]

    output_files: Dict[str, Path] = {}
    output_files["shells"] = export_shells_tsv(shells, out_p / "shells.tsv")

    if not shells:
        logger.warning(
            "No encapsulin shell genes found (E-value <= %e). Profiles used: %s.",
            config.evalue_cutoff,
            ", ".join(p["acc"] or p["name"] for p in profile_info),
        )
        output_files["tsv"] = export_candidates_tsv([], out_p / "candidates.tsv")
        output_files["json"] = export_candidates_json([], out_p / "candidates.json")
        summary = generate_run_summary(
            genomes, [], [], config, profiles=profile_info, clustering="not run"
        )
        output_files["summary"] = save_run_summary(summary, out_p / "run_summary.json")
        return RunResult(genomes[0], [], [], output_files, summary, genomes)

    # 4. Extract Neighborhood
    logger.info(
        "Step 4: Extracting genomic neighborhoods (+/- %d genes)...", config.window_size
    )
    neighbor_cands: List[NeighborCandidate] = []
    for genome, group in zip(genomes, shells_by_genome, strict=True):
        neighbor_cands.extend(
            extract_neighborhood(
                genome=genome,
                shell_hits=group,
                window_size=config.window_size,
                max_operon_gap_bp=config.max_operon_gap_bp,
            )
        )

    # 4b. Screen each shell locus for prophage context.
    #
    # Encapsulin shells share the HK97 fold with phage capsids, and
    # Andreas & Giessen (2021) found their Family 1 and 2 datasets significantly
    # contaminated with phage capsid proteins. See cargofinder.prophage.
    marker_profiles = find_phage_marker_profiles(_profile_root(profile_dir, out_p))
    prophage_windows = collect_prophage_windows(
        genomes, shells_by_genome, prophage_window
    )
    prophage_evidence = assess_prophage_context(
        prophage_windows,
        marker_profiles,
        evalue_cutoff=config.evalue_cutoff,
    )
    # Rewrite shells.tsv now that the prophage columns can be filled in.
    output_files["shells"] = export_shells_tsv(
        shells, out_p / "shells.tsv", prophage_evidence
    )

    if drop_prophage_loci:
        flagged = {
            shell_id
            for shell_id, ev in prophage_evidence.items()
            if ev.is_probable_prophage
        }
        if flagged:
            before = len(neighbor_cands)
            neighbor_cands = [
                c for c in neighbor_cands if c.shell_gene_id not in flagged
            ]
            shells = [s for s in shells if s.gene_id not in flagged]
            shells_by_genome = [
                [s for s in group if s.gene_id not in flagged]
                for group in shells_by_genome
            ]
            logger.warning(
                "--drop-prophage-loci removed %d shell locus/loci and %d candidate(s)",
                len(flagged),
                before - len(neighbor_cands),
            )
            output_files["shells"] = export_shells_tsv(
                shells, out_p / "shells.tsv", prophage_evidence
            )

    # 4c. Recognise known cargo classes among the candidates.
    #
    # Andreas & Giessen (2021) assigned cargo by identity, not position alone.
    # See cargofinder.cargo for what a class match does and does not mean.
    cargo_profiles = find_cargo_class_profiles(_profile_root(profile_dir, out_p))
    scoring_weights = config.weights
    cargo_classes: Dict[str, CargoClassHit] = {}
    if cargo_profiles:
        by_candidate_id = {c.candidate_id: c.candidate_gene for c in neighbor_cands}
        # Classify each distinct gene once, then map back to candidate ids, since
        # one gene can be a candidate for two overlapping shells.
        unique: Dict[str, Gene] = {}
        for cand in neighbor_cands:
            unique.setdefault(cand.candidate_gene.gene_id, cand.candidate_gene)
        by_gene = classify_cargo_candidates(
            list(unique.values()),
            cargo_profiles,
            evalue_cutoff=config.evalue_cutoff,
        )
        for cand_id, gene in by_candidate_id.items():
            hit = by_gene.get(gene.gene_id)
            if hit is not None:
                cargo_classes[cand_id] = hit
    else:
        # No profiles: drop the term and rescale, so totals stay on [0, 1] and
        # remain comparable with runs that did have the profiles.
        scoring_weights = config.weights.without_cargo_class()

    # 5. Cluster Candidate Sequences
    if len(shells) >= 2:
        logger.info("Step 5: Clustering neighbor sequences via MMseqs2...")
        cluster_res = cluster_candidates(
            candidates=neighbor_cands,
            min_seq_id=config.mmseqs_min_seq_id,
            coverage=config.mmseqs_coverage,
            mmseqs_bin=mmseqs_bin,
        )
        clustering = f"MMseqs2 across {len(shells)} shells"
        mmseqs_version = get_mmseqs_version(mmseqs_bin)
    else:
        logger.info(
            "Step 5: One shell in this run, so recurrence cannot be measured; "
            "skipping MMseqs2. Every candidate gets the baseline recurrence score."
        )
        cluster_res = singleton_clusters(neighbor_cands)
        clustering = "skipped (one shell: recurrence is not informative)"
        mmseqs_version = "not used"

    # 6. Score and Rank
    logger.info("Step 6: Computing terminal features and scoring candidates...")
    motifs: Optional[Dict[str, str]] = None
    if motif_file is not None:
        motifs = load_motif_patterns(motif_file)

    scored = score_candidates(
        candidates=neighbor_cands,
        cluster_result=cluster_res,
        weights=scoring_weights,
        window_size=config.window_size,
        tail_length=config.tail_length,
        motifs=motifs,
        termini=config.termini,
        cargo_classes=cargo_classes,
    )

    for item in scored:
        ev = prophage_evidence.get(item.shell_gene_id)
        if ev is not None:
            item.prophage_marker_count = ev.marker_count
            item.prophage_marker_families = ev.distinct_family_count
            item.is_probable_prophage_locus = ev.is_probable_prophage

    # 7. Export Reports
    logger.info("Step 7: Writing output reports and diagrams...")
    output_files["tsv"] = export_candidates_tsv(scored, out_p / "candidates.tsv")
    output_files["json"] = export_candidates_json(scored, out_p / "candidates.json")

    summary = generate_run_summary(
        genomes,
        shells,
        scored,
        config,
        profiles=profile_info,
        mmseqs_version=mmseqs_version,
        clustering=clustering,
    )
    output_files["summary"] = save_run_summary(summary, out_p / "run_summary.json")

    # 8. Generate Visualizations
    if generate_plots:
        plots_dir = out_p / "plots"
        plots_dir.mkdir(parents=True, exist_ok=True)
        for genome, group in zip(genomes, shells_by_genome, strict=True):
            for s in group:
                plot_file = plots_dir / f"neighborhood_{safe_filename(s.gene_id)}.png"
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
    return RunResult(
        genomes[0],
        shells,
        scored,
        output_files,
        summary,
        genomes,
        prophage_evidence,
        cargo_classes,
    )
