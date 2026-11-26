"""Transparent additive scoring engine for candidate cargo discovery.

Biological and design principles:
1. No black-box machine learning in v0.1: Every candidate row must provide full
   breakdown into interpretable biological sub-scores.
2. Synteny/Adjacency: Distance decay reflects high operon linkage in prokaryotes.
3. Operon likelihood: Co-directional genes within short intergenic spacing (<=150 bp).
4. TP consensus similarity: BLOSUM62 similarity of the terminal slice to the
   published Family 1 targeting-peptide consensus SDGSLGIGSLKRS (Giessen & Silver
   2017 Nat Microbiol 2:17029; used as a blastp query to assign Family 1 cargo by
   Andreas & Giessen 2021 Nat Commun 12:4748).
5. Terminal peptide composition: Hydrophobic terminus, flexible Gly/Ser linker,
   optional motifs. This component is a cargofinder heuristic with no published
   basis and is down-weighted accordingly.
6. Cargo class identity: whether the candidate matches a Pfam family marking a
   known encapsulin cargo class. Andreas & Giessen (2021) assigned cargo this
   way. See cargofinder.cargo for what a match does and does not mean.
7. Cross-system recurrence: Homologous candidate families found next to independent
   shells across loci/genomes receive an evolutionary conservation bonus.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

from cargofinder.cargo import CargoClassHit
from cargofinder.cluster import ClusterResult
from cargofinder.config import (
    NTERM_EXTENSION_FULL_AA,
    NTERM_EXTENSION_MIN_AA,
    TARGETING_MODE_BY_PROFILE,
    TARGETING_MODE_C_TERMINAL,
    TARGETING_MODE_DEFAULT,
    TARGETING_MODE_N_TERMINAL,
    ScoringWeights,
)
from cargofinder.features import PeptideFeatures, compute_peptide_features
from cargofinder.neighborhood import NeighborCandidate

logger = logging.getLogger(__name__)


@dataclass
class ScoredCandidate:
    """A ranked candidate cargo gene with full score component breakdown."""

    candidate_id: str
    shell_gene_id: str
    contig_id: str
    rank: int
    score_total: float
    score_distance: float
    score_operon: float
    score_targeting: float
    score_cargo_class: float
    score_tail_features: float
    score_recurrence: float
    gene_distance: int
    abs_gene_distance: int
    strand_relative_distance: int
    bp_distance: int
    same_strand: bool
    is_likely_operon: bool
    cluster_id: str
    cluster_size: int
    cluster_shell_count: int
    candidate: NeighborCandidate
    features: PeptideFeatures
    rank_in_shell: int = 0
    best_terminus: str = "c"
    # Prophage context is evidence about the shell locus, not about this
    # candidate being cargo, so it is reported but never added to score_total.
    # Populated by the pipeline after scoring; see cargofinder.prophage.
    targeting_mode: str = TARGETING_MODE_C_TERMINAL
    n_terminal_extension_aa: int = 0
    cargo_class: str = ""
    cargo_class_pfam: str = ""
    cargo_class_specific: bool = False
    prophage_marker_count: int = 0
    prophage_marker_families: int = 0
    is_probable_prophage_locus: bool = False

    def to_dict(self) -> Dict[str, object]:
        """Convert scored candidate to a flat dictionary for TSV/JSON reporting."""
        d: Dict[str, object] = {
            "rank": self.rank,
            "rank_in_shell": self.rank_in_shell,
            "candidate_id": self.candidate_id,
            "shell_gene_id": self.shell_gene_id,
            "genome_id": self.candidate.genome_id,
            "contig_id": self.contig_id,
            "shell_profile": self.candidate.shell_hit.profile_name,
            "shell_evalue": self.candidate.shell_hit.evalue,
            "shell_bitscore": round(self.candidate.shell_hit.bitscore, 2),
            "score_total": round(self.score_total, 4),
            "score_distance": round(self.score_distance, 4),
            "score_operon": round(self.score_operon, 4),
            "score_targeting": round(self.score_targeting, 4),
            "targeting_mode": self.targeting_mode,
            "n_terminal_extension_aa": self.n_terminal_extension_aa,
            "score_cargo_class": round(self.score_cargo_class, 4),
            "cargo_class": self.cargo_class,
            "cargo_class_pfam": self.cargo_class_pfam,
            "cargo_class_specific": self.cargo_class_specific,
            "score_tail_features": round(self.score_tail_features, 4),
            "score_recurrence": round(self.score_recurrence, 4),
            "best_terminus": self.best_terminus,
            "gene_distance": self.gene_distance,
            "abs_gene_distance": self.abs_gene_distance,
            "strand_relative_distance": self.strand_relative_distance,
            "bp_distance": self.bp_distance,
            "same_strand": self.same_strand,
            "is_likely_operon": self.is_likely_operon,
            "candidate_start": self.candidate.candidate_start,
            "candidate_end": self.candidate.candidate_end,
            "candidate_strand": self.candidate.candidate_strand,
            "candidate_length_aa": self.candidate.candidate_gene.length_aa,
            "candidate_product": self.candidate.candidate_gene.product or "",
            "cluster_id": self.cluster_id,
            "cluster_size": self.cluster_size,
            "cluster_shell_count": self.cluster_shell_count,
            "prophage_marker_count": self.prophage_marker_count,
            "prophage_marker_families": self.prophage_marker_families,
            "is_probable_prophage_locus": self.is_probable_prophage_locus,
        }
        # Add peptide features
        d.update(self.features.to_dict())
        return d


def compute_distance_score(abs_distance: int, window_size: int = 5) -> float:
    """Linear distance decay score in [0.0, 1.0].

    abs_distance = 1 -> 1.0
    abs_distance = window_size -> 1.0 / window_size
    abs_distance > window_size -> 0.0
    """
    if abs_distance <= 0:
        return 1.0
    if abs_distance > window_size:
        return 0.0
    return (window_size - abs_distance + 1) / window_size


def compute_operon_score(same_strand: bool, is_likely_operon: bool) -> float:
    """Operon continuity score in [0.0, 1.0].

    - 1.0: co-directional and unbroken tight operon gap (<= 150 bp)
    - 0.5: same strand but broken operon or larger gap
    - 0.0: opposite strand
    """
    if is_likely_operon:
        return 1.0
    if same_strand:
        return 0.5
    return 0.0


def compute_tail_features_score(
    features: PeptideFeatures,
    has_motifs: bool = False,
) -> float:
    """Normalized C-terminal peptide score in [0.0, 1.0].

    Components:
    - Hydrophobic terminus fraction: weight 0.60 (or 0.48 with motifs)
    - Flexible Gly/Ser linker: weight 0.40 (or 0.32 with motifs)
    - Motif match: weight 0.20 (if motifs provided)
    """
    base = _terminus_base_score(
        features.c_hydrophobic_fraction, features.c_gly_ser_fraction
    )
    if has_motifs and features.motif_matches:
        s_motif = 1.0 if any(features.motif_matches.values()) else 0.0
        return 0.80 * base + 0.20 * s_motif
    return base


def _terminus_base_score(hydrophobic_fraction: float, gly_ser_fraction: float) -> float:
    """Hydrophobic and Gly/Ser score; 35% Gly/Ser achieves the maximum of 1.0."""
    return 0.60 * hydrophobic_fraction + 0.40 * min(1.0, gly_ser_fraction / 0.35)


def score_termini(
    features: PeptideFeatures,
    has_motifs: bool = False,
    termini: str = "c",
) -> Tuple[float, str]:
    """Return the tail score and the terminus it came from ('c' or 'n').

    With termini='both' the better terminus counts. Motifs are matched against the
    C-terminus only. N-termini need care: signal peptides are also hydrophobic.
    """
    c_score = compute_tail_features_score(features, has_motifs=has_motifs)
    if termini == "c":
        return c_score, "c"
    n_score = _terminus_base_score(
        features.n_hydrophobic_fraction, features.n_gly_ser_fraction
    )
    return (n_score, "n") if n_score > c_score else (c_score, "c")


def targeting_mode_for_profile(profile_acc: str) -> str:
    """Which targeting mechanism applies at a shell locus, from its Pfam family.

    Family 1 recruits cargo with a short C-terminal targeting peptide; Family 2
    uses a large N-terminal domain instead, and has no short targeting peptide
    at all (Andreas & Giessen 2021 Nat Commun 12:4748). Scoring a Family 2
    candidate against the Family 1 C-terminal consensus therefore measures
    nothing. Unknown families fall back to the C-terminal mode, which is a
    documented default rather than a claim; every row records the mode used.
    """
    bare = (profile_acc or "").split(".")[0]
    return TARGETING_MODE_BY_PROFILE.get(bare, TARGETING_MODE_DEFAULT)


def score_c_terminal_tp(features: PeptideFeatures, termini: str = "c") -> float:
    """TP consensus similarity for the scored terminus/termini."""
    if termini == "c":
        return features.c_tp_consensus_similarity
    return max(
        features.c_tp_consensus_similarity,
        features.n_tp_consensus_similarity,
    )


def score_n_terminal_extension(extension_aa: int) -> float:
    """Score an unannotated N-terminal extension, the Family 2 targeting signal.

    Ramps linearly from 0 at NTERM_EXTENSION_MIN_AA to 1 at
    NTERM_EXTENSION_FULL_AA. Those bounds come from the truncation series of
    Nichols et al. (2021) eLife 10:e59288 on the S. elongatus CyD cargo: the
    full 225-residue N-terminal domain loaded best, the first 100 residues still
    worked but less well, a 29-residue fragment carrying a conserved motif did
    not work, and deleting the domain abolished loading. See config for the
    quoted sentences. The linear shape between the anchors is a judgment call.

    Returns 0.0 when the extension cannot be measured, which is the same value
    as "no extension". That is a known ambiguity: it requires a recognised
    domain to measure the extension against, so novel cargo of an unknown class
    gets nothing here.
    """
    if extension_aa <= NTERM_EXTENSION_MIN_AA:
        return 0.0
    span = NTERM_EXTENSION_FULL_AA - NTERM_EXTENSION_MIN_AA
    if span <= 0:
        return 1.0
    return min(1.0, (extension_aa - NTERM_EXTENSION_MIN_AA) / span)


def compute_recurrence_score(shell_count: int) -> float:
    """Cluster recurrence score in [0.0, 1.0].

    - 1 shell: 0.2 (baseline single locus)
    - 2 shells: 0.6
    - >= 3 shells: 1.0
    """
    if shell_count <= 1:
        return 0.2
    return min(1.0, round(0.2 + 0.4 * (shell_count - 1), 6))


def score_candidates(
    candidates: List[NeighborCandidate],
    cluster_result: ClusterResult,
    weights: Optional[ScoringWeights] = None,
    window_size: int = 5,
    tail_length: int = 20,
    motifs: Optional[Dict[str, str]] = None,
    termini: str = "c",
    cargo_classes: Optional[Dict[str, CargoClassHit]] = None,
) -> List[ScoredCandidate]:
    """Score and rank all neighbor candidates.

    Parameters
    ----------
    candidates : list of NeighborCandidate
        Candidate genes around encapsulin shells.
    cluster_result : ClusterResult
        MMseqs2 clustering result.
    weights : ScoringWeights, optional
        Weights for additive scoring model. Defaults to standard weights.
    window_size : int
        Window size used for distance normalization.
    tail_length : int
        Tail length for peptide feature extraction.
    motifs : dict of str -> str, optional
        Targeting peptide regex motifs.
    termini : str
        'c' scores the C-terminus only; 'both' uses the better of the two termini.
    cargo_classes : dict of str -> CargoClassHit, optional
        Known cargo class calls keyed by candidate_id. When omitted, the cargo
        class component scores 0 for every candidate; the caller is expected to
        pass weights from `ScoringWeights.without_cargo_class()` in that case so
        totals stay on the same scale.

    Returns
    -------
    list of ScoredCandidate
        Ranked candidates sorted descending by score_total.
    """
    if weights is None:
        weights = ScoringWeights()
    weights.validate()

    has_motifs = bool(motifs)
    scored: List[ScoredCandidate] = []

    for cand in candidates:
        cid = cand.candidate_id
        cluster_id = cluster_result.get_cluster_id(cid)
        cluster_size = cluster_result.get_cluster_size(cid)
        shell_count = cluster_result.get_shell_recurrence(cid)

        # Compute peptide features
        feats = compute_peptide_features(
            cand.candidate_protein_seq,
            tail_length=tail_length,
            motifs=motifs,
        )

        s_dist = compute_distance_score(cand.abs_gene_distance, window_size=window_size)
        s_op = compute_operon_score(cand.same_strand, cand.is_likely_operon)
        s_tail, best_terminus = score_termini(
            feats, has_motifs=has_motifs, termini=termini
        )
        cargo_hit = (cargo_classes or {}).get(cid)
        s_cargo = cargo_hit.score if cargo_hit is not None else 0.0

        mode = targeting_mode_for_profile(cand.shell_hit.profile_acc)
        nterm_ext = cargo_hit.n_terminal_extension_aa if cargo_hit else 0
        if mode == TARGETING_MODE_N_TERMINAL:
            s_target = score_n_terminal_extension(nterm_ext)
        else:
            s_target = score_c_terminal_tp(feats, termini=termini)
        s_rec = compute_recurrence_score(shell_count)

        total = (
            weights.w_distance * s_dist
            + weights.w_operon * s_op
            + weights.w_targeting * s_target
            + weights.w_cargo_class * s_cargo
            + weights.w_tail_features * s_tail
            + weights.w_recurrence * s_rec
        )

        scored.append(
            ScoredCandidate(
                candidate_id=cid,
                shell_gene_id=cand.shell_gene_id,
                contig_id=cand.contig_id,
                rank=0,  # assigned after sorting
                score_total=total,
                score_distance=s_dist,
                score_operon=s_op,
                score_targeting=s_target,
                targeting_mode=mode,
                n_terminal_extension_aa=nterm_ext,
                score_cargo_class=s_cargo,
                cargo_class=cargo_hit.cargo_class if cargo_hit else "",
                cargo_class_pfam=cargo_hit.pfam_id if cargo_hit else "",
                cargo_class_specific=(
                    cargo_hit.encapsulin_specific if cargo_hit else False
                ),
                score_tail_features=s_tail,
                score_recurrence=s_rec,
                gene_distance=cand.gene_distance,
                abs_gene_distance=cand.abs_gene_distance,
                strand_relative_distance=cand.strand_relative_distance,
                bp_distance=cand.bp_distance,
                same_strand=cand.same_strand,
                is_likely_operon=cand.is_likely_operon,
                cluster_id=cluster_id,
                cluster_size=cluster_size,
                cluster_shell_count=shell_count,
                candidate=cand,
                features=feats,
                best_terminus=best_terminus,
            )
        )

    # Sort descending by score_total (with ties broken by distance, then candidate_id)
    scored.sort(key=lambda x: (-x.score_total, x.abs_gene_distance, x.candidate_id))

    # Assign 1-based ranks, across the run and within each shell's neighborhood
    per_shell: Dict[str, int] = {}
    for rank_idx, item in enumerate(scored, start=1):
        item.rank = rank_idx
        per_shell[item.shell_gene_id] = per_shell.get(item.shell_gene_id, 0) + 1
        item.rank_in_shell = per_shell[item.shell_gene_id]

    logger.info("Scored and ranked %d candidate cargo genes", len(scored))
    return scored
