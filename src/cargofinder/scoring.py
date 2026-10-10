"""Transparent additive scoring engine for candidate cargo discovery.

Biological and design principles:
1. No black-box machine learning in v0.1: Every candidate row must provide full
   breakdown into interpretable biological sub-scores.
2. Synteny/Adjacency: Distance decay reflects high operon linkage in prokaryotes.
3. Operon likelihood: Co-directional genes within short intergenic spacing (<=150 bp).
4. Terminal peptide: Hydrophobic terminus, flexible Gly/Ser linker, low complexity.
5. Cross-system recurrence: Homologous candidate families found next to independent
   shells across loci/genomes receive an evolutionary conservation bonus.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Dict, List, Optional

from cargofinder.cluster import ClusterResult
from cargofinder.config import ScoringWeights
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

    def to_dict(self) -> Dict[str, object]:
        """Convert scored candidate to a flat dictionary for TSV/JSON reporting."""
        d: Dict[str, object] = {
            "rank": self.rank,
            "candidate_id": self.candidate_id,
            "shell_gene_id": self.shell_gene_id,
            "contig_id": self.contig_id,
            "score_total": round(self.score_total, 4),
            "score_distance": round(self.score_distance, 4),
            "score_operon": round(self.score_operon, 4),
            "score_tail_features": round(self.score_tail_features, 4),
            "score_recurrence": round(self.score_recurrence, 4),
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
    """Normalized terminal peptide score in [0.0, 1.0].

    Components:
    - Hydrophobic terminus fraction: weight 0.40 (or 0.50 if no motifs)
    - Flexible Gly/Ser linker: weight 0.30 (or 0.35 if no motifs)
    - Low complexity score: weight 0.10 (or 0.15 if no motifs)
    - Motif match bonus: weight 0.20 (if motifs provided)
    """
    # Scale Gly/Ser: 35% Gly/Ser achieves maximum 1.0
    s_gly_ser = min(1.0, features.c_gly_ser_fraction / 0.35)
    s_hydro = features.c_hydrophobic_fraction
    s_low_comp = features.c_low_complexity_score

    if has_motifs and features.motif_matches:
        any_match = any(features.motif_matches.values())
        s_motif = 1.0 if any_match else 0.0
        score = 0.40 * s_hydro + 0.30 * s_gly_ser + 0.10 * s_low_comp + 0.20 * s_motif
    else:
        score = 0.50 * s_hydro + 0.35 * s_gly_ser + 0.15 * s_low_comp

    return max(0.0, min(1.0, score))


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
        s_tail = compute_tail_features_score(feats, has_motifs=has_motifs)
        s_rec = compute_recurrence_score(shell_count)

        total = (
            weights.w_distance * s_dist
            + weights.w_operon * s_op
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
            )
        )

    # Sort descending by score_total (with ties broken by distance, then candidate_id)
    scored.sort(key=lambda x: (-x.score_total, x.abs_gene_distance, x.candidate_id))

    # Assign 1-based ranks
    for rank_idx, item in enumerate(scored, start=1):
        item.rank = rank_idx

    logger.info("Scored and ranked %d candidate cargo genes", len(scored))
    return scored
