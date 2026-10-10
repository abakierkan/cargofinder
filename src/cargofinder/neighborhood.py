"""Genomic neighborhood extraction around encapsulin shell genes."""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import List

from cargofinder.io import Gene, GenomeData
from cargofinder.shells import ShellHit

logger = logging.getLogger(__name__)


@dataclass
class NeighborCandidate:
    """A gene located in the genomic vicinity (+/- N window) of an encapsulin shell."""

    candidate_id: str
    shell_gene_id: str
    contig_id: str
    gene_distance: int
    abs_gene_distance: int
    strand_relative_distance: int
    bp_distance: int
    same_strand: bool
    is_likely_operon: bool
    candidate_gene: Gene
    shell_hit: ShellHit

    @property
    def candidate_start(self) -> int:
        return self.candidate_gene.start

    @property
    def candidate_end(self) -> int:
        return self.candidate_gene.end

    @property
    def candidate_strand(self) -> int:
        return self.candidate_gene.strand

    @property
    def candidate_protein_seq(self) -> str:
        return self.candidate_gene.protein_seq


def _compute_intergenic_distance(g1: Gene, g2: Gene) -> int:
    """Compute nucleotide distance between two genes on the same contig.

    g1 is assumed to precede g2 in genomic start coordinates (g1.start <= g2.start).
    Returns positive value if there is a gap, negative if overlapping.
    """
    if g1.end < g2.start:
        return g2.start - g1.end - 1
    # Overlapping genes
    return -(g1.end - g2.start + 1)


def _check_operon_continuity(
    contig_genes: List[Gene],
    start_idx: int,
    end_idx: int,
    required_strand: int,
    max_gap_bp: int,
) -> bool:
    """Check if all genes between start_idx and end_idx form an unbroken operon.

    Conditions:
    1. All genes in the interval have strand == required_strand.
    2. Every pairwise intergenic gap between adjacent genes is <= max_gap_bp.
    """
    i1, i2 = min(start_idx, end_idx), max(start_idx, end_idx)
    for i in range(i1, i2 + 1):
        if contig_genes[i].strand != required_strand:
            return False

    for i in range(i1, i2):
        gap = _compute_intergenic_distance(contig_genes[i], contig_genes[i + 1])
        if gap > max_gap_bp:
            return False

    return True


def extract_neighborhood(
    genome: GenomeData,
    shell_hits: List[ShellHit],
    window_size: int = 5,
    max_operon_gap_bp: int = 150,
) -> List[NeighborCandidate]:
    """Extract candidate cargo genes in the genomic window around each shell hit.

    Parameters
    ----------
    genome : GenomeData
        Annotated genome containing all genes.
    shell_hits : list of ShellHit
        Identified encapsulin shell genes.
    window_size : int
        Number of genes upstream and downstream to include (default 5).
    max_operon_gap_bp : int
        Maximum intergenic distance in nucleotides between consecutive co-directional
        genes to consider them part of the same operon (default 150 bp).

    Returns
    -------
    list of NeighborCandidate
        Candidate cargo genes extracted across all shells in the genome.
    """
    if not shell_hits:
        logger.debug("No shell hits provided; returning empty neighborhood list")
        return []

    candidates: List[NeighborCandidate] = []

    for shell in shell_hits:
        contig_genes = genome.genes_on_contig(shell.contig_id)
        if not contig_genes:
            logger.warning(
                "Contig %s containing shell %s has no genes in genome record",
                shell.contig_id,
                shell.gene_id,
            )
            continue

        # Find the index of the shell gene on this contig
        shell_idx = -1
        for idx, g in enumerate(contig_genes):
            if g.gene_id == shell.gene_id:
                shell_idx = idx
                break

        if shell_idx == -1:
            logger.warning(
                "Shell gene %s not found on contig %s",
                shell.gene_id,
                shell.contig_id,
            )
            continue

        n_genes = len(contig_genes)
        min_idx = max(0, shell_idx - window_size)
        max_idx = min(n_genes - 1, shell_idx + window_size)

        for neighbor_idx in range(min_idx, max_idx + 1):
            if neighbor_idx == shell_idx:
                # Do not evaluate shell as its own candidate
                continue

            neighbor_gene = contig_genes[neighbor_idx]
            gene_distance = neighbor_idx - shell_idx
            abs_dist = abs(gene_distance)

            # Strand-relative distance:
            # Positive = downstream of shell transcription
            # Negative = upstream of shell transcription
            strand_rel_dist = gene_distance * shell.strand

            same_strand = neighbor_gene.strand == shell.strand

            # Direct physical intergenic distance to shell
            if neighbor_gene.start <= shell.gene.start:
                bp_dist = _compute_intergenic_distance(neighbor_gene, shell.gene)
            else:
                bp_dist = _compute_intergenic_distance(shell.gene, neighbor_gene)

            # Operon likelihood: same strand and uninterrupted by opposite
            # strand genes or large intergenic gaps
            is_likely_operon = False
            if same_strand:
                is_likely_operon = _check_operon_continuity(
                    contig_genes=contig_genes,
                    start_idx=shell_idx,
                    end_idx=neighbor_idx,
                    required_strand=shell.strand,
                    max_gap_bp=max_operon_gap_bp,
                )

            cand = NeighborCandidate(
                candidate_id=neighbor_gene.gene_id,
                shell_gene_id=shell.gene_id,
                contig_id=shell.contig_id,
                gene_distance=gene_distance,
                abs_gene_distance=abs_dist,
                strand_relative_distance=strand_rel_dist,
                bp_distance=bp_dist,
                same_strand=same_strand,
                is_likely_operon=is_likely_operon,
                candidate_gene=neighbor_gene,
                shell_hit=shell,
            )
            candidates.append(cand)

    logger.info(
        "Extracted %d neighborhood candidates for %d shell(s) (window=%d)",
        len(candidates),
        len(shell_hits),
        window_size,
    )
    return candidates
