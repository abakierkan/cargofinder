"""Recognise prophage context around a shell hit.

Why this module exists
----------------------
Encapsulin shell proteins share the HK97 fold with bacteriophage major capsid
proteins; both sit in Pfam clan CL0373. A profile search for an encapsulin shell
therefore also finds phage capsids, and no E-value threshold separates them,
because the similarity is real.

Andreas & Giessen (2021) Nat Commun 12:4748, doi:10.1038/s41467-021-25071-y,
hit exactly this problem at scale: their Family 1 and Family 2 datasets "were
significantly contaminated with bacteriophage capsid proteins". They resolved it
with context rather than sequence, excluding hits whose genome neighborhoods
contained phage-associated Pfam domains.

This module does the same thing locally: it scans the genes already extracted
around each shell hit against a small set of phage marker profiles (terminase,
prohead protease, integrase, head-tail connector) and reports what it finds.

What it deliberately does not do
--------------------------------
It does not change any candidate's score, and it does not delete anything by
default. A prophage flag is evidence about the *shell* — whether this locus is an
encapsulin at all — not evidence about whether a given neighbour is cargo. Mixing
it into the additive cargo score would conflate the two. The flag is reported in
`shells.tsv` and on every candidate row, and `--drop-prophage-loci` removes
flagged loci if you want that.

The markers are necessary-but-not-sufficient in both directions: a real
encapsulin can sit near a genuine prophage, and a prophage can lack all four
markers within the window. Treat the flag as a reason to look, not a verdict.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Dict, Iterable, List, Sequence

import pyhmmer
from pyhmmer import easel

from cargofinder.config import PHAGE_MARKER_FLAG_THRESHOLD, PHAGE_MARKER_SUBDIR
from cargofinder.io import Gene
from cargofinder.shells import _as_text, load_hmm_profiles

if TYPE_CHECKING:
    from cargofinder.io import GenomeData
    from cargofinder.shells import ShellHit

logger = logging.getLogger(__name__)


@dataclass
class ProphageEvidence:
    """Phage marker genes found in one shell's neighborhood."""

    shell_gene_id: str
    marker_families: Dict[str, str] = field(default_factory=dict)
    """Maps gene_id -> phage marker profile name that it matched."""

    @property
    def marker_count(self) -> int:
        """Number of neighbouring genes that matched any phage marker."""
        return len(self.marker_families)

    @property
    def distinct_family_count(self) -> int:
        """Number of distinct phage marker families represented."""
        return len(set(self.marker_families.values()))

    @property
    def is_probable_prophage(self) -> bool:
        """True when at least PHAGE_MARKER_FLAG_THRESHOLD distinct families hit."""
        return self.distinct_family_count >= PHAGE_MARKER_FLAG_THRESHOLD

    @property
    def marker_summary(self) -> str:
        """Deterministic, human-readable summary for a table cell."""
        if not self.marker_families:
            return ""
        pairs = sorted(self.marker_families.items())
        return ";".join(f"{gene}:{fam}" for gene, fam in pairs)


def find_phage_marker_profiles(profile_dir: Path | str) -> List[Path]:
    """Locate phage marker profiles in the `phage_markers` subdirectory.

    Returns an empty list when the subdirectory is absent, so that a profile
    directory fetched before this feature existed keeps working: the run simply
    reports no prophage screening rather than failing.
    """
    sub = Path(profile_dir) / PHAGE_MARKER_SUBDIR
    if not sub.is_dir():
        logger.warning(
            "No %s/ subdirectory in %s, so prophage screening is skipped. "
            "Run 'cargofinder fetch-profiles' again to download the phage "
            "marker profiles.",
            PHAGE_MARKER_SUBDIR,
            profile_dir,
        )
        return []
    return sorted(list(sub.glob("*.hmm")) + list(sub.glob("*.hmm.gz")))


def scan_genes_for_phage_markers(
    genes: Sequence[Gene],
    marker_profile_paths: Iterable[Path | str],
    evalue_cutoff: float = 1e-5,
) -> Dict[str, str]:
    """Scan proteins against phage marker profiles.

    Parameters
    ----------
    genes : sequence of Gene
        Genes to scan, typically the union of all shell neighborhood windows.
    marker_profile_paths : iterable of Path or str
        Phage marker HMM files.
    evalue_cutoff : float
        Maximum full-sequence E-value for a marker call.

    Returns
    -------
    dict of str -> str
        Maps gene_id to the name of the best-scoring phage marker profile.
        Genes with no marker hit are absent from the mapping.
    """
    paths = [Path(p) for p in marker_profile_paths]
    if not paths or not genes:
        return {}

    hmms = load_hmm_profiles(paths)
    alphabet = easel.Alphabet.amino()

    digital_seqs: List[easel.DigitalSequence] = []
    for g in genes:
        if not g.protein_seq:
            continue
        text_seq = easel.TextSequence(
            name=g.gene_id.encode("utf-8"),
            sequence=g.protein_seq,
        )
        digital_seqs.append(text_seq.digitize(alphabet))

    if not digital_seqs:
        return {}

    # Keep the best-scoring marker per gene so the result is deterministic
    # regardless of profile iteration order.
    best: Dict[str, tuple[float, str]] = {}

    for top_hits in pyhmmer.hmmer.hmmsearch(hmms, digital_seqs):
        marker_name = _as_text(getattr(top_hits.query, "name", ""))
        for hit in top_hits:
            if hit.evalue > evalue_cutoff:
                continue
            gene_id = _as_text(hit.name)
            prev = best.get(gene_id)
            if prev is None or hit.evalue < prev[0]:
                best[gene_id] = (hit.evalue, marker_name)

    return {gene_id: name for gene_id, (_, name) in sorted(best.items())}


def collect_prophage_windows(
    genomes: Sequence["GenomeData"],
    shells_by_genome: Sequence[Sequence["ShellHit"]],
    prophage_window: int,
) -> Dict[str, List[Gene]]:
    """Collect the genes within +/- `prophage_window` genes of each shell.

    This window is deliberately wider than the cargo window. A prophage spans
    far more than the handful of genes that could plausibly be cargo, so a +/-5
    gene cargo window rarely contains two distinct phage marker families even
    for an obvious prophage capsid. Andreas & Giessen (2021) used a 20-ORF
    genome neighborhood when screening out phage contamination, while using a
    5000 bp radius to assign cargo; this mirrors that split.

    Contigs are treated as linear, consistent with `extract_neighborhood`, so a
    window does not wrap around the origin of a circular replicon.

    Returns
    -------
    dict of str -> list of Gene
        Maps shell gene_id to the genes around it, excluding the shell itself.
    """
    windows: Dict[str, List[Gene]] = {}
    for genome, shells in zip(genomes, shells_by_genome, strict=True):
        for shell in shells:
            contig_genes = genome.genes_on_contig(shell.contig_id)
            index = next(
                (i for i, g in enumerate(contig_genes) if g.gene_id == shell.gene_id),
                None,
            )
            if index is None:
                windows[shell.gene_id] = []
                continue
            lo = max(0, index - prophage_window)
            hi = min(len(contig_genes), index + prophage_window + 1)
            windows[shell.gene_id] = [
                g for g in contig_genes[lo:hi] if g.gene_id != shell.gene_id
            ]
    return windows


def assess_prophage_context(
    candidates_by_shell: Dict[str, Sequence[Gene]],
    marker_profile_paths: Iterable[Path | str],
    evalue_cutoff: float = 1e-5,
) -> Dict[str, ProphageEvidence]:
    """Assess prophage context for each shell from its neighborhood genes.

    Profiles are loaded once and every neighborhood gene is scanned in a single
    hmmsearch pass, then results are attributed back to each shell. A gene shared
    by two overlapping windows therefore counts for both.

    Parameters
    ----------
    candidates_by_shell : dict of str -> sequence of Gene
        Maps shell gene_id to the genes in its window.
    marker_profile_paths : iterable of Path or str
        Phage marker HMM files. If empty, every shell gets empty evidence.
    evalue_cutoff : float
        Maximum full-sequence E-value for a marker call.

    Returns
    -------
    dict of str -> ProphageEvidence
        One entry per shell, including shells with no markers found.
    """
    paths = [Path(p) for p in marker_profile_paths]
    evidence = {
        shell_id: ProphageEvidence(shell_gene_id=shell_id)
        for shell_id in candidates_by_shell
    }
    if not paths:
        return evidence

    # Deduplicate genes across overlapping windows before the search.
    unique_genes: Dict[str, Gene] = {}
    for genes in candidates_by_shell.values():
        for g in genes:
            unique_genes.setdefault(g.gene_id, g)

    marker_by_gene = scan_genes_for_phage_markers(
        list(unique_genes.values()),
        paths,
        evalue_cutoff=evalue_cutoff,
    )
    if not marker_by_gene:
        logger.info("No phage markers found in any shell neighborhood")
        return evidence

    for shell_id, genes in candidates_by_shell.items():
        for g in genes:
            marker = marker_by_gene.get(g.gene_id)
            if marker:
                evidence[shell_id].marker_families[g.gene_id] = marker

    flagged = sorted(
        shell_id for shell_id, ev in evidence.items() if ev.is_probable_prophage
    )
    if flagged:
        logger.warning(
            "%d shell locus/loci show probable prophage context "
            "(>=%d distinct phage marker families in the window): %s. "
            "Encapsulin shells share the HK97 fold with phage capsids, so these "
            "may not be encapsulins. Use --drop-prophage-loci to exclude them.",
            len(flagged),
            PHAGE_MARKER_FLAG_THRESHOLD,
            ", ".join(flagged),
        )
    return evidence
