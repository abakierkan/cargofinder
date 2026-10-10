"""Encapsulin shell detection using pyhmmer (HMMER3 Python bindings).

Encapsulin shells are self-assembling icosahedral or quasi-spherical protein
nanocompartments based on the viral HK97-fold. In Family 1 (Linocin_M18 / PF04454),
the monomer typically ranges between ~240-320 amino acids.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List

import pyhmmer
from pyhmmer import easel, plan7

from cargofinder.io import Gene, GenomeData

logger = logging.getLogger(__name__)


@dataclass
class ShellHit:
    """Represents a confident encapsulin shell gene identified by HMMER."""

    gene_id: str
    profile_name: str
    profile_acc: str
    evalue: float
    bitscore: float
    domain_evalue: float
    domain_count: int
    gene: Gene

    @property
    def contig_id(self) -> str:
        return self.gene.contig_id

    @property
    def start(self) -> int:
        return self.gene.start

    @property
    def end(self) -> int:
        return self.gene.end

    @property
    def strand(self) -> int:
        return self.gene.strand


def load_hmm_profiles(profile_paths: Iterable[Path | str]) -> List[plan7.HMM]:
    """Load HMMER3 profile models from one or more file paths.

    Parameters
    ----------
    profile_paths : iterable of Path or str
        Paths to .hmm or .hmm.gz files.

    Returns
    -------
    list of plan7.HMM
        Loaded HMM objects.
    """
    loaded: List[plan7.HMM] = []
    for path_str in profile_paths:
        p = Path(path_str)
        if not p.is_file():
            raise FileNotFoundError(f"HMM profile file not found: {p}")
        with plan7.HMMFile(str(p)) as hmm_file:
            for hmm in hmm_file:
                loaded.append(hmm)

    if not loaded:
        raise ValueError("No HMM profile models could be loaded from specified paths.")

    logger.debug("Loaded %d HMM profile model(s)", len(loaded))
    return loaded


def detect_shells(
    genome: GenomeData,
    profile_paths: Iterable[Path | str],
    evalue_cutoff: float = 1e-5,
    min_length_aa: int = 150,
) -> List[ShellHit]:
    """Scan genome protein sequences against encapsulin shell HMM profiles.

    Parameters
    ----------
    genome : GenomeData
        Annotated genome data with protein sequences.
    profile_paths : iterable of Path or str
        Paths to shell profile HMM files.
    evalue_cutoff : float
        Maximum full-sequence E-value threshold for considering a hit significant.
    min_length_aa : int
        Minimum protein length in amino acids. Encapsulin shell monomers
        contain an HK97 core fold that is rarely under 150 aa.

    Returns
    -------
    list of ShellHit
        Detected shell candidates, deduplicated by best E-value if a gene matches
        multiple profiles.
    """
    if len(genome.genes) == 0:
        logger.warning(
            "No genes present in genome %s to scan for shells", genome.genome_id
        )
        return []

    hmms = load_hmm_profiles(profile_paths)
    alphabet = easel.Alphabet.amino()

    # Convert genome proteins to Easel digital sequences
    digital_seqs: List[easel.DigitalSequence] = []
    for g in genome.genes:
        # Skip severely truncated ORFs below length threshold
        if g.length_aa < min_length_aa:
            continue
        text_seq = easel.TextSequence(
            name=g.gene_id.encode("utf-8"),
            sequence=g.protein_seq,
        )
        digital_seqs.append(text_seq.digitize(alphabet))

    if not digital_seqs:
        logger.warning(
            "No genes met the minimum length of %d aa in genome %s",
            min_length_aa,
            genome.genome_id,
        )
        return []

    # Map gene_id -> best ShellHit to avoid duplicate reports across profiles
    best_hits_by_gene: Dict[str, ShellHit] = {}

    # Run hmmsearch via pyhmmer
    for top_hits in pyhmmer.hmmer.hmmsearch(hmms, digital_seqs):
        hmm_name = (getattr(top_hits.query, "name", "") if top_hits.query else "") or ""
        hmm_acc = (
            getattr(top_hits.query, "accession", "") if top_hits.query else ""
        ) or ""

        for hit in top_hits:
            if hit.evalue > evalue_cutoff:
                continue

            hit_name = (
                hit.name.decode("utf-8")
                if isinstance(hit.name, bytes)
                else str(hit.name)
            )
            gene_obj = genome.get_gene(hit_name)
            if gene_obj is None:
                continue

            best_dom_evalue = (
                min(dom.i_evalue for dom in hit.domains) if hit.domains else hit.evalue
            )

            cand = ShellHit(
                gene_id=hit_name,
                profile_name=hmm_name,
                profile_acc=hmm_acc,
                evalue=hit.evalue,
                bitscore=hit.score,
                domain_evalue=best_dom_evalue,
                domain_count=len(hit.domains),
                gene=gene_obj,
            )

            existing = best_hits_by_gene.get(hit_name)
            if existing is None or cand.evalue < existing.evalue:
                best_hits_by_gene[hit_name] = cand

    hits = sorted(
        best_hits_by_gene.values(),
        key=lambda h: (h.gene.contig_id, h.gene.start),
    )

    logger.info(
        "Detected %d encapsulin shell candidate(s) in genome %s",
        len(hits),
        genome.genome_id,
    )
    return hits
