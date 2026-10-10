"""Gene calling using Pyrodigal (Python bindings for Prodigal)."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict, List, Optional

import pyrodigal

from cargofinder.io import Gene, GenomeData, read_fasta_sequences

logger = logging.getLogger(__name__)

# Minimum sequence length in nucleotides required for Prodigal training mode
PRODIGAL_MIN_TRAIN_NT = 20000


def call_genes_pyrodigal(
    contigs: Dict[str, str],
    genome_id: str = "genome",
    meta: Optional[bool] = None,
) -> GenomeData:
    """Call genes from contig DNA sequences using Pyrodigal.

    Parameters
    ----------
    contigs : dict of str -> str
        Contig IDs mapped to uppercase nucleotide sequences.
    genome_id : str
        Identifier for the genome.
    meta : bool, optional
        Whether to force metagenomic / short-contig mode. If None, auto-detected:
        uses trained mode if total sequence length >= 20,000 nt, else meta mode.

    Returns
    -------
    GenomeData
        Genome container with called genes.
    """
    total_nt = sum(len(seq) for seq in contigs.values())

    if meta is None:
        use_meta = total_nt < PRODIGAL_MIN_TRAIN_NT
    else:
        use_meta = meta

    if use_meta:
        logger.info(
            "Using Pyrodigal meta mode (total genome length %d nt < %d nt)",
            total_nt,
            PRODIGAL_MIN_TRAIN_NT,
        )
        finder = pyrodigal.GeneFinder(meta=True)
    else:
        logger.info(
            "Training Pyrodigal on %d nt of sequence across %d contig(s)",
            total_nt,
            len(contigs),
        )
        finder = pyrodigal.GeneFinder(meta=False)
        # Train across all concatenated or primary contigs
        # We join contigs with 50 'N's to prevent chimeric ORF calling across contigs
        joined_seq = "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN".join(
            contigs.values()
        )
        finder.train(joined_seq.encode("ascii"))

    all_genes: List[Gene] = []

    for contig_id, seq in contigs.items():
        encoded = seq.encode("ascii")
        preds = finder.find_genes(encoded)

        for idx, pred in enumerate(preds):
            gene_id = f"{contig_id}_{idx + 1:04d}"
            # pred.begin and pred.end are 1-based inclusive
            # pred.strand is +1 or -1
            aa_seq = pred.translate()
            if aa_seq.endswith("*"):
                aa_seq = aa_seq[:-1]

            all_genes.append(
                Gene(
                    gene_id=gene_id,
                    contig_id=contig_id,
                    start=pred.begin,
                    end=pred.end,
                    strand=pred.strand,
                    gene_index=idx,
                    protein_seq=aa_seq,
                    product=f"Pyrodigal predicted protein {gene_id}",
                )
            )

    logger.info("Predicted %d genes across %d contigs", len(all_genes), len(contigs))
    return GenomeData(genome_id=genome_id, genes=all_genes, contigs=contigs)


def call_genes_from_genome_fasta(
    fasta_path: Path | str,
    genome_id: Optional[str] = None,
    meta: Optional[bool] = None,
) -> GenomeData:
    """Convenience wrapper to load a FASTA file and run Pyrodigal gene calling."""
    p = Path(fasta_path)
    if genome_id is None:
        genome_id = p.stem
    contigs = read_fasta_sequences(p)
    return call_genes_pyrodigal(contigs, genome_id=genome_id, meta=meta)
