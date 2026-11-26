"""BLOSUM62 and an ungapped similarity scan against a reference peptide.

Why this module exists
----------------------
`cargofinder` v0.1 scored terminal peptides only by amino-acid composition
(hydrophobic fraction and Gly/Ser fraction). That was an invented heuristic, and
the v0.1 benchmark showed it added almost nothing over adjacency alone.

Published work instead compares candidate termini against a *reference targeting
peptide*. Andreas & Giessen (2021) identified Family 1 cargo proteins by taking
every protein encoded within 5000 bp of a Family 1 shell gene and running
`blastp` against the Family 1 targeting-peptide consensus `SDGSLGIGSLKRS`, with
parameters adjusted for short queries (Nat Commun 12:4748,
doi:10.1038/s41467-021-25071-y, Methods, "Genome neighborhood analysis"). That
consensus is from Giessen & Silver (2017) Nat Microbiol 2:17029,
doi:10.1038/nmicrobiol.2017.29.

This module reproduces the comparison in a small, deterministic, dependency-free
form: an ungapped BLOSUM62 window scan, normalised against the reference's
self-score. It is not BLAST. It reports no E-value and makes no claim of
statistical significance; it is a bounded similarity feature in [0, 1] whose
value can be read directly off the output table.

BLOSUM62 is from Henikoff S & Henikoff JG (1992) PNAS 89:10915-10919,
doi:10.1073/pnas.89.22.10915, as distributed by NCBI.
"""

from __future__ import annotations

from typing import Dict, Tuple

# Standard NCBI BLOSUM62, in the conventional residue order.
_BLOSUM62_ORDER = "ARNDCQEGHILKMFPSTWYV"

_BLOSUM62_ROWS = (
    "  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0",
    " -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3",
    " -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3",
    " -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3",
    "  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1",
    " -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2",
    " -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2",
    "  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3",
    " -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3",
    " -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3",
    " -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1",
    " -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2",
    " -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1",
    " -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1",
    " -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2",
    "  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2",
    "  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0",
    " -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3",
    " -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1",
    "  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4",
)


def _build_blosum62() -> Dict[Tuple[str, str], int]:
    matrix: Dict[Tuple[str, str], int] = {}
    for row_aa, row in zip(_BLOSUM62_ORDER, _BLOSUM62_ROWS, strict=True):
        values = [int(v) for v in row.split()]
        if len(values) != len(_BLOSUM62_ORDER):
            raise AssertionError(
                f"BLOSUM62 row for {row_aa!r} has {len(values)} entries, "
                f"expected {len(_BLOSUM62_ORDER)}"
            )
        for col_aa, value in zip(_BLOSUM62_ORDER, values, strict=True):
            matrix[(row_aa, col_aa)] = value
    return matrix


BLOSUM62: Dict[Tuple[str, str], int] = _build_blosum62()

# The Family 1 targeting-peptide consensus used as a BLAST query by
# Andreas & Giessen (2021); sequence from Giessen & Silver (2017).
FAMILY1_TP_CONSENSUS = "SDGSLGIGSLKRS"


def score_pair(a: str, b: str) -> int:
    """BLOSUM62 score for one residue pair; unknown residues score 0."""
    return BLOSUM62.get((a, b), 0)


def self_score(peptide: str) -> int:
    """Sum of BLOSUM62 diagonal scores for a peptide (its maximum achievable score)."""
    return sum(score_pair(aa, aa) for aa in peptide)


def best_ungapped_similarity(
    peptide: str,
    reference: str = FAMILY1_TP_CONSENSUS,
) -> float:
    """Normalised best ungapped BLOSUM62 similarity of `peptide` to `reference`.

    The reference is slid along the peptide one position at a time, including
    partial overlaps at both ends, and every alignment is scored with BLOSUM62.
    The best raw score is divided by the reference's self-score, so an exact
    full-length match of the reference returns 1.0.

    Negative-scoring alignments are clamped to 0.0, so the result is in [0, 1].
    This is a similarity feature, not a statistical test: there is no E-value and
    no significance claim. Short peptides can only ever reach a fraction of 1.0,
    because they cannot overlap the whole reference.

    Parameters
    ----------
    peptide : str
        Peptide to score, typically a terminal slice of a candidate protein.
    reference : str
        Reference peptide. Defaults to the published Family 1 TP consensus.

    Returns
    -------
    float
        Normalised similarity in [0.0, 1.0].
    """
    pep = peptide.strip().upper()
    ref = reference.strip().upper()
    if not pep or not ref:
        return 0.0

    denominator = self_score(ref)
    if denominator <= 0:
        return 0.0

    best = 0
    # Offsets allow the reference to hang off either end of the peptide, so that
    # a consensus-like stretch at the very terminus is still detected.
    for offset in range(-len(ref) + 1, len(pep)):
        total = 0
        for k, ref_aa in enumerate(ref):
            pos = offset + k
            if 0 <= pos < len(pep):
                total += score_pair(pep[pos], ref_aa)
        if total > best:
            best = total

    return min(1.0, best / denominator)
