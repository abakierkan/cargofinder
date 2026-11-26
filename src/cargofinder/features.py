"""Terminal peptide feature extraction and physicochemical profiling.

Encapsulin cargo recruitment:
In canonical encapsulin systems (Family 1), cargo enzymes possess short terminal
targeting peptides (typically 10-35 residues at the C-terminus, occasionally N-terminus)
that bind hydrophobic pockets on the interior surface of the shell.

Key physicochemical signatures:
1. Hydrophobic anchor: A short hydrophobic patch that docks into the shell lumen pocket
   (quantified using the Kyte-Doolittle scale, Kyte & Doolittle 1982 JMB).
2. Gly/Ser-rich flexible linker: Small flexible residues connecting the
   globular cargo core to the targeting peptide.
3. Similarity to the published Family 1 targeting-peptide consensus
   SDGSLGIGSLKRS (Giessen & Silver 2017 Nat Microbiol 2:17029; used as a blastp
   query by Andreas & Giessen 2021 Nat Commun 12:4748). See substitution.py.
4. User-configurable sequence motif regexes.

Sequence complexity (Shannon entropy, low-complexity score) is reported for
inspection but is not used in the score: we found no evidence for either direction.
"""

from __future__ import annotations

import json
import math
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Optional, Tuple

from cargofinder.substitution import (
    FAMILY1_TP_CONSENSUS,
    best_ungapped_similarity,
)

# Kyte & Doolittle (1982) J. Mol. Biol. 157:105-132
KYTE_DOOLITTLE: Dict[str, float] = {
    "I": 4.5,
    "V": 4.2,
    "L": 3.8,
    "F": 2.8,
    "C": 2.5,
    "M": 1.9,
    "A": 1.8,
    "G": -0.4,
    "T": -0.7,
    "S": -0.8,
    "W": -0.9,
    "Y": -1.3,
    "P": -1.6,
    "H": -3.2,
    "E": -3.5,
    "Q": -3.5,
    "D": -3.5,
    "N": -3.5,
    "K": -3.9,
    "R": -4.5,
}

HYDROPHOBIC_RESIDUES = frozenset("AVILMFWP")
POSITIVE_RESIDUES = frozenset("KR")
NEGATIVE_RESIDUES = frozenset("DE")


@dataclass
class PeptideFeatures:
    """Computed physicochemical and motif features for a protein's termini."""

    c_term_seq: str
    n_term_seq: str
    tail_length: int
    c_kd_hydrophobicity: float
    c_hydrophobic_fraction: float
    c_gly_ser_fraction: float
    c_shannon_entropy: float
    c_low_complexity_score: float
    c_net_charge: float
    c_tp_consensus_similarity: float
    n_kd_hydrophobicity: float
    n_gly_ser_fraction: float
    n_hydrophobic_fraction: float
    n_tp_consensus_similarity: float
    motif_matches: Dict[str, bool] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, float | str | int | bool]:
        """Flatten features into a dictionary for reporting."""
        d: Dict[str, float | str | int | bool] = {
            "c_term_seq": self.c_term_seq,
            "n_term_seq": self.n_term_seq,
            "tail_length": self.tail_length,
            "c_kd_hydrophobicity": round(self.c_kd_hydrophobicity, 4),
            "c_hydrophobic_fraction": round(self.c_hydrophobic_fraction, 4),
            "c_gly_ser_fraction": round(self.c_gly_ser_fraction, 4),
            "c_shannon_entropy": round(self.c_shannon_entropy, 4),
            "c_low_complexity_score": round(self.c_low_complexity_score, 4),
            "c_net_charge": round(self.c_net_charge, 2),
            "c_tp_consensus_similarity": round(self.c_tp_consensus_similarity, 4),
            "n_kd_hydrophobicity": round(self.n_kd_hydrophobicity, 4),
            "n_gly_ser_fraction": round(self.n_gly_ser_fraction, 4),
            "n_hydrophobic_fraction": round(self.n_hydrophobic_fraction, 4),
            "n_tp_consensus_similarity": round(self.n_tp_consensus_similarity, 4),
        }
        for name, matched in self.motif_matches.items():
            d[f"motif_{name}"] = matched
        return d


def extract_terminal_peptides(seq: str, tail_length: int = 20) -> Tuple[str, str]:
    """Extract N-terminal and C-terminal peptide sequences.

    Parameters
    ----------
    seq : str
        Full protein amino acid sequence.
    tail_length : int
        Number of residues to slice from each terminus (default 20).

    Returns
    -------
    tuple of (str, str)
        (n_terminal_slice, c_terminal_slice)
    """
    if tail_length < 1:
        raise ValueError(f"Tail length must be at least 1 residue (got {tail_length})")
    clean = seq.strip().upper().rstrip("*")
    if not clean:
        return "", ""
    if len(clean) <= tail_length:
        return clean, clean
    return clean[:tail_length], clean[-tail_length:]


def compute_hydrophobicity_kd(peptide: str) -> float:
    """Compute mean Kyte-Doolittle hydropathy for a peptide string."""
    if not peptide:
        return 0.0
    scores = [KYTE_DOOLITTLE.get(aa, 0.0) for aa in peptide]
    return sum(scores) / len(scores)


def compute_hydrophobic_fraction(peptide: str) -> float:
    """Compute fraction of aliphatic/aromatic hydrophobic residues (AVILMFWP)."""
    if not peptide:
        return 0.0
    count = sum(1 for aa in peptide if aa in HYDROPHOBIC_RESIDUES)
    return count / len(peptide)


def compute_gly_ser_fraction(peptide: str) -> float:
    """Compute fraction of Glycine and Serine (hallmarks of flexible linkers)."""
    if not peptide:
        return 0.0
    count = sum(1 for aa in peptide if aa in ("G", "S"))
    return count / len(peptide)


def compute_shannon_entropy(peptide: str) -> float:
    """Compute Shannon entropy H = -sum(p_i * log2(p_i)) over amino acid frequencies."""
    if not peptide:
        return 0.0
    length = len(peptide)
    counts: Dict[str, int] = {}
    for aa in peptide:
        counts[aa] = counts.get(aa, 0) + 1

    entropy = 0.0
    for c in counts.values():
        p = c / length
        if p > 0.0:
            entropy -= p * math.log2(p)
    return entropy


def compute_low_complexity_score(peptide: str) -> float:
    """Compute a normalized low-complexity score between 0.0 and 1.0.

    1.0 indicates maximum repetitiveness / low complexity (e.g. homopolymer).
    0.0 indicates maximal diverse distribution across 20 amino acids.
    """
    if not peptide or len(peptide) <= 1:
        return 1.0
    h = compute_shannon_entropy(peptide)
    # Maximum possible entropy for this peptide length across 20 possible amino acids
    max_h = math.log2(min(len(peptide), 20))
    if max_h <= 0.0:
        return 1.0
    norm_entropy = min(1.0, h / max_h)
    return 1.0 - norm_entropy


def compute_net_charge(peptide: str) -> float:
    """Approximate net charge at neutral pH (+1 for K/R, -1 for D/E)."""
    if not peptide:
        return 0.0
    pos = sum(1 for aa in peptide if aa in POSITIVE_RESIDUES)
    neg = sum(1 for aa in peptide if aa in NEGATIVE_RESIDUES)
    return float(pos - neg)


def compute_tp_consensus_similarity(
    peptide: str,
    reference: str = FAMILY1_TP_CONSENSUS,
) -> float:
    """Normalised BLOSUM62 similarity of a terminal slice to a reference TP.

    Defaults to the published Family 1 targeting-peptide consensus
    SDGSLGIGSLKRS (Giessen & Silver 2017; used as a blastp query by
    Andreas & Giessen 2021). See `substitution.best_ungapped_similarity`.

    Caveat worth knowing before you read this number: the consensus is itself
    Gly/Ser-rich, so generic Gly/Ser linkers score around 0.3 as background. A
    plain GGSGGSGGSGG control peptide, which Altenburg et al. (2021) showed does
    not bind the shell pocket at all, scores about as high as some native TPs.
    This feature does not separate a real anchor from a flexible linker on its
    own; the anchor motifs in data/published_tp_motifs.json are what address
    that, and neither is proof of encapsulation.

    Returns
    -------
    float
        Similarity in [0.0, 1.0]; 1.0 is an exact full-length match.
    """
    return best_ungapped_similarity(peptide, reference=reference)


def load_motif_patterns(motif_file: Path | str) -> Dict[str, str]:
    """Load regex motif patterns from a JSON file.

    Format expected:
    {
      "motif_name": {
        "regex": "...",
        "citation": "...",
        "description": "..."
      }
    }
    or simple mapping:
    {
      "motif_name": "REGEX_PATTERN"
    }

    Keys beginning with an underscore are treated as documentation and skipped,
    so a motif file can carry notes without creating an output column.
    """
    p = Path(motif_file)
    if not p.is_file():
        raise FileNotFoundError(f"Motif file not found: {motif_file}")

    data = json.loads(p.read_text(encoding="utf-8"))
    patterns: Dict[str, str] = {}
    for k, v in data.items():
        # Keys starting with an underscore are documentation, not motifs, so a
        # motif file can carry provenance notes without adding output columns.
        if k.startswith("_"):
            continue
        if isinstance(v, str):
            patterns[k] = v
        elif isinstance(v, dict) and "regex" in v:
            patterns[k] = v["regex"]
        else:
            raise ValueError(f"Malformed motif definition for '{k}' in {motif_file}")

    # Validate all patterns compile
    for name, pat in patterns.items():
        try:
            re.compile(pat)
        except re.error as exc:
            raise ValueError(
                f"Invalid regex for motif '{name}' ({pat!r}): {exc}"
            ) from exc

    return patterns


def match_motifs(peptide: str, patterns: Dict[str, str]) -> Dict[str, bool]:
    """Test peptide against compiled regex patterns."""
    matches: Dict[str, bool] = {}
    for name, pat in patterns.items():
        matches[name] = bool(re.search(pat, peptide))
    return matches


def compute_peptide_features(
    sequence: str,
    tail_length: int = 20,
    motifs: Optional[Dict[str, str]] = None,
) -> PeptideFeatures:
    """Compute complete suite of terminal peptide features for a candidate protein.

    Parameters
    ----------
    sequence : str
        Protein sequence.
    tail_length : int
        Tail slice length in residues (default 20).
    motifs : dict of str -> str, optional
        User-supplied regex patterns to scan against the C-terminal tail.

    Returns
    -------
    PeptideFeatures
        Full physicochemical and motif feature set.
    """
    n_tail, c_tail = extract_terminal_peptides(sequence, tail_length=tail_length)

    c_kd = compute_hydrophobicity_kd(c_tail)
    c_hydro = compute_hydrophobic_fraction(c_tail)
    c_gly_ser = compute_gly_ser_fraction(c_tail)
    c_entropy = compute_shannon_entropy(c_tail)
    c_low_comp = compute_low_complexity_score(c_tail)
    c_charge = compute_net_charge(c_tail)
    c_tp_sim = compute_tp_consensus_similarity(c_tail)

    n_kd = compute_hydrophobicity_kd(n_tail)
    n_gly_ser = compute_gly_ser_fraction(n_tail)
    n_hydro = compute_hydrophobic_fraction(n_tail)
    n_tp_sim = compute_tp_consensus_similarity(n_tail)

    motif_matches = match_motifs(c_tail, motifs) if motifs else {}

    return PeptideFeatures(
        c_term_seq=c_tail,
        n_term_seq=n_tail,
        tail_length=tail_length,
        c_kd_hydrophobicity=c_kd,
        c_hydrophobic_fraction=c_hydro,
        c_gly_ser_fraction=c_gly_ser,
        c_shannon_entropy=c_entropy,
        c_low_complexity_score=c_low_comp,
        c_net_charge=c_charge,
        c_tp_consensus_similarity=c_tp_sim,
        n_kd_hydrophobicity=n_kd,
        n_gly_ser_fraction=n_gly_ser,
        n_hydrophobic_fraction=n_hydro,
        n_tp_consensus_similarity=n_tp_sim,
        motif_matches=motif_matches,
    )
