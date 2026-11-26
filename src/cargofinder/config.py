"""Configuration constants, biological parameters, and default weights for cargofinder.

Encapsulin background:
Encapsulins are prokaryotic protein nanocompartments sharing the HK97-like
fold (Pfam clan CL0373). They self-assemble into hollow shells (24-42 nm)
enclosing specific cargo enzymes.

Families and v0.1 scope:
- Family 1 (Linocin_M18 / PF04454): The canonical and most widespread encapsulin
  family. Well characterized in Thermotoga maritima, Mycobacterium tuberculosis,
  Myxococcus xanthus, and Rhodococcus jostii. Cargo enzymes include ferritin-like
  proteins (Flp) and dye-decolorizing peroxidases (DyP). Covered by default in v0.1.
- Family 2A (PF19307, "Type 2A encapsulin shell protein SrpI-like", IPR045641):
  fetchable with `fetch-profiles --extra PF19307`, not scanned by default.
- PF05065 ("Phage capsid family", IPR054612): a general phage capsid family, so it
  is expected to hit prophages too. Fetchable, not scanned by default.
- PF08967 ("Putative type 4B encapsulin shell protein", ENCP4): fetchable, not
  scanned by default.

Names and accessions were checked against the InterPro API on 2026-10-31
(https://www.ebi.ac.uk/interpro/api/entry/pfam/<accession>).
"""

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence


@dataclass(frozen=True)
class ShellProfileDef:
    """Definition of an encapsulin shell HMM profile."""

    family: str
    pfam_id: str
    description: str
    default_evalue: float = 1e-5


# Default profiles scanned by cargofinder v0.1
DEFAULT_SHELL_PROFILES: List[ShellProfileDef] = [
    ShellProfileDef(
        family="Family 1 (Linocin_M18)",
        pfam_id="PF04454",
        description=(
            "Canonical Family 1 encapsulin shell protein "
            "(Pfam PF04454, InterPro IPR007544)"
        ),
        default_evalue=1e-5,
    )
]

# Additional profiles that can be fetched on request
EXTENDED_SHELL_PROFILES: Dict[str, ShellProfileDef] = {
    "PF19307": ShellProfileDef(
        family="Family 2A",
        pfam_id="PF19307",
        description="Type 2A encapsulin shell protein SrpI-like (InterPro IPR045641)",
        default_evalue=1e-5,
    ),
    "PF05065": ShellProfileDef(
        family="Phage capsid family",
        pfam_id="PF05065",
        description="Phage capsid family (InterPro IPR054612); expect prophage hits",
        default_evalue=1e-5,
    ),
    "PF08967": ShellProfileDef(
        family="Family 4B",
        pfam_id="PF08967",
        description="Putative type 4B encapsulin shell protein (ENCP4)",
        default_evalue=1e-5,
    ),
}


def select_profiles(extra: Sequence[str] = ()) -> List[ShellProfileDef]:
    """Return the default profiles plus the requested extended profiles."""
    unknown = sorted(set(extra) - set(EXTENDED_SHELL_PROFILES))
    if unknown:
        raise ValueError(
            f"Unknown extra profile(s) {unknown}. "
            f"Available: {sorted(EXTENDED_SHELL_PROFILES)}"
        )
    return DEFAULT_SHELL_PROFILES + [EXTENDED_SHELL_PROFILES[e] for e in extra]


@dataclass
class ScoringWeights:
    """Additive scoring component weights.

    Score = w_distance * S_distance
          + w_operon * S_operon
          + w_tp_consensus * S_tp_consensus
          + w_tail_features * S_tail_features
          + w_recurrence * S_recurrence

    Rationale:
    - Adjacency/Synteny: Cargo enzymes are strongly enriched in the immediate genomic
      vicinity of the shell (often +/- 1 to 3 genes).
    - Operon likelihood: Genes on the same strand with short intergenic distance
      (e.g., < 150 bp) are frequently co-transcribed.
    - TP consensus similarity: BLOSUM62 similarity of the terminal slice to the
      published Family 1 targeting-peptide consensus SDGSLGIGSLKRS (Giessen &
      Silver 2017; used as a blastp query by Andreas & Giessen 2021). This is the
      one tail signal with a literature basis, so it outweighs raw composition.
    - Targeting peptide composition: Weakly conserved terminal extension
      (hydrophobic/amphipathic patch, Gly/Ser flexibility). This component is a
      cargofinder heuristic with no published basis, and the v0.1 benchmark showed
      it added almost nothing over adjacency; it is deliberately down-weighted.
    - Recurrence: The same gene family appearing next to encapsulins across multiple
      independent loci/genomes provides strong evolutionary support.
    """

    w_distance: float = 0.30
    w_operon: float = 0.20
    w_tp_consensus: float = 0.20
    w_tail_features: float = 0.15
    w_recurrence: float = 0.15

    def validate(self) -> None:
        """Validate that weights are non-negative and sum to 1.0."""
        weights = [
            self.w_distance,
            self.w_operon,
            self.w_tp_consensus,
            self.w_tail_features,
            self.w_recurrence,
        ]
        if any(w < 0 for w in weights):
            raise ValueError(f"Weights must be non-negative: {self}")
        total = sum(weights)
        if not (0.99 <= total <= 1.01):
            raise ValueError(f"Weights must sum to 1.0 (got {total:.4f}): {self}")


def load_weights(weights_arg: str) -> ScoringWeights:
    """Load scoring weights from a JSON file path or an inline JSON string.

    Keys left out keep their default value; unknown keys are rejected so that a
    typo cannot silently fall back to a default.
    """
    path = Path(weights_arg)
    try:
        if path.is_file():
            data = json.loads(path.read_text(encoding="utf-8"))
        else:
            data = json.loads(weights_arg)
    except json.JSONDecodeError as exc:
        raise ValueError(
            f"Failed to parse --weights as JSON file or JSON string: {exc}"
        ) from exc

    known = set(ScoringWeights.__dataclass_fields__)
    unknown = sorted(set(data) - known)
    if unknown:
        raise ValueError(
            f"Unknown scoring weight(s) {unknown}. Valid keys: {sorted(known)}"
        )
    weights = ScoringWeights(**{k: float(v) for k, v in data.items()})
    weights.validate()
    return weights


@dataclass
class RunConfig:
    """Global execution configuration for cargofinder."""

    window_size: int = 5
    tail_length: int = 20
    max_operon_gap_bp: int = 150
    evalue_cutoff: float = 1e-5
    min_shell_length_aa: int = 150
    weights: ScoringWeights = field(default_factory=ScoringWeights)
    profile_dir: Optional[Path] = None
    mmseqs_min_seq_id: float = 0.3
    mmseqs_coverage: float = 0.8
    termini: str = "c"

    def validate(self) -> None:
        """Reject parameter values that would produce a meaningless run."""
        if self.window_size < 1:
            raise ValueError(f"Window must be at least 1 gene (got {self.window_size})")
        if self.tail_length < 1:
            raise ValueError(
                f"Tail length must be at least 1 residue (got {self.tail_length})"
            )
        if self.max_operon_gap_bp < 0:
            raise ValueError(
                f"Max operon gap cannot be negative (got {self.max_operon_gap_bp})"
            )
        if self.evalue_cutoff <= 0:
            raise ValueError(f"E-value must be positive (got {self.evalue_cutoff})")
        if self.termini not in ("c", "both"):
            raise ValueError(f"Termini must be 'c' or 'both' (got {self.termini!r})")
        self.weights.validate()
