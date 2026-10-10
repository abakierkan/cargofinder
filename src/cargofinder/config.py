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
- Family 2 (PF19307): Associated with specialized biosynthetic and metabolic cargo.
  Configurable by user in v0.1, but not included in default scan profile.
- Family 3 (PF05065 / Phage capsid-like): Found in biosynthetic gene clusters;
  shares phage capsid homology. Out-of-scope for v0.1 default scan.
- Family 4 (PF08967 / DUF1884): Minimalist A-domain encapsulins in thermophiles.
  Out-of-scope for v0.1 default scan.
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional


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

# Additional families documented for future expansion or user configuration
EXTENDED_SHELL_PROFILES: Dict[str, ShellProfileDef] = {
    "PF19307": ShellProfileDef(
        family="Family 2",
        pfam_id="PF19307",
        description="Family 2 encapsulin shell protein (InterPro IPR041443)",
        default_evalue=1e-5,
    ),
    "PF05065": ShellProfileDef(
        family="Family 3",
        pfam_id="PF05065",
        description="Phage capsid family / Family 3 encapsulin-associated",
        default_evalue=1e-5,
    ),
    "PF08967": ShellProfileDef(
        family="Family 4",
        pfam_id="PF08967",
        description="DUF1884 / Family 4 minimalist encapsulin",
        default_evalue=1e-5,
    ),
}


@dataclass
class ScoringWeights:
    """Additive scoring component weights.

    Score = w_distance * S_distance
          + w_operon * S_operon
          + w_tail_features * S_tail_features
          + w_recurrence * S_recurrence

    Rationale:
    - Adjacency/Synteny: Cargo enzymes are strongly enriched in the immediate genomic
      vicinity of the shell (often +/- 1 to 3 genes).
    - Operon likelihood: Genes on the same strand with short intergenic distance
      (e.g., < 150 bp) are frequently co-transcribed.
    - Targeting peptide: Weakly conserved terminal extension (typically C-terminal,
      hydrophobic/amphipathic patch, Gly/Ser flexibility).
    - Recurrence: The same gene family appearing next to encapsulins across multiple
      independent loci/genomes provides strong evolutionary support.
    """

    w_distance: float = 0.35
    w_operon: float = 0.25
    w_tail_features: float = 0.25
    w_recurrence: float = 0.15

    def validate(self) -> None:
        """Validate that weights are non-negative and sum to 1.0."""
        weights = [
            self.w_distance,
            self.w_operon,
            self.w_tail_features,
            self.w_recurrence,
        ]
        if any(w < 0 for w in weights):
            raise ValueError(f"Weights must be non-negative: {self}")
        total = sum(weights)
        if not (0.99 <= total <= 1.01):
            raise ValueError(f"Weights must sum to 1.0 (got {total:.4f}): {self}")


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
