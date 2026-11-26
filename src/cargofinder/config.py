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
- Family 2 (PF19307, "Type 2A encapsulin shell protein SrpI-like", IPR045641):
  scanned by default since v0.2. Family 2 is the largest group in the
  Andreas & Giessen (2021) census (3,523 systems against Family 1's 2,383), so
  scanning Family 1 alone missed most encapsulin systems. Despite being named
  for Family 2A, this profile also detects Family 2B shells in practice: in
  M. xanthus it finds MXAN_6248 and MXAN_6249 at E ~ 1e-136, which RefSeq
  annotates as "family 2B encapsulin nanocompartment shell protein" and which
  are the adjacent shell-gene pair that Andreas & Giessen describe as the
  hallmark of 2B operons. Family 2 cargo is recruited by an N-terminal domain,
  so TARGETING_MODE_BY_PROFILE routes these loci to N-terminal scoring.
- PF05065 ("Phage capsid family", IPR054612): a general phage capsid family, so it
  is expected to hit prophages too. Fetchable, not scanned by default.
- PF08967 ("Putative type 4B encapsulin shell protein", ENCP4): fetchable, not
  scanned by default.

Names and accessions were checked against the InterPro API on 2026-10-31
(https://www.ebi.ac.uk/interpro/api/entry/pfam/<accession>).
"""

import json
import re
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
    ),
    ShellProfileDef(
        family="Family 2 (SrpI-like)",
        pfam_id="PF19307",
        description=(
            "Type 2A encapsulin shell protein SrpI-like (Pfam PF19307, "
            "InterPro IPR045641); also detects Family 2B shells in practice"
        ),
        default_evalue=1e-5,
    ),
]

# Additional profiles that can be fetched on request
EXTENDED_SHELL_PROFILES: Dict[str, ShellProfileDef] = {
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


# How each shell family recruits its cargo, which decides what "targeting signal"
# means for candidates at that locus.
#
# Andreas & Giessen (2021) Nat Commun 12:4748, doi:10.1038/s41467-021-25071-y:
# "In classical encapsulins (Family 1), encapsulation is mediated by short
# C-terminal peptide sequences referred to as targeting peptides (TPs) or
# cargo-loading peptides (CLPs), while for Family 2 systems, larger N-terminal
# protein domains are proposed to mediate encapsulation."
#
# Scoring every candidate against a C-terminal Family 1 consensus is therefore
# wrong at a Family 2 locus, where no short targeting peptide exists at all:
# "While no specific targeting peptide has been identified in CD systems, the
# unannotated N-terminal domain has been shown to be responsible for mediating
# encapsulation."
TARGETING_MODE_C_TERMINAL = "c_terminal_tp"
TARGETING_MODE_N_TERMINAL = "n_terminal_extension"

#: Maps shell Pfam accession to the targeting mode used for its candidates.
TARGETING_MODE_BY_PROFILE: Dict[str, str] = {
    "PF04454": TARGETING_MODE_C_TERMINAL,  # Family 1
    "PF19307": TARGETING_MODE_N_TERMINAL,  # Family 2A
}

#: Mode used when the shell profile is not in the map above. Family 4B (PF08967)
#: and any user-supplied profile land here. This is a documented fallback, not a
#: claim about those families; every output row records the mode actually used.
TARGETING_MODE_DEFAULT = TARGETING_MODE_C_TERMINAL

# N-terminal extension length, in residues, bounding the Family 2 targeting
# score. Anchored to the truncation series in Nichols RJ et al. (2021) eLife
# 10:e59288, doi:10.7554/eLife.59288, for the S. elongatus SrpI cysteine
# desulfurase (CyD): "The entire 225 amino acid N-terminal domain from CyD fused
# to sfGFP (CyD 1-225-sfGFP) yielded the highest loading signal. Targeting with
# the first 100 amino acids of the NTD also functioned, albeit not as
# efficiently"; while "CyD 155-183-sfGFP, containing the conserved motif of
# 'AASPYYFLDG', was not sufficient to sequester cargo", and a mutant lacking the
# whole NTD "was not measurably encapsulated".
#
# So an extension of about 30 residues is known not to work and 225 is known to
# work best. The score ramps linearly between them; 100 residues lands at ~0.36,
# consistent with "functioned, albeit not as efficiently". The linear shape
# between the anchors is a judgment call.
NTERM_EXTENSION_MIN_AA = 30
NTERM_EXTENSION_FULL_AA = 225


@dataclass(frozen=True)
class CargoClassDef:
    """A Pfam family that marks a known encapsulin cargo class."""

    pfam_id: str
    cargo_class: str
    description: str
    #: True only for families defined specifically as encapsulin cargo. False for
    #: broad enzyme families that merely contain cargo among much else.
    encapsulin_specific: bool = False


# Cargo class profiles.
#
# Why this exists: Andreas & Giessen (2021) Nat Commun 12:4748
# (doi:10.1038/s41467-021-25071-y) assigned cargo not by position alone but by
# identity, building "HMM templates representative of each cargo class" and
# classifying the hits "based on their Pfam or Interpro annotation". The cargo
# classes below are the ones that paper names.
#
# Family 1 ("at least 6 operon types defined by the presence of 6 distinct and
# conserved cargo proteins"): DyP peroxidases are "the most abundant cargo class
# in Family 1" (1,505 systems), and ferritin-like proteins (FLPs) "comprise the
# second largest set", clustering into families that "largely correspond to the
# following Pfam families: Flp, rubrerythrin, CoQ7, Mn catalase, IMEF,
# hemerythrin, and bacterioferritin".
#
# Family 2: "The most commonly enriched genes associated with Family 2
# encapsulins encode for cysteine desulfurases (CD), terpene cyclases (TC),
# polyprenyl transferases (PT), and xylulose kinases (XK)". The paper gives the
# cysteine desulfurase domain accession directly as PF00266.
#
# Deliberately NOT included: serine O-acetyltransferases and rhodaneses. They are
# the most enriched accessory components of Family 2A operons, but that paper
# states plainly that such components "are not cargo proteins".
#
# Accessions and names were checked against the InterPro API on 2026-11-26
# (https://www.ebi.ac.uk/interpro/api/entry/pfam/<accession>).
CARGO_CLASS_PROFILES: List[CargoClassDef] = [
    # -- Encapsulin-cargo-specific families. Only these two are defined by Pfam
    # -- as encapsulin cargo; a hit here is far more meaningful than the rest.
    CargoClassDef(
        "PF22277",
        "encapsulated_ferritin",
        "EncFtn-like: encapsulated ferritin-like ferroxidase (R. rubrum EncFtn)",
        encapsulin_specific=True,
    ),
    CargoClassDef(
        "PF24309",
        "imef_flp",
        "IMEF_Flp: IMEF cargo Flp protein of a type 1 encapsulin nanocompartment",
        encapsulin_specific=True,
    ),
    # -- Dye-decolorizing peroxidase, the most abundant Family 1 cargo class.
    CargoClassDef("PF04261", "dyp_peroxidase", "Dyp_perox_N: DyP-type peroxidase"),
    CargoClassDef("PF20628", "dyp_peroxidase", "Dyp_perox_C: DyP-type peroxidase"),
    CargoClassDef("PF21105", "dyp_peroxidase", "DyP_N: DyP dimeric a+b barrel domain"),
    # -- Ferritin-like protein superfamily, the second largest Family 1 set.
    CargoClassDef("PF00210", "ferritin_like", "Ferritin: ferritin-like domain"),
    CargoClassDef(
        "PF13668",
        "ferritin_like",
        "Ferritin_2: ferritins, DPS family and bacterioferritins",
    ),
    CargoClassDef("PF02915", "rubrerythrin", "Rubrerythrin"),
    CargoClassDef(
        "PF05067", "mn_catalase", "Mn_catalase: manganese containing catalase"
    ),
    CargoClassDef("PF03232", "coq7", "COQ7: ubiquinone biosynthesis protein COQ7"),
    CargoClassDef(
        "PF01814", "hemerythrin", "Hemerythrin: hemerythrin HHE cation binding domain"
    ),
    # -- Family 2 cargo classes.
    CargoClassDef(
        "PF00266",
        "cysteine_desulfurase",
        "Aminotran_5: SufS-like cysteine desulfurase domain, named as PF00266 "
        "for Family 2 CD cargo by Andreas & Giessen (2021)",
    ),
    CargoClassDef(
        "PF03936", "terpene_cyclase", "Terpene_synth_C: terpene synthase metal binding"
    ),
    CargoClassDef(
        "PF00348", "polyprenyl_transferase", "polyprenyl_synt: polyprenyl synthetase"
    ),
    CargoClassDef(
        "PF00370", "xylulose_kinase", "FGGY_N: FGGY carbohydrate kinase, N-terminal"
    ),
]

# Subdirectory of the profile directory that holds cargo class profiles, so that
# find_profile_files() does not pick them up as shell profiles.
CARGO_CLASS_SUBDIR = "cargo_classes"

# Score given to a candidate matching a cargo class profile. These are judgment
# calls, not fitted values. An encapsulin-specific family scores full marks; a
# broad enzyme family scores less because most of its members are not cargo.
CARGO_CLASS_SCORE_SPECIFIC = 1.0
CARGO_CLASS_SCORE_GENERAL = 0.6


# Phage marker profiles used to recognise prophage context around a shell hit.
#
# Why this exists: Andreas & Giessen (2021) Nat Commun 12:4748
# (doi:10.1038/s41467-021-25071-y) report that their Family 1 and Family 2
# encapsulin datasets "were significantly contaminated with bacteriophage capsid
# proteins", and they removed hits whose genome neighborhoods contained
# phage-associated Pfam domains. Encapsulin shells share the HK97 fold with phage
# capsids (Pfam clan CL0373), so a sequence-only shell search cannot separate the
# two; the surrounding genes can.
#
# That paper lists the Pfam domains it filtered on as
# "PF0860, PF03354, PF04586, PF00589, PF05135". Four of those are used below.
# "PF0860" is not a well-formed Pfam accession (they have five digits), and the
# nearest valid accession, PF00860, is Xan_ur_permease, a xanthine/uracil
# permease family with no phage association. We could not determine what was
# intended, so it is deliberately omitted rather than guessed at. Flagged for
# review: see docs/methods.md.
#
# Accessions and names were checked against the InterPro API on 2026-11-26
# (https://www.ebi.ac.uk/interpro/api/entry/pfam/<accession>).
PHAGE_MARKER_PROFILES: List[ShellProfileDef] = [
    ShellProfileDef(
        family="Phage marker",
        pfam_id="PF03354",
        description="TerL_ATPase: terminase large subunit, ATPase domain",
        default_evalue=1e-5,
    ),
    ShellProfileDef(
        family="Phage marker",
        pfam_id="PF04586",
        description="Peptidase_S78: caudovirus prohead serine protease",
        default_evalue=1e-5,
    ),
    ShellProfileDef(
        family="Phage marker",
        pfam_id="PF00589",
        description="Phage_integrase: phage integrase family",
        default_evalue=1e-5,
    ),
    ShellProfileDef(
        family="Phage marker",
        pfam_id="PF05135",
        description="Phage_connect_1: phage gp6-like head-tail connector protein",
        default_evalue=1e-5,
    ),
]

# Subdirectory of the profile directory that holds the phage marker profiles, so
# that find_profile_files() does not pick them up as shell profiles.
PHAGE_MARKER_SUBDIR = "phage_markers"

# A shell locus is flagged as probable prophage context at or above this many
# distinct phage marker families in its window. One marker is reported but not
# flagged, because single phage-like genes are common in bacterial genomes
# (e.g. a lone integrase next to an unrelated gene).
PHAGE_MARKER_FLAG_THRESHOLD = 2


#: A syntactically valid Pfam accession: PF followed by five digits.
PFAM_ACCESSION_RE = re.compile(r"^PF\d{5}$")


def select_profiles(extra: Sequence[str] = ()) -> List[ShellProfileDef]:
    """Return the default profiles plus the requested extra profiles.

    Accessions in EXTENDED_SHELL_PROFILES come with a curated description. Any
    other syntactically valid Pfam accession is accepted too and fetched as-is,
    so that families cargofinder does not know about can still be scanned. Such
    a profile carries no claim that it detects an encapsulin: that is the
    caller's judgement, and the run summary records what was used.

    Note that Family 2B and Family 3 encapsulins have no Pfam entry at all, so
    covering them needs a custom HMM supplied through `--profiles`, not this.
    """
    profiles = list(DEFAULT_SHELL_PROFILES)
    for acc in extra:
        known = EXTENDED_SHELL_PROFILES.get(acc)
        if known is not None:
            profiles.append(known)
            continue
        if not PFAM_ACCESSION_RE.match(acc):
            raise ValueError(
                f"{acc!r} is not a Pfam accession (expected PF followed by five "
                f"digits, e.g. PF19307). Curated choices: "
                f"{sorted(EXTENDED_SHELL_PROFILES)}. To use a profile that is "
                f"not in Pfam, put the .hmm file in a directory and pass "
                f"--profiles <dir>."
            )
        profiles.append(
            ShellProfileDef(
                family=f"User-requested ({acc})",
                pfam_id=acc,
                description=(
                    f"{acc}, requested by the user. cargofinder makes no claim "
                    f"that this family detects encapsulin shells."
                ),
            )
        )
    return profiles


@dataclass
class ScoringWeights:
    """Additive scoring component weights.

    Score = w_distance * S_distance
          + w_operon * S_operon
          + w_targeting * S_targeting
          + w_cargo_class * S_cargo_class
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
    - Cargo class: Whether the candidate is a known encapsulin cargo class, by Pfam
      profile. Andreas & Giessen (2021) assigned cargo by identity this way, and it
      is the strongest single signal; but broad enzyme families contain many
      non-cargo proteins, and genuinely novel cargo scores zero, so it does not
      dominate. See cargofinder.cargo.
    - Recurrence: The same gene family appearing next to encapsulins across multiple
      independent loci/genomes provides strong evolutionary support.
    """

    w_distance: float = 0.25
    w_operon: float = 0.15
    w_targeting: float = 0.15
    w_cargo_class: float = 0.30
    w_tail_features: float = 0.05
    w_recurrence: float = 0.10

    def validate(self) -> None:
        """Validate that weights are non-negative and sum to 1.0."""
        weights = [
            self.w_distance,
            self.w_operon,
            self.w_targeting,
            self.w_cargo_class,
            self.w_tail_features,
            self.w_recurrence,
        ]
        if any(w < 0 for w in weights):
            raise ValueError(f"Weights must be non-negative: {self}")
        total = sum(weights)
        if not (0.99 <= total <= 1.01):
            raise ValueError(f"Weights must sum to 1.0 (got {total:.4f}): {self}")

    def without_cargo_class(self) -> "ScoringWeights":
        """Return these weights with the cargo class term removed and rescaled.

        Used when no cargo class profiles are available: rather than scoring
        every candidate 0 for that component, which would depress all totals and
        make runs incomparable, the remaining weights are rescaled to sum to 1.
        """
        remaining = (
            self.w_distance
            + self.w_operon
            + self.w_targeting
            + self.w_tail_features
            + self.w_recurrence
        )
        if remaining <= 0:
            raise ValueError(
                "Cannot drop the cargo class term: all other weights are zero."
            )
        factor = 1.0 / remaining
        return ScoringWeights(
            w_distance=self.w_distance * factor,
            w_operon=self.w_operon * factor,
            w_targeting=self.w_targeting * factor,
            w_cargo_class=0.0,
            w_tail_features=self.w_tail_features * factor,
            w_recurrence=self.w_recurrence * factor,
        )


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
