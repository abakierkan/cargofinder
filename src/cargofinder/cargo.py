"""Recognise known encapsulin cargo classes among neighborhood candidates.

Why this module exists
----------------------
v0.1 ranked neighbours by position and terminal-peptide properties alone, and so
had no notion of *what a gene is*. That is not how cargo is assigned in the
literature. Andreas & Giessen (2021) Nat Commun 12:4748,
doi:10.1038/s41467-021-25071-y, built "HMM templates representative of each
cargo class" and classified hits "based on their Pfam or Interpro annotation".
Identity is the strongest single piece of evidence they used, and it was the
main thing cargofinder did not do.

This module scans candidate proteins against the Pfam families that mark the
cargo classes that paper names (see `config.CARGO_CLASS_PROFILES`) and reports,
for each candidate, which class it matched.

What a match does and does not mean
-----------------------------------
Two of the families, PF22277 (EncFtn-like) and PF24309 (IMEF cargo Flp), are
defined by Pfam as encapsulin cargo proteins. A hit there is strong evidence.

Every other family is a broad enzyme family. A cysteine desulfurase domain
(PF00266) or a ferritin-like domain (PF00210) occurs throughout bacterial
genomes in proteins that have nothing to do with encapsulins. A hit says "this
is the right *kind* of enzyme to be cargo", which is a prior, not proof. The two
tiers are scored differently for exactly this reason, and the matched class is
printed on every output row so the distinction stays visible.

A corresponding blind spot: this can only recognise cargo classes that are
already known. Genuinely novel cargo scores zero here, so a candidate with no
class match is not thereby excluded. That is the main reason the component does
not dominate the score.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence

import pyhmmer
from pyhmmer import easel

from cargofinder.config import (
    CARGO_CLASS_PROFILES,
    CARGO_CLASS_SCORE_GENERAL,
    CARGO_CLASS_SCORE_SPECIFIC,
    CARGO_CLASS_SUBDIR,
    CargoClassDef,
)
from cargofinder.io import Gene
from cargofinder.shells import _as_text, load_hmm_profiles

logger = logging.getLogger(__name__)

#: Maps Pfam accession (without version) to its cargo class definition.
_BY_ACCESSION: Dict[str, CargoClassDef] = {d.pfam_id: d for d in CARGO_CLASS_PROFILES}
#: Maps the Pfam short name as it appears in the HMM NAME field.
_PROFILE_NAME_HINT: Dict[str, CargoClassDef] = {}


@dataclass(frozen=True)
class CargoClassHit:
    """A candidate protein matched to a known cargo class."""

    gene_id: str
    cargo_class: str
    pfam_id: str
    profile_name: str
    evalue: float
    encapsulin_specific: bool
    #: 1-based start of the earliest recognised domain on the protein. The
    #: residues before it are the unannotated N-terminal extension, which is
    #: what mediates encapsulation in Family 2 systems.
    first_domain_start: int = 1
    #: Full length of the candidate protein, for computing the extension.
    protein_length: int = 0

    @property
    def n_terminal_extension_aa(self) -> int:
        """Residues N-terminal to the earliest recognised domain."""
        return max(0, self.first_domain_start - 1)

    @property
    def score(self) -> float:
        """Score contribution in [0, 1]; see config for the two tiers."""
        return (
            CARGO_CLASS_SCORE_SPECIFIC
            if self.encapsulin_specific
            else CARGO_CLASS_SCORE_GENERAL
        )


def find_cargo_class_profiles(profile_dir: Path | str) -> List[Path]:
    """Locate cargo class profiles in the `cargo_classes` subdirectory.

    Returns an empty list when the subdirectory is absent, so a profile
    directory fetched before this feature existed keeps working: the run simply
    reports no cargo class evidence rather than failing.
    """
    sub = Path(profile_dir) / CARGO_CLASS_SUBDIR
    if not sub.is_dir():
        logger.warning(
            "No %s/ subdirectory in %s, so cargo class recognition is skipped and "
            "its score component is dropped. Run 'cargofinder fetch-profiles' "
            "again to download the cargo class profiles.",
            CARGO_CLASS_SUBDIR,
            profile_dir,
        )
        return []
    return sorted(list(sub.glob("*.hmm")) + list(sub.glob("*.hmm.gz")))


def _definition_for(profile_name: str, profile_acc: str) -> Optional[CargoClassDef]:
    """Resolve an HMM's accession or name back to its cargo class definition."""
    # Pfam accessions in HMM headers carry a version, e.g. "PF00266.23".
    bare = profile_acc.split(".")[0]
    if bare in _BY_ACCESSION:
        return _BY_ACCESSION[bare]
    return _PROFILE_NAME_HINT.get(profile_name)


def classify_cargo_candidates(
    genes: Sequence[Gene],
    cargo_profile_paths: Iterable[Path | str],
    evalue_cutoff: float = 1e-5,
) -> Dict[str, CargoClassHit]:
    """Assign each gene to a known cargo class, where one matches.

    All profiles are loaded once and every gene is scanned in a single
    hmmsearch pass. Where a gene matches several profiles, the best E-value
    wins, except that an encapsulin-specific family always outranks a broad one:
    a protein matching both PF22277 and the generic ferritin domain is reported
    as encapsulated ferritin, which is the more informative call.

    Parameters
    ----------
    genes : sequence of Gene
        Candidate genes to classify.
    cargo_profile_paths : iterable of Path or str
        Cargo class HMM files.
    evalue_cutoff : float
        Maximum full-sequence E-value for a class call.

    Returns
    -------
    dict of str -> CargoClassHit
        Genes with no cargo class match are absent from the mapping.
    """
    paths = [Path(p) for p in cargo_profile_paths]
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

    best: Dict[str, CargoClassHit] = {}
    # The N-terminal extension must be measured against the EARLIEST domain
    # found by ANY cargo class profile, not the best-scoring one. A DyP
    # peroxidase matches Dyp_perox_C around residue 140, but Dyp_perox_N starts
    # near residue 10, so using the best hit alone would report a spurious
    # ~140-residue "extension" for an ordinary two-domain protein.
    earliest: Dict[str, int] = {}
    lengths = {g.gene_id: len(g.protein_seq) for g in genes}

    for top_hits in pyhmmer.hmmer.hmmsearch(hmms, digital_seqs):
        profile_name = _as_text(getattr(top_hits.query, "name", ""))
        profile_acc = _as_text(getattr(top_hits.query, "accession", ""))
        definition = _definition_for(profile_name, profile_acc)
        if definition is None:
            logger.warning(
                "Profile %s (%s) in the cargo class directory is not a known "
                "cargo class and was ignored",
                profile_name,
                profile_acc or "no accession",
            )
            continue

        for hit in top_hits:
            if hit.evalue > evalue_cutoff:
                continue
            gene_id = _as_text(hit.name)
            starts = [
                d.alignment.target_from for d in hit.domains if d.alignment is not None
            ]
            cand = CargoClassHit(
                gene_id=gene_id,
                cargo_class=definition.cargo_class,
                pfam_id=definition.pfam_id,
                profile_name=profile_name,
                evalue=hit.evalue,
                encapsulin_specific=definition.encapsulin_specific,
                first_domain_start=min(starts) if starts else 1,
                protein_length=lengths.get(gene_id, 0),
            )
            prev = best.get(gene_id)
            if prev is None or _outranks(cand, prev):
                best[gene_id] = cand
            earliest[gene_id] = min(
                earliest.get(gene_id, cand.first_domain_start),
                cand.first_domain_start,
            )

    # Fold the across-profile earliest domain start back into each kept hit.
    for gene_id, hit in list(best.items()):
        first = earliest.get(gene_id, hit.first_domain_start)
        if first != hit.first_domain_start:
            best[gene_id] = replace(hit, first_domain_start=first)

    if best:
        counts: Dict[str, int] = {}
        for hit in best.values():
            counts[hit.cargo_class] = counts.get(hit.cargo_class, 0) + 1
        logger.info(
            "Recognised %d candidate(s) as known cargo classes: %s",
            len(best),
            ", ".join(f"{k} x{v}" for k, v in sorted(counts.items())),
        )
    return dict(sorted(best.items()))


def _outranks(new: CargoClassHit, old: CargoClassHit) -> bool:
    """An encapsulin-specific call beats a general one; otherwise best E-value."""
    if new.encapsulin_specific != old.encapsulin_specific:
        return new.encapsulin_specific
    return new.evalue < old.evalue
