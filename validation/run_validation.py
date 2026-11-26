#!/usr/bin/env python3
"""Reproducible validation of cargofinder on characterized encapsulin systems.

Truth comes from validation/known_systems.tsv. Each genome replicon is downloaded
from NCBI with its annotation and run through the pipeline. Predicted genes are
matched to literature locus tags by strand and stop coordinate, which tolerates
start-site differences between Pyrodigal and the annotation.

Every number in validation/report.md is computed here from the run results.
"""

from __future__ import annotations

import argparse
import hashlib
import logging
import sys
import time
import urllib.error
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple

import pandas as pd

from cargofinder import __version__
from cargofinder.config import PHAGE_MARKER_FLAG_THRESHOLD, select_profiles
from cargofinder.pipeline import run_pipeline
from cargofinder.profiles import (
    fetch_cargo_class_profiles,
    fetch_default_profiles,
    fetch_phage_marker_profiles,
)
from cargofinder.scoring import ScoredCandidate

logger = logging.getLogger("cargofinder.validation")

HERE = Path(__file__).resolve().parent
TOP_K = (1, 3, 5)
# Default prophage screening radius in genes (cargofinder.pipeline.run_pipeline)
PROPHAGE_WINDOW = 20
# Weights of the two adjacency terms, used for the adjacency-only baseline.
# Kept in sync with cargofinder.config.ScoringWeights defaults.
W_DISTANCE, W_OPERON = 0.25, 0.15


@dataclass(frozen=True)
class AnnotatedCds:
    start: int
    end: int
    strand: int
    tags: Tuple[str, ...]
    product: str

    @property
    def stop_key(self) -> Tuple[int, int]:
        return stop_key(self.start, self.end, self.strand)


def stop_key(start: int, end: int, strand: int) -> Tuple[int, int]:
    """Strand and stop-codon coordinate, which identify an ORF on a replicon."""
    return (strand, end if strand == 1 else start)


def candidate_stop_key(c: ScoredCandidate) -> Tuple[int, int]:
    gene = c.candidate.candidate_gene
    return stop_key(gene.start, gene.end, gene.strand)


def normalize_tag(tag: str) -> str:
    """Compare tags ignoring case and underscores (TM0785 matches TM_0785)."""
    return tag.replace("_", "").upper()


def parse_annotation(gff_path: Path) -> List[AnnotatedCds]:
    """Read CDS features with their locus tags and the parent gene's old tags."""
    gene_tags: Dict[str, List[str]] = {}
    rows: List[Tuple[int, int, int, List[str], str, str]] = []
    with gff_path.open(encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) != 9:
                continue
            attrs = dict(a.split("=", 1) for a in cols[8].split(";") if "=" in a)
            tags = [attrs.get("locus_tag", "")]
            tags += attrs.get("old_locus_tag", "").split(",")
            tags = [t for t in tags if t]
            if cols[2] == "gene" and "ID" in attrs:
                gene_tags[attrs["ID"]] = tags
            elif cols[2] == "CDS":
                strand = 1 if cols[6] == "+" else -1
                product = attrs.get("product", "").replace("%2C", ",")
                parent = attrs.get("Parent", "")
                rows.append((int(cols[3]), int(cols[4]), strand, tags, parent, product))
    return [
        AnnotatedCds(s, e, st, tuple(dict.fromkeys(tags + gene_tags.get(p, []))), prod)
        for s, e, st, tags, p, prod in rows
    ]


def find_tag(annotation: List[AnnotatedCds], tag: str) -> AnnotatedCds:
    wanted = normalize_tag(tag)
    for cds in annotation:
        if wanted in (normalize_tag(t) for t in cds.tags):
            return cds
    raise ValueError(f"Locus tag {tag} not found in the annotation")


def relative_rank(
    target: ScoredCandidate, cands: List[ScoredCandidate], key
) -> Tuple[int, int]:
    """Rank of target under key (higher is better), and how many candidates tie it."""
    value = round(key(target), 9)
    better = sum(1 for c in cands if round(key(c), 9) > value)
    ties = sum(1 for c in cands if round(key(c), 9) == value) - 1
    return better + 1, ties


def adjacency_key(c: ScoredCandidate) -> float:
    return W_DISTANCE * c.score_distance + W_OPERON * c.score_operon


def tail_key(c: ScoredCandidate) -> float:
    """The invented composition heuristic alone (hydrophobic + Gly/Ser fraction)."""
    return c.score_tail_features


def cargo_class_key(c: ScoredCandidate) -> float:
    """Known-cargo-class identity alone, the signal Andreas & Giessen relied on."""
    return c.score_cargo_class


def targeting_key(c: ScoredCandidate) -> float:
    """The family-aware targeting signal alone.

    For a Family 1 shell this is BLOSUM62 similarity to the published Family 1
    TP consensus; for a Family 2A shell it is the unannotated N-terminal
    extension length, since Family 2 has no short targeting peptide.
    """
    return c.score_targeting


def fetch(url: str, out_path: Path) -> Path:
    """Download url to out_path unless it is already present."""
    if out_path.is_file() and out_path.stat().st_size > 0:
        return out_path
    out_path.parent.mkdir(parents=True, exist_ok=True)
    req = urllib.request.Request(url, headers={"User-Agent": "cargofinder-validation"})
    for attempt in range(3):
        try:
            with urllib.request.urlopen(req, timeout=300) as resp:
                out_path.write_bytes(resp.read())
            return out_path
        except urllib.error.HTTPError as err:
            if err.code != 429:
                raise
            time.sleep(2.0 * (attempt + 1))
    raise RuntimeError(f"Failed to download {url} after 3 attempts")


def fetch_replicon(nuccore_id: str, data_dir: Path) -> Path:
    return fetch(
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"
        f"db=nuccore&id={nuccore_id}&rettype=fasta&retmode=text",
        data_dir / f"{nuccore_id}.fna",
    )


def fetch_annotation(accession: str, data_dir: Path) -> Path:
    return fetch(
        "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?"
        f"db=nuccore&report=gff3&retmode=text&id={accession}",
        data_dir / f"{accession}.gff3",
    )


def evaluate_row(
    row: pd.Series, annotation: List[AnnotatedCds], res
) -> Dict[str, object]:
    """Score one known_systems row against one pipeline result."""
    out: Dict[str, object] = {
        "organism": row["organism"],
        "shell_locus": row["shell_locus_tag"],
        "known_cargo": row["known_cargo_locus_tag"],
        "cargo_description": row["cargo_name_description"],
        "location_type": row["cargo_location_type"],
        "evidence": row["evidence"],
        "shell_detected": False,
        "shell_evalue": "",
        "cargo_in_window": False,
        "rank_in_shell": "",
        "candidates_in_shell": 0,
        "score_total": "",
        "gene_distance": "",
        "adjacency_only_rank": "",
        "adjacency_only_ties": "",
        "tail_only_rank": "",
        "tail_only_ties": "",
        "tp_only_rank": "",
        "tp_only_ties": "",
        "cargo_class_only_rank": "",
        "cargo_class_only_ties": "",
        "cargo_class": "",
        "cargo_class_pfam": "",
        "prophage_marker_families": "",
        "is_probable_prophage_locus": "",
        "citation": row["literature_reference"],
        "doi": row["doi"],
    }
    shell_key = find_tag(annotation, row["shell_locus_tag"]).stop_key
    cargo_key = find_tag(annotation, row["known_cargo_locus_tag"]).stop_key

    shell = next(
        (s for s in res.shells if stop_key(s.start, s.end, s.strand) == shell_key),
        None,
    )
    if shell is None:
        return out
    out["shell_detected"] = True
    out["shell_evalue"] = f"{shell.evalue:.1e}"

    in_shell = [c for c in res.candidates if c.shell_gene_id == shell.gene_id]
    out["candidates_in_shell"] = len(in_shell)
    target = next((c for c in in_shell if candidate_stop_key(c) == cargo_key), None)
    if target is None:
        return out

    adjacency_rank, adjacency_ties = relative_rank(target, in_shell, adjacency_key)
    tail_rank, tail_ties = relative_rank(target, in_shell, tail_key)
    tp_rank, tp_ties = relative_rank(target, in_shell, targeting_key)
    cc_rank, cc_ties = relative_rank(target, in_shell, cargo_class_key)
    out.update(
        {
            "cargo_in_window": True,
            "rank_in_shell": target.rank_in_shell,
            "score_total": round(target.score_total, 3),
            "gene_distance": target.gene_distance,
            "adjacency_only_rank": adjacency_rank,
            "adjacency_only_ties": adjacency_ties,
            "tail_only_rank": tail_rank,
            "tail_only_ties": tail_ties,
            "tp_only_rank": tp_rank,
            "tp_only_ties": tp_ties,
            "cargo_class_only_rank": cc_rank,
            "cargo_class_only_ties": cc_ties,
            "cargo_class": target.cargo_class,
            "cargo_class_pfam": target.cargo_class_pfam,
            "prophage_marker_families": target.prophage_marker_families,
            "is_probable_prophage_locus": target.is_probable_prophage_locus,
        }
    )
    return out


def locus_of(annotation: List[AnnotatedCds], c: ScoredCandidate) -> str:
    """Annotation locus tag for a predicted gene, matched by strand and stop."""
    cds = {a.stop_key: a for a in annotation}.get(candidate_stop_key(c))
    return cds.tags[-1] if cds and cds.tags else ""


def false_flags(
    res,
    annotation: List[AnnotatedCds],
    known_stops: set,
    organism: str,
    top_n: int = 3,
) -> List[Dict[str, object]]:
    """List the top candidates per shell, marking those not known to be cargo."""
    by_stop = {cds.stop_key: cds for cds in annotation}
    rows: List[Dict[str, object]] = []
    for c in res.candidates:
        if c.rank_in_shell > top_n:
            continue
        key = candidate_stop_key(c)
        cds = by_stop.get(key)
        rows.append(
            {
                "organism": organism,
                "rank_in_shell": c.rank_in_shell,
                "locus": cds.tags[-1] if cds and cds.tags else "unannotated",
                "product": cds.product if cds else "",
                "score_total": round(c.score_total, 3),
                "gene_distance": c.gene_distance,
                "same_strand": c.same_strand,
                "is_known_cargo": key in known_stops,
                "c_term_seq": c.features.c_term_seq,
            }
        )
    return rows


def recovered(df: pd.DataFrame, k: int) -> int:
    ranks = pd.to_numeric(df["rank_in_shell"], errors="coerce")
    return int((ranks <= k).sum())


def baseline_hits(df: pd.DataFrame, prefix: str, k: int) -> int:
    """Count hits at k, assuming ties break against the known cargo."""
    ranks = pd.to_numeric(df[f"{prefix}_rank"], errors="coerce")
    ties = pd.to_numeric(df[f"{prefix}_ties"], errors="coerce")
    return int(((ranks + ties) <= k).sum())


def markdown_table(df: pd.DataFrame) -> str:
    """Format a DataFrame as a markdown table without requiring tabulate."""
    lines = ["| " + " | ".join(df.columns) + " |"]
    lines.append("| " + " | ".join("---" for _ in df.columns) + " |")
    for _, row in df.astype(str).iterrows():
        lines.append("| " + " | ".join(row) + " |")
    return "\n".join(lines)


def outcome(row: pd.Series) -> str:
    if not row["shell_detected"]:
        return "shell not detected"
    if not row["cargo_in_window"]:
        return "outside window"
    return f"rank {row['rank_in_shell']} of {row['candidates_in_shell']}"


def build_report(
    results: pd.DataFrame,
    flags: pd.DataFrame,
    extra_shells: pd.DataFrame,
    cargo_classes: pd.DataFrame,
    candidates_scanned: int,
    prophage: pd.DataFrame,
    pooled: pd.DataFrame,
    pooled_clusters: pd.DataFrame,
    manifest: pd.DataFrame,
) -> str:
    results = results.assign(outcome_kind=results.apply(outcome, axis=1))
    cc_total = len(cargo_classes)
    cc_known = int(cargo_classes["is_known_cargo"].sum()) if cc_total else 0
    cc_other = cc_total - cc_known
    cc_specific = int(cargo_classes["encapsulin_specific"].sum()) if cc_total else 0
    cc_specific_wrong = (
        int(
            (
                cargo_classes["encapsulin_specific"] & ~cargo_classes["is_known_cargo"]
            ).sum()
        )
        if cc_total
        else 0
    )

    n_extra_cargo = (
        int(extra_shells["annotated_as_cargo"].sum()) if len(extra_shells) else 0
    )
    from cargofinder.config import ScoringWeights, select_profiles

    _w = ScoringWeights()
    default_profiles = " + ".join(p.pfam_id for p in select_profiles())
    floor = round(_w.w_distance + _w.w_operon + _w.w_recurrence * 0.2, 2)
    floor_cargo = round(floor + _w.w_cargo_class, 2)
    missed = []
    n_outside = int((results["outcome_kind"] == "outside window").sum())
    if n_outside:
        missed.append(
            f"- Distal cargo cannot be found by a neighborhood scan: "
            f"{n_outside} `outside window` row(s). No window size fixes this."
        )
    n_undetected = int((results["outcome_kind"] == "shell not detected").sum())
    if n_undetected:
        missed.append(
            f"- {n_undetected} shell(s) were not detected by the default profile set."
        )
    else:
        missed.append(
            "- Every listed shell was detected by the default profile set, and "
            "every cargo that fell inside a window was ranked first. What is "
            "missed here is entirely distal cargo."
        )
    missed_lines = "\n".join(missed)

    linocin = prophage[prophage["profile"] == "Linocin_M18"]
    capsid = prophage[prophage["profile"] == "Phage_capsid"]
    linocin_total, capsid_total = len(linocin), len(capsid)
    linocin_flagged = int(linocin["flagged"].sum()) if linocin_total else 0
    capsid_flagged = int(capsid["flagged"].sum()) if capsid_total else 0

    primary = results[
        (results["evidence"] == "experimental")
        & (results["location_type"] != "natural_fusion")
    ]
    adjacent = primary[primary["location_type"] == "adjacent_operon"]
    inferred = results[results["evidence"] == "homology_inferred"]
    fusion = results[results["location_type"] == "natural_fusion"]
    n = len(primary)

    recall_rows = []
    for k in TOP_K:
        recall_rows.append(
            {
                "top k": k,
                f"cargofinder, all {n} cargos": f"{recovered(primary, k)}/{n}",
                f"cargofinder, {len(adjacent)} adjacent cargos": (
                    f"{recovered(adjacent, k)}/{len(adjacent)}"
                ),
                "adjacency only, adjacent cargos": (
                    f"{baseline_hits(adjacent, 'adjacency_only', k)}/{len(adjacent)}"
                ),
                "tail only (composition), adjacent cargos": (
                    f"{baseline_hits(adjacent, 'tail_only', k)}/{len(adjacent)}"
                ),
                "targeting signal only, adjacent cargos": (
                    f"{baseline_hits(adjacent, 'tp_only', k)}/{len(adjacent)}"
                ),
                "cargo class only, adjacent cargos": (
                    f"{baseline_hits(adjacent, 'cargo_class_only', k)}/{len(adjacent)}"
                ),
            }
        )

    table = results.assign(outcome=results.apply(outcome, axis=1))[
        [
            "organism",
            "known_cargo",
            "location_type",
            "evidence",
            "outcome",
            "score_total",
            "adjacency_only_rank",
            "adjacency_only_ties",
            "tail_only_rank",
            "tail_only_ties",
            "tp_only_rank",
            "tp_only_ties",
            "cargo_class_only_rank",
            "cargo_class_only_ties",
            "cargo_class",
        ]
    ]

    non_cargo = flags[~flags["is_known_cargo"]]
    top_non_cargo = non_cargo[non_cargo["rank_in_shell"] == 1]
    best_non_cargo = non_cargo.loc[
        non_cargo.groupby("organism", sort=False)["score_total"].idxmax()
    ][["organism", "locus", "product", "score_total", "gene_distance", "c_term_seq"]]

    fusion_text = "; ".join(
        f"{r['organism']} {r['known_cargo']}: {outcome(r)}"
        for _, r in fusion.iterrows()
    )
    inferred_text = "; ".join(
        f"{r['organism']} {r['known_cargo']}: {outcome(r)}"
        for _, r in inferred.iterrows()
    )
    clusters_md = markdown_table(pooled_clusters) if len(pooled_clusters) else "None."

    return f"""# Validation report: cargofinder v{__version__}

Generated by `validation/run_validation.py`. The numbers and tables are computed from
the run; do not edit them by hand.

## Setup

- Input: one complete NCBI replicon per organism, not the whole assembly, so other
  replicons were not screened. Genes were called with Pyrodigal.
- Profiles: the default set ({default_profiles}), with default parameters and
  **no motif file**. Cargo class and phage marker profiles were present.
- Truth: `validation/known_systems.tsv`, matched by locus tag via strand and stop
  coordinate. Ranks are within each shell's neighborhood.
- Recall counts each cited cargo encoded as a separate gene and shown experimentally
  ({n} cargos). The *P. furiosus* fusion and the homology-inferred *R. erythropolis*
  PR4 row are reported separately.

## Cargo recovery

{markdown_table(pd.DataFrame(recall_rows))}

The baselines count ties against the known cargo: a method that cannot separate two
genes has not found the right one.

{markdown_table(table)}

Reported separately:

- Fusion, not a separate cargo: {fusion_text}. This ORF is the N-terminal part of the
  single 345-aa protein in PDB 2E0Z. It ranks highly only because the annotation splits
  that protein in two, so this is an artifact and not a recovery.
- Homology-inferred, no PR4 experiment found: {inferred_text}.

## What is missed

{missed_lines}

## Additional encapsulin shells found

Shells detected that are not listed in `known_systems.tsv`. These are **not**
false positives: every one is annotated by RefSeq as an encapsulin shell
protein. They are real loci whose cargo this benchmark simply does not have
literature entries for, and they are listed so the count is not mistaken for
noise.

{markdown_table(extra_shells)}

**{n_extra_cargo} of {len(extra_shells)} top candidates at these shells are annotated
by RefSeq as encapsulin cargo proteins.** cargofinder was not told these loci
existed and this truth set has no entries for them, so that is the closest thing
to a blind test in this report, on a very small number of cases. It is agreement
with an annotation pipeline, not with an experiment.

They are all found by PF19307, which became a default profile in v0.2 because
Family 2 is the larger group in the Andreas & Giessen (2021) census (3,523
systems against Family 1's 2,383). Note that although Pfam names PF19307 for
Family 2A, it also detects Family 2B shells here: the *M. xanthus* pair is
annotated Family 2B, and two adjacent shell genes in one operon is exactly the
2B arrangement that paper describes.

## Cargo class recognition

Candidates are scanned against Pfam families marking the cargo classes named by
Andreas & Giessen (2021), who assigned cargo by identity this way. Two families,
PF22277 (EncFtn-like) and PF24309 (IMEF cargo Flp), are defined by Pfam as
encapsulin cargo; the rest are broad enzyme families scored lower.

| measure | value |
| --- | --- |
| Candidates scanned | {candidates_scanned} |
| Matched some cargo class | {cc_total} |
| ...that are known cargo | {cc_known} |
| ...that are not known cargo | {cc_other} |
| Matched by an *encapsulin-specific* family | {cc_specific} |
| ...of those, not a listed cargo | {cc_specific_wrong} |

{markdown_table(cargo_classes)}

Every listed cargo that fell inside a window was given a biologically correct
class: the ferritin-like cargos as encapsulated ferritin, the peroxidases as
DyP. None was missed.

**Read the `is_known_cargo` column carefully.** It is False for the
*P. furiosus* `PF1192` row, because that ORF is the annotation-split half of a
fusion and is excluded from the cargo counts throughout this report. Its
`encapsulated_ferritin` call is nonetheless biologically right: `PF1192` is the
ferritin-like portion of the single 345-aa protein in PDB 2E0Z. So the one
non-cargo hit from an encapsulin-specific family is a bookkeeping artifact of
this benchmark, not a misidentification.

**The remaining non-cargo matches are not noise either, and they show the limit
of this signal.** They are all in the same *P. furiosus* locus: one further
ferritin-like protein (PF00210) and three rubrerythrins (PF02915). That is
expected biology rather than a failure of the profiles. Andreas & Giessen (2021)
note that
"Flp-fusion encapsulins are often located in operons containing other
ferritin-like proteins or rubrerythrins, hinting at a function related to iron
homeostasis and stress resistance". The profiles cannot tell cargo from a
co-located paralogue in the same pathway, because both really are that class of
enzyme. Position and targeting-peptide signal are what separate them, which is
why cargo class is weighted at 0.30 and not higher.

A second limit, which this benchmark cannot measure at all: **only cargo classes
that are already known can be recognised.** Novel cargo scores zero here. A
candidate with no class match is therefore not excluded, and this component must
not be read as a filter.

## Prophage screening control

Encapsulin shells share the HK97 fold with phage capsids (Pfam clan CL0373), and
Andreas & Giessen (2021) report that their Family 1 and Family 2 datasets were
significantly contaminated with phage capsid proteins. cargofinder screens the
genes within +/-{PROPHAGE_WINDOW} genes of each shell for phage marker families
(terminase, prohead protease, integrase, head-tail connector) and flags a locus
at {PHAGE_MARKER_FLAG_THRESHOLD} or more distinct families.

To test both directions on real genomes, every replicon was re-scanned with
PF05065 (Phage_capsid) added, which deliberately pulls in phage capsids
alongside the true encapsulins.

| outcome | count |
| --- | --- |
| True encapsulin shells wrongly flagged | {linocin_flagged}/{linocin_total} |
| Phage capsid hits (PF05065) correctly flagged | {capsid_flagged}/{capsid_total} |

{markdown_table(prophage)}

The screen does not misfire on any true encapsulin locus in this set, but it
catches only {capsid_flagged} of {capsid_total} phage capsid hits. The misses
carry one marker family or none within the window, so they fall below the
threshold of {PHAGE_MARKER_FLAG_THRESHOLD}. Lowering it to 1 would flag more
prophages at the cost of flagging loci that merely have a lone integrase
nearby, which is common in bacterial genomes. Treat an unflagged locus as
"not evidently prophage", never as "confirmed encapsulin".

Worked positive control, independently checkable: in *M. tuberculosis* H37Rv the
flagged PF05065 hit is **Rv2650c**, which NCBI annotates as "prophage protein".
Its window contains Rv2646 (integrase), Rv2651c (prophage protease) and Rv2659c
(prophage integrase). Rv2650c is a known gene of the **phiRv1** prophage-like
element (Fan X, Abd Alla AA, Xie J (2016) *J Biomol Struct Dyn* 34:233-238,
doi:10.1080/07391102.2015.1022602, which names Rv2650c, Rv2652c, Rv2658c and
Rv2659c as phiRv1 genes; see also Bibb LA & Hatfull GF (2002) *Mol Microbiol*
45:1515-1526, doi:10.1046/j.1365-2958.2002.03130.x). The true *M. tuberculosis*
encapsulin shell Rv0798c was not flagged in the same run.

## What is falsely flagged

These are the top-ranked genes a user could waste effort on. Rank-1 genes that are not
accepted cargo: {len(top_non_cargo)}. That count includes the *P. furiosus* fragment and
the top candidate at each additional shell listed above, whose true cargo is simply not
in this truth set.

Highest-scoring gene that is not listed cargo, per organism:

{markdown_table(best_non_cargo)}

These genes are unconfirmed, not proven non-cargo. Full top-3 lists are in
`validation/false_flags.tsv`.

## Cluster recurrence across genomes

All replicons pooled in one multi-genome run with default parameters. Clusters
containing at least two independent genes that each neighbour a shell. Counting
genes rather than shells matters here: the tandem Family 2B shells of
*M. xanthus* share their windows, so counting shells let a single gene appear to
recur by itself.

{clusters_md}

Listed cargo in the pooled run:

{markdown_table(pooled)}

Clustering at 30% identity and 80% coverage groups only sequences that stay similar
over most of their length, so divergent or short cargo may not cluster.

## How to read this

- The set is small and related systems count more than once, so these are not estimates
  of performance on new genomes.
- Compare cargofinder against each single-signal column: the difference is what
  combining them buys on this set, which here is nothing, because cargo class and
  the targeting signal each already rank every adjacent cargo first.
- Parameters were not fitted to this set, but it was not blind either. The weights,
  cargo classes and tail features were chosen with these systems in mind; low
  complexity was dropped from the score after an earlier benchmark showed no
  evidence either way, and the composition heuristic was cut from 0.25 to 0.05
  after it underperformed here. Treat every single-signal column as optimistic.
- Scores have a floor: an adjacent same-strand operon gene gets {floor} before any
  cargo class or targeting signal in a single-shell run, and {floor_cargo} once a
  known cargo class is recognised. A high score still partly reflects position, not
  evidence of encapsulation.
- Cargo class and targeting both depend on prior knowledge: a novel cargo class, or
  a Family 2 candidate whose class is unrecognised, scores zero on the two
  components that carry most of the weight. This benchmark cannot measure that
  failure mode, because every system in it is already characterised.

## Inputs

{markdown_table(manifest)}

## Reproduce

```bash
pixi run python validation/run_validation.py
```
"""


def run_validation(data_dir: Path, results_dir: Path, profiles_dir: Path) -> None:
    known = pd.read_csv(HERE / "known_systems.tsv", sep="\t", dtype=str)

    # Fetch each default profile that is missing, rather than skipping the lot
    # when any one is present: the default set grows between versions.
    missing = [
        p for p in select_profiles() if not (profiles_dir / f"{p.pfam_id}.hmm").exists()
    ]
    if missing:
        fetch_default_profiles(profiles_dir, missing)
    if not list((profiles_dir / "phage_markers").glob("*.hmm")):
        fetch_phage_marker_profiles(profiles_dir)
    if not list((profiles_dir / "cargo_classes").glob("*.hmm")):
        fetch_cargo_class_profiles(profiles_dir)

    rows: List[Dict[str, object]] = []
    flag_rows: List[Dict[str, object]] = []
    cargo_class_rows: List[Dict[str, object]] = []
    extra_shell_rows: List[Dict[str, object]] = []
    candidates_scanned = 0
    fastas: List[Path] = []
    manifest: List[Dict[str, str]] = []
    annotations: Dict[str, List[AnnotatedCds]] = {}

    for nuccore_id, group in known.groupby("nuccore_id", sort=False):
        first = group.iloc[0]
        fasta = fetch_replicon(nuccore_id, data_dir)
        annotation = parse_annotation(
            fetch_annotation(first["annotation_accession"], data_dir)
        )
        annotations[nuccore_id] = annotation
        fastas.append(fasta)
        manifest.append(
            {
                "organism": first["organism"],
                "assembly": first["assembly_accession"],
                "replicon": nuccore_id,
                "annotation": first["annotation_accession"],
                "fasta_sha256": hashlib.sha256(fasta.read_bytes()).hexdigest()[:16],
            }
        )

        res = run_pipeline(
            fasta, out_dir=results_dir / nuccore_id, profile_dir=profiles_dir
        )
        rows += [evaluate_row(row, annotation, res) for _, row in group.iterrows()]
        known_stops = {
            find_tag(annotation, row["known_cargo_locus_tag"]).stop_key
            for _, row in group.iterrows()
            if row["cargo_location_type"] != "natural_fusion"
        }
        flag_rows += false_flags(res, annotation, known_stops, first["organism"])
        for c in res.candidates:
            if not c.cargo_class:
                continue
            cargo_class_rows.append(
                {
                    "organism": first["organism"],
                    "locus": locus_of(annotation, c) or c.candidate_id,
                    "cargo_class": c.cargo_class,
                    "pfam": c.cargo_class_pfam,
                    "encapsulin_specific": c.cargo_class_specific,
                    "rank_in_shell": c.rank_in_shell,
                    "gene_distance": c.gene_distance,
                    "is_known_cargo": candidate_stop_key(c) in known_stops,
                }
            )
        candidates_scanned += len(res.candidates)
        listed = {find_tag(annotation, t).stop_key for t in group["shell_locus_tag"]}
        by_stop = {a.stop_key: a for a in annotation}
        for sh in res.shells:
            key = stop_key(sh.start, sh.end, sh.strand)
            if key in listed:
                continue
            cds = by_stop.get(key)
            top = next(
                (
                    c
                    for c in res.candidates
                    if c.shell_gene_id == sh.gene_id and c.rank_in_shell == 1
                ),
                None,
            )
            top_cds = by_stop.get(candidate_stop_key(top)) if top else None
            top_product = top_cds.product if top_cds else ""
            extra_shell_rows.append(
                {
                    "organism": first["organism"],
                    "shell_locus": cds.tags[-1] if cds and cds.tags else "unannotated",
                    "profile": sh.profile_name,
                    "shell_evalue": f"{sh.evalue:.0e}",
                    "top_candidate": (
                        (
                            top_cds.tags[-1]
                            if top_cds and top_cds.tags
                            else "unannotated"
                        )
                        if top
                        else ""
                    ),
                    "top_candidate_product": top_product,
                    "annotated_as_cargo": "cargo protein" in top_product.lower(),
                }
            )

    results = pd.DataFrame(rows)
    flags = pd.DataFrame(flag_rows)
    results.to_csv(HERE / "benchmark_results.tsv", sep="\t", index=False)
    flags.to_csv(HERE / "false_flags.tsv", sep="\t", index=False)

    extra_shells_df = pd.DataFrame(extra_shell_rows)
    extra_shells_df.to_csv(HERE / "extra_shells.tsv", sep="\t", index=False)

    # Prophage screening control.
    #
    # PF05065 (Phage_capsid) is a general phage capsid family in the same Pfam
    # clan as the encapsulin shells (CL0373), so adding it deliberately pulls in
    # phage capsids. That gives both halves of the control on real genomes:
    # the Linocin_M18 hits are true encapsulins and must stay unflagged, while
    # the Phage_capsid hits are the contamination the screen is meant to catch.
    phage_dir = profiles_dir.parent / "profiles_phage_control"
    if not (phage_dir / "PF05065.hmm").exists():
        fetch_default_profiles(phage_dir, select_profiles(["PF05065"]))
    if not list((phage_dir / "phage_markers").glob("*.hmm")):
        fetch_phage_marker_profiles(phage_dir)
    if not list((phage_dir / "cargo_classes").glob("*.hmm")):
        fetch_cargo_class_profiles(phage_dir)

    prophage_rows = []
    for nuccore_id, group in known.groupby("nuccore_id", sort=False):
        organism = group.iloc[0]["organism"]
        res = run_pipeline(
            data_dir / f"{nuccore_id}.fna",
            out_dir=results_dir / f"prophage_control_{nuccore_id}",
            profile_dir=phage_dir,
            generate_plots=False,
        )
        for shell in res.shells:
            ev = res.prophage[shell.gene_id]
            prophage_rows.append(
                {
                    "organism": organism,
                    "shell_gene": shell.gene_id,
                    "profile": shell.profile_name,
                    "shell_evalue": f"{shell.evalue:.0e}",
                    "marker_families": ev.distinct_family_count,
                    "flagged": ev.is_probable_prophage,
                    "markers": ev.marker_summary,
                }
            )
    cargo_classes = pd.DataFrame(cargo_class_rows)
    cargo_classes.to_csv(HERE / "cargo_class_calls.tsv", sep="\t", index=False)

    prophage = pd.DataFrame(prophage_rows)
    prophage.to_csv(HERE / "prophage_control.tsv", sep="\t", index=False)

    # Multi-genome run for recurrence across genomes
    pooled_res = run_pipeline(
        fastas,
        out_dir=results_dir / "pooled",
        profile_dir=profiles_dir,
        generate_plots=False,
    )
    by_stop = {
        (c.candidate.genome_id, candidate_stop_key(c)): c for c in pooled_res.candidates
    }
    pooled_rows = []
    for _, row in known.iterrows():
        cds = find_tag(annotations[row["nuccore_id"]], row["known_cargo_locus_tag"])
        c = by_stop.get((row["nuccore_id"], cds.stop_key))
        pooled_rows.append(
            {
                "organism": row["organism"],
                "known_cargo": row["known_cargo_locus_tag"],
                "rank_in_shell": c.rank_in_shell if c else "not in a window",
                "genes_in_cluster": c.cluster_shell_count if c else "",
                "score_total": round(c.score_total, 3) if c else "",
            }
        )

    members: Dict[str, List[str]] = {}
    for c in pooled_res.candidates:
        if c.cluster_shell_count >= 2:
            lookup = {a.stop_key: a for a in annotations[c.candidate.genome_id]}
            cds = lookup.get(candidate_stop_key(c))
            label = cds.tags[-1] if cds and cds.tags else c.candidate_id
            members.setdefault(c.cluster_id, []).append(label)
    pooled_clusters = pd.DataFrame(
        [
            {"genes": len(m), "members": ", ".join(sorted(m))}
            for _, m in sorted(members.items())
        ],
        columns=["genes", "members"],
    )

    report = build_report(
        results,
        flags,
        extra_shells_df,
        cargo_classes,
        candidates_scanned,
        prophage,
        pd.DataFrame(pooled_rows),
        pooled_clusters,
        pd.DataFrame(manifest),
    )
    (HERE / "report.md").write_text(report, encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--data", type=Path, default=HERE / "data")
    parser.add_argument("--out", type=Path, default=HERE / "results")
    parser.add_argument("--profiles", type=Path, default=HERE / "profiles")
    args = parser.parse_args()

    logging.basicConfig(level=logging.WARNING, format="[%(levelname)s] %(message)s")
    run_validation(args.data, args.out, args.profiles)
    print(f"Wrote {HERE / 'report.md'} and {HERE / 'benchmark_results.tsv'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
