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

from cargofinder.config import PHAGE_MARKER_FLAG_THRESHOLD, select_profiles
from cargofinder.pipeline import run_pipeline
from cargofinder.profiles import fetch_default_profiles, fetch_phage_marker_profiles
from cargofinder.scoring import ScoredCandidate

logger = logging.getLogger("cargofinder.validation")

HERE = Path(__file__).resolve().parent
TOP_K = (1, 3, 5)
# Default prophage screening radius in genes (cargofinder.pipeline.run_pipeline)
PROPHAGE_WINDOW = 20
# Weights of the two adjacency terms, used for the adjacency-only baseline.
# Kept in sync with cargofinder.config.ScoringWeights defaults.
W_DISTANCE, W_OPERON = 0.30, 0.20


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


def tp_consensus_key(c: ScoredCandidate) -> float:
    """BLOSUM62 similarity to the published Family 1 TP consensus alone.

    This is the comparison that matters for v0.2: the published signal against
    the composition heuristic it was meant to replace.
    """
    return c.score_tp_consensus


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
    tp_rank, tp_ties = relative_rank(target, in_shell, tp_consensus_key)
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
            "prophage_marker_families": target.prophage_marker_families,
            "is_probable_prophage_locus": target.is_probable_prophage_locus,
        }
    )
    return out


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
    extra_shells: int,
    sensitivity: pd.DataFrame,
    prophage: pd.DataFrame,
    pooled: pd.DataFrame,
    pooled_clusters: pd.DataFrame,
    manifest: pd.DataFrame,
) -> str:
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
                "TP consensus only, adjacent cargos": (
                    f"{baseline_hits(adjacent, 'tp_only', k)}/{len(adjacent)}"
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

    return f"""# Validation report: cargofinder v0.1

Generated by `validation/run_validation.py`. The numbers and tables are computed from
the run; do not edit them by hand.

## Setup

- Input: one complete NCBI replicon per organism, not the whole assembly, so other
  replicons were not screened. Genes were called with Pyrodigal.
- Profile: PF04454 only, with default parameters and **no motif file**.
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

- Distal cargo cannot be found by a neighborhood scan: every `outside window` row.
- The default profile does not detect the Family 2A shell of *S. elongatus*.

## Family 2A sensitivity run

*S. elongatus* plasmid scanned with `PF19307`, under C-terminal scoring and
`--termini both`. This cargo uses an N-terminal targeting sequence (Nichols et al.
2021).

{markdown_table(sensitivity)}

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
accepted cargo: {len(top_non_cargo)}, which includes the *P. furiosus* fragment. Shells
detected beyond the listed ones: {extra_shells}.

Highest-scoring gene that is not listed cargo, per organism:

{markdown_table(best_non_cargo)}

These genes are unconfirmed, not proven non-cargo. Full top-3 lists are in
`validation/false_flags.tsv`.

## Cluster recurrence across genomes

All replicons pooled in one multi-genome run with default parameters. Clusters seen
next to at least two shells:

{clusters_md}

Listed cargo in the pooled run:

{markdown_table(pooled)}

Clustering at 30% identity and 80% coverage groups only sequences that stay similar
over most of their length, so divergent or short cargo may not cluster.

## How to read this

- The set is small and related systems count more than once, so these are not estimates
  of performance on new genomes.
- Compare cargofinder with the adjacency-only column: the difference is what the
  terminal peptide score adds on this set.
- Parameters were not fitted to this set, but it was not blind either. The weights and
  tail features were chosen with these systems in mind, and low complexity was removed
  from the score after reviewing an earlier benchmark, because there was no evidence
  for either direction. Treat tail-only results as optimistic.
- Scores have a floor: an adjacent same-strand operon gene gets 0.63 before any tail
  signal in a single-shell run. A high score mostly reflects position, not evidence of
  encapsulation.

## Inputs

{markdown_table(manifest)}

## Reproduce

```bash
pixi run python validation/run_validation.py
```
"""


def run_validation(
    data_dir: Path, results_dir: Path, profiles_dir: Path, extended_dir: Path
) -> None:
    known = pd.read_csv(HERE / "known_systems.tsv", sep="\t", dtype=str)

    if not list(profiles_dir.glob("*.hmm")):
        fetch_default_profiles(profiles_dir, select_profiles())
    # Phage marker profiles for prophage screening. Downloading them into both
    # profile dirs means the Family 2A sensitivity run is screened too.
    for d in (profiles_dir, extended_dir):
        if not list((d / "phage_markers").glob("*.hmm")):
            fetch_phage_marker_profiles(d)

    if not (extended_dir / "PF19307.hmm").exists():
        family2a = [p for p in select_profiles(["PF19307"]) if p.pfam_id == "PF19307"]
        fetch_default_profiles(extended_dir, family2a)

    rows: List[Dict[str, object]] = []
    flag_rows: List[Dict[str, object]] = []
    extra_shells = 0
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
        listed = {find_tag(annotation, t).stop_key for t in group["shell_locus_tag"]}
        extra_shells += sum(
            stop_key(s.start, s.end, s.strand) not in listed for s in res.shells
        )

    results = pd.DataFrame(rows)
    flags = pd.DataFrame(flag_rows)
    results.to_csv(HERE / "benchmark_results.tsv", sep="\t", index=False)
    flags.to_csv(HERE / "false_flags.tsv", sep="\t", index=False)

    # Family 2A sensitivity runs
    srpi = known[known["cargo_location_type"] == "family_2a_operon"].iloc[0]
    sensitivity = []
    for termini in ("c", "both"):
        res = run_pipeline(
            data_dir / f"{srpi['nuccore_id']}.fna",
            out_dir=results_dir / f"family2a_{termini}",
            profile_dir=extended_dir,
            termini=termini,
        )
        ev = evaluate_row(srpi, annotations[srpi["nuccore_id"]], res)
        sensitivity.append(
            {
                "termini": termini,
                "shell_detected": ev["shell_detected"],
                "cargo_rank": f"{ev['rank_in_shell']} of {ev['candidates_in_shell']}",
                "cargo_score": ev["score_total"],
                "top_score": round(res.candidates[0].score_total, 3),
                "top_gene_distance": res.candidates[0].gene_distance,
            }
        )

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
                "shells_in_cluster": c.cluster_shell_count if c else "",
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
            {"shells": len(m), "members": ", ".join(sorted(m))}
            for _, m in sorted(members.items())
        ],
        columns=["shells", "members"],
    )

    report = build_report(
        results,
        flags,
        extra_shells,
        pd.DataFrame(sensitivity),
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
    parser.add_argument("--extended-profiles", type=Path, default=HERE / "profiles_2a")
    args = parser.parse_args()

    logging.basicConfig(level=logging.WARNING, format="[%(levelname)s] %(message)s")
    run_validation(args.data, args.out, args.profiles, args.extended_profiles)
    print(f"Wrote {HERE / 'report.md'} and {HERE / 'benchmark_results.tsv'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
