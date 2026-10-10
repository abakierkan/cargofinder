"""Unit tests for candidate reporting and run summary generation."""

import json

import pandas as pd

from cargofinder.config import RunConfig
from cargofinder.io import Gene, GenomeData
from cargofinder.neighborhood import NeighborCandidate
from cargofinder.reporting import (
    export_candidates_json,
    export_candidates_tsv,
    generate_run_summary,
    save_run_summary,
)
from cargofinder.scoring import ScoredCandidate
from cargofinder.shells import ShellHit


def _build_mock_scored_candidate(
    cand_id: str, rank: int, score: float
) -> ScoredCandidate:
    cand_gene = Gene(cand_id, "ctg1", 100, 250, 1, 0, "MVKTY")
    shell_gene = Gene("shell", "ctg1", 300, 500, 1, 1, "M" * 200)
    shell_hit = ShellHit(
        "shell", "Linocin", "PF04454", 1e-40, 200.0, 1e-40, 1, shell_gene
    )
    neighbor = NeighborCandidate(
        candidate_id=cand_id,
        shell_gene_id="shell",
        contig_id="ctg1",
        gene_distance=1,
        abs_gene_distance=1,
        strand_relative_distance=1,
        bp_distance=50,
        same_strand=True,
        is_likely_operon=True,
        candidate_gene=cand_gene,
        shell_hit=shell_hit,
    )
    from cargofinder.features import compute_peptide_features

    feats = compute_peptide_features("MVKTY")
    return ScoredCandidate(
        candidate_id=cand_id,
        shell_gene_id="shell",
        contig_id="ctg1",
        rank=rank,
        score_total=score,
        score_distance=1.0,
        score_operon=1.0,
        score_tail_features=0.5,
        score_recurrence=0.2,
        gene_distance=1,
        abs_gene_distance=1,
        strand_relative_distance=1,
        bp_distance=50,
        same_strand=True,
        is_likely_operon=True,
        cluster_id="rep1",
        cluster_size=1,
        cluster_shell_count=1,
        candidate=neighbor,
        features=feats,
    )


def test_export_tsv_and_json(tmp_path):
    c1 = _build_mock_scored_candidate("c1", 1, 0.95)
    c2 = _build_mock_scored_candidate("c2", 2, 0.75)

    tsv_path = tmp_path / "candidates.tsv"
    export_candidates_tsv([c1, c2], tsv_path)
    assert tsv_path.exists()

    df = pd.read_csv(tsv_path, sep="\t")
    assert len(df) == 2
    assert list(df["candidate_id"]) == ["c1", "c2"]
    assert "score_total" in df.columns
    assert "score_distance" in df.columns

    json_path = tmp_path / "candidates.json"
    export_candidates_json([c1, c2], json_path, metadata={"genome": "test_genome"})
    assert json_path.exists()
    payload = json.loads(json_path.read_text())
    assert "candidates" in payload
    assert len(payload["candidates"]) == 2
    assert payload["metadata"]["genome"] == "test_genome"


def test_run_summary_generation(tmp_path):
    c1 = _build_mock_scored_candidate("c1", 1, 0.95)
    shell_gene = Gene("shell", "ctg1", 300, 500, 1, 0, "M" * 200)
    shell_hit = ShellHit(
        "shell", "Linocin", "PF04454", 1e-40, 200.0, 1e-40, 1, shell_gene
    )
    gd = GenomeData("test_genome", [shell_gene, c1.candidate.candidate_gene])
    cfg = RunConfig()

    summary = generate_run_summary(gd, [shell_hit], [c1], cfg)
    assert summary["counts"]["shells_detected"] == 1
    assert summary["counts"]["candidates_extracted"] == 1
    assert "cargofinder_version" in summary
    assert "environment" in summary

    summary_file = tmp_path / "summary.json"
    save_run_summary(summary, summary_file)
    assert summary_file.exists()
    txt_file = tmp_path / "summary.txt"
    assert txt_file.exists()
    assert "Shells Detected:     1" in txt_file.read_text()
