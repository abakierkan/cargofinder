"""Unit tests for transparent additive candidate scoring."""

import math
from pathlib import Path

from cargofinder.cluster import ClusterResult
from cargofinder.config import ScoringWeights
from cargofinder.features import compute_peptide_features, load_motif_patterns
from cargofinder.io import Gene, load_protein_and_gff3
from cargofinder.neighborhood import NeighborCandidate, extract_neighborhood
from cargofinder.scoring import (
    compute_distance_score,
    compute_operon_score,
    compute_recurrence_score,
    compute_tail_features_score,
    score_candidates,
    score_termini,
)
from cargofinder.shells import ShellHit, detect_shells

FIXTURES_DIR = Path(__file__).parent / "data"
PROFILE_PATH = FIXTURES_DIR / "PF04454_test.hmm"
MOTIF_FILE = Path(__file__).parent.parent / "src/cargofinder/data/example_motifs.json"


def test_hand_computed_score_components():
    # Distance scores (window 5)
    assert compute_distance_score(1, 5) == 1.0
    assert compute_distance_score(2, 5) == 0.8
    assert compute_distance_score(3, 5) == 0.6
    assert compute_distance_score(4, 5) == 0.4
    assert compute_distance_score(5, 5) == 0.2
    assert compute_distance_score(6, 5) == 0.0

    # Operon scores
    assert compute_operon_score(same_strand=True, is_likely_operon=True) == 1.0
    assert compute_operon_score(same_strand=True, is_likely_operon=False) == 0.5
    assert compute_operon_score(same_strand=False, is_likely_operon=False) == 0.0

    # Recurrence scores
    assert math.isclose(compute_recurrence_score(1), 0.2)
    assert math.isclose(compute_recurrence_score(2), 0.6)
    assert math.isclose(compute_recurrence_score(3), 1.0)
    assert math.isclose(compute_recurrence_score(10), 1.0)


def _build_mock_candidate(
    cand_id: str,
    dist: int,
    same_strand: bool,
    is_operon: bool,
    seq: str,
) -> NeighborCandidate:
    cand_gene = Gene(cand_id, "ctg1", 100, 300, 1 if same_strand else -1, 0, seq)
    shell_gene = Gene("shell", "ctg1", 400, 600, 1, 1, "M" * 200)
    shell_hit = ShellHit(
        "shell", "Linocin", "PF04454", 1e-40, 200.0, 1e-40, 1, shell_gene
    )
    return NeighborCandidate(
        candidate_id=cand_id,
        shell_gene_id="shell",
        contig_id="ctg1",
        gene_distance=dist,
        abs_gene_distance=abs(dist),
        strand_relative_distance=dist,
        bp_distance=50,
        same_strand=same_strand,
        is_likely_operon=is_operon,
        candidate_gene=cand_gene,
        shell_hit=shell_hit,
    )


def test_toy_case_known_ranking():
    # Candidate A: Ideal cargo (adjacent +1, in operon, matching motif,
    # recurrent across 2 shells)
    cand_a = _build_mock_candidate(
        cand_id="cand_A",
        dist=1,
        same_strand=True,
        is_operon=True,
        seq="M" * 50 + "EEETSGGSENTGGDLGIRKL",
    )
    # Candidate B: Moderate (distance 3, same strand, no operon,
    # non-matching tail, 1 shell)
    cand_b = _build_mock_candidate(
        cand_id="cand_B",
        dist=3,
        same_strand=True,
        is_operon=False,
        seq="M" * 50 + "AAAAAAAAAAGGGGGGGGGG",
    )
    # Candidate C: Poor (distance 5, opposite strand, no operon, 1 shell)
    cand_c = _build_mock_candidate(
        cand_id="cand_C",
        dist=5,
        same_strand=False,
        is_operon=False,
        seq="M" * 50 + "DDDDDDDDDDDDDDDDDDDD",
    )

    clu = ClusterResult(
        clusters={"cand_A": "repA", "cand_B": "repB", "cand_C": "repC"},
        cluster_sizes={"repA": 2, "repB": 1, "repC": 1},
        cluster_shell_counts={"repA": 2, "repB": 1, "repC": 1},
        cluster_genome_counts={"repA": 2, "repB": 1, "repC": 1},
    )

    motifs = load_motif_patterns(MOTIF_FILE)
    weights = ScoringWeights(
        w_distance=0.25,
        w_operon=0.15,
        w_tp_consensus=0.15,
        w_cargo_class=0.30,
        w_tail_features=0.05,
        w_recurrence=0.10,
    )
    scored = score_candidates(
        [cand_c, cand_a, cand_b], clu, weights=weights, motifs=motifs
    )

    assert len(scored) == 3
    # Check strict ordering: A > B > C
    assert scored[0].candidate_id == "cand_A"
    assert scored[0].rank == 1
    assert scored[1].candidate_id == "cand_B"
    assert scored[1].rank == 2
    assert scored[2].candidate_id == "cand_C"
    assert scored[2].rank == 3

    assert scored[0].score_total > scored[1].score_total > scored[2].score_total

    # Verify every row retains its individual score components
    d0 = scored[0].to_dict()
    assert "score_distance" in d0
    assert "score_operon" in d0
    assert "score_tp_consensus" in d0
    assert "score_cargo_class" in d0
    assert "cargo_class" in d0
    assert "score_tail_features" in d0
    assert "score_recurrence" in d0
    assert math.isclose(
        d0["score_total"],
        sum(
            [
                weights.w_distance * scored[0].score_distance,
                weights.w_operon * scored[0].score_operon,
                weights.w_tp_consensus * scored[0].score_tp_consensus,
                weights.w_cargo_class * scored[0].score_cargo_class,
                weights.w_tail_features * scored[0].score_tail_features,
                weights.w_recurrence * scored[0].score_recurrence,
            ]
        ),
        rel_tol=1e-4,
    )


def test_score_thermotoga_slice_fixture():
    faa = FIXTURES_DIR / "thermotoga_slice.faa"
    gff = FIXTURES_DIR / "thermotoga_slice.gff3"
    gd = load_protein_and_gff3(faa, gff, genome_id="t_maritima")
    hits = detect_shells(gd, [PROFILE_PATH])
    cands = extract_neighborhood(gd, hits, window_size=5)

    clu = ClusterResult(
        clusters={c.candidate_id: c.candidate_id for c in cands},
        cluster_sizes={c.candidate_id: 1 for c in cands},
        cluster_shell_counts={c.candidate_id: 1 for c in cands},
        cluster_genome_counts={c.candidate_id: 1 for c in cands},
    )

    motifs = load_motif_patterns(MOTIF_FILE)
    scored = score_candidates(cands, clu, motifs=motifs)

    # Known cargo TM0786 homolog in slice is thermotoga_slice_0009
    top_cand = scored[0]
    assert top_cand.candidate_id == "thermotoga_slice_0009"
    assert top_cand.rank == 1
    assert top_cand.abs_gene_distance == 1
    assert top_cand.is_likely_operon is True
    assert top_cand.features.motif_matches["T_maritima_Flp_TP"] is True


def test_tail_score_hand_computed():
    # 10 A (hydrophobic) + 7 G (Gly/Ser) + 3 D in a 20-residue tail
    feats = compute_peptide_features("M" * 30 + "A" * 10 + "G" * 7 + "D" * 3)
    # 0.60 * 10/20 + 0.40 * min(1, (7/20) / 0.35) = 0.30 + 0.40
    assert math.isclose(compute_tail_features_score(feats), 0.70)


def test_tail_score_ignores_sequence_complexity():
    homopolymer = compute_peptide_features("M" * 30 + "D" * 20)
    assert compute_tail_features_score(homopolymer) == 0.0


def test_motif_takes_a_fifth_of_tail_score():
    motifs = {"m": "DDD$"}
    hit = compute_peptide_features("M" * 30 + "A" * 17 + "DDD", motifs=motifs)
    miss = compute_peptide_features("M" * 30 + "A" * 20, motifs=motifs)
    assert math.isclose(
        compute_tail_features_score(hit, has_motifs=True), 0.80 * 0.6 * 17 / 20 + 0.20
    )
    assert math.isclose(compute_tail_features_score(miss, has_motifs=True), 0.48)


def test_score_termini_uses_better_end_only_when_asked():
    feats = compute_peptide_features("A" * 20 + "K" * 40 + "D" * 20)
    assert score_termini(feats, termini="c") == (0.0, "c")
    assert score_termini(feats, termini="both") == (0.6, "n")


def test_rank_in_shell_is_per_shell():
    a1 = _build_mock_candidate("a1", 1, True, True, "M" * 50)
    a2 = _build_mock_candidate("a2", 3, True, False, "M" * 50)
    b1 = _build_mock_candidate("b1", 2, True, True, "M" * 50)
    b1.shell_gene_id = "shell_b"
    ids = ["a1", "a2", "b1"]
    clu = ClusterResult(
        {i: i for i in ids},
        {i: 1 for i in ids},
        {i: 1 for i in ids},
        {i: 1 for i in ids},
    )
    scored = {c.candidate_id: c for c in score_candidates([a1, a2, b1], clu)}
    assert (scored["a1"].rank, scored["a1"].rank_in_shell) == (1, 1)
    assert (scored["b1"].rank, scored["b1"].rank_in_shell) == (2, 1)
    assert (scored["a2"].rank, scored["a2"].rank_in_shell) == (3, 2)
