"""Unit and integration tests for MMseqs2 candidate clustering."""

from unittest.mock import patch

import pytest

from cargofinder.cluster import (
    ClusterResult,
    cluster_candidates,
    find_mmseqs_binary,
)
from cargofinder.io import Gene
from cargofinder.neighborhood import NeighborCandidate
from cargofinder.shells import ShellHit


def _make_dummy_candidate(
    cand_id: str,
    shell_id: str,
    seq: str,
    contig: str = "ctg1",
) -> NeighborCandidate:
    cand_gene = Gene(cand_id, contig, 100, 200, 1, 0, seq)
    shell_gene = Gene(shell_id, contig, 300, 500, 1, 1, "M" * 200)
    shell_hit = ShellHit(
        shell_id, "Linocin", "PF04454", 1e-40, 100.0, 1e-40, 1, shell_gene
    )
    return NeighborCandidate(
        candidate_id=cand_id,
        shell_gene_id=shell_id,
        contig_id=contig,
        gene_distance=1,
        abs_gene_distance=1,
        strand_relative_distance=1,
        bp_distance=50,
        same_strand=True,
        is_likely_operon=True,
        candidate_gene=cand_gene,
        shell_hit=shell_hit,
    )


def test_cluster_empty_and_single():
    # Empty
    res_empty = cluster_candidates([])
    assert len(res_empty.clusters) == 0

    # Single candidate does not require running external MMseqs2
    c1 = _make_dummy_candidate("c1", "s1", "MVKTY" * 10)
    res_single = cluster_candidates([c1])
    assert res_single.get_cluster_id("c1") == "c1"
    assert res_single.get_cluster_size("c1") == 1
    assert res_single.get_shell_recurrence("c1") == 1


def test_find_mmseqs_binary_missing():
    with (
        patch("shutil.which", return_value=None),
        patch("pathlib.Path.is_file", return_value=False),
    ):
        with pytest.raises(RuntimeError, match="MMseqs2 executable .* was not found"):
            find_mmseqs_binary()


def test_find_mmseqs_binary_invalid_custom_path():
    with pytest.raises(
        FileNotFoundError, match="Specified MMseqs2 executable does not exist"
    ):
        find_mmseqs_binary(custom_path="/nonexistent/path/mmseqs")


def test_cluster_result_methods():
    res = ClusterResult(
        clusters={"c1": "repA", "c2": "repA", "c3": "repB"},
        cluster_sizes={"repA": 2, "repB": 1},
        cluster_shell_counts={"repA": 2, "repB": 1},
        cluster_genome_counts={"repA": 2, "repB": 1},
    )
    assert res.get_cluster_id("c1") == "repA"
    assert res.get_cluster_size("c1") == 2
    assert res.get_shell_recurrence("c1") == 2
    assert res.get_shell_recurrence("c3") == 1


@pytest.mark.external
def test_cluster_candidates_with_real_mmseqs():
    # Check if binary is present in environment before running
    try:
        binary = find_mmseqs_binary()
    except RuntimeError:
        pytest.skip("MMseqs2 binary not found in environment")

    # Construct two identical candidates next to two DIFFERENT shells
    seq_homolog = "MVKTY" * 20
    seq_distinct = "AILVF" * 20

    c1 = _make_dummy_candidate("c1", "shell_1", seq_homolog, contig="genome1")
    c2 = _make_dummy_candidate("c2", "shell_2", seq_homolog, contig="genome2")
    c3 = _make_dummy_candidate("c3", "shell_1", seq_distinct, contig="genome1")

    res = cluster_candidates([c1, c2, c3], mmseqs_bin=binary, min_seq_id=0.8)
    assert len(res.clusters) == 3

    # c1 and c2 should share the same cluster ID
    assert res.get_cluster_id("c1") == res.get_cluster_id("c2")
    # c3 should be in a distinct cluster
    assert res.get_cluster_id("c3") != res.get_cluster_id("c1")

    # Homolog cluster should have size 2 and recurrence across 2 distinct shells!
    assert res.get_cluster_size("c1") == 2
    assert res.get_shell_recurrence("c1") == 2
    assert res.get_shell_recurrence("c3") == 1
