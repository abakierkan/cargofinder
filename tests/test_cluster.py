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


def _cand(gene_id, shell_id, genome="g1"):
    from cargofinder.io import Gene
    from cargofinder.neighborhood import NeighborCandidate
    from cargofinder.shells import ShellHit

    gene = Gene(gene_id, "c1", 100, 700, 1, 0, "M" * 200)
    shell_gene = Gene(shell_id, "c1", 900, 1500, 1, 1, "M" * 200)
    hit = ShellHit(shell_id, "Linocin_M18", "PF04454", 1e-40, 1.0, 1e-40, 1, shell_gene)
    return NeighborCandidate(
        candidate_id=f"{shell_id}|{gene_id}",
        shell_gene_id=shell_id,
        contig_id="c1",
        gene_distance=1,
        abs_gene_distance=1,
        strand_relative_distance=1,
        bp_distance=50,
        same_strand=True,
        is_likely_operon=True,
        candidate_gene=gene,
        shell_hit=hit,
        genome_id=genome,
    )


def test_tandem_shells_do_not_manufacture_recurrence():
    """One gene next to two adjacent shells is one gene, not two occurrences.

    Tandem shells are real: M. xanthus Family 2B operons encode two adjacent
    shell proteins, so their windows overlap and share neighbours.
    """
    from cargofinder.cluster import _count_independent_genes

    same_gene_two_shells = [_cand("geneA", "shell1"), _cand("geneA", "shell2")]
    assert _count_independent_genes(same_gene_two_shells) == 1


def test_distinct_genes_still_count_as_recurrence():
    from cargofinder.cluster import _count_independent_genes

    members = [_cand("geneA", "shell1"), _cand("geneB", "shell2")]
    assert _count_independent_genes(members) == 2


def test_same_gene_id_in_different_genomes_counts_twice():
    """Gene ids are only unique within a genome."""
    from cargofinder.cluster import _count_independent_genes

    members = [_cand("gene_0001", "s1", "genomeA"), _cand("gene_0001", "s2", "genomeB")]
    assert _count_independent_genes(members) == 2
