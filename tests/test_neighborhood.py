"""Unit tests for genomic neighborhood window extraction and operon inference."""

from pathlib import Path

from cargofinder.io import Gene, GenomeData, load_protein_and_gff3
from cargofinder.neighborhood import (
    _check_operon_continuity,
    _compute_intergenic_distance,
    extract_neighborhood,
)
from cargofinder.shells import ShellHit, detect_shells

FIXTURES_DIR = Path(__file__).parent / "data"
PROFILE_PATH = FIXTURES_DIR / "PF04454_test.hmm"


def test_intergenic_distance_gap_and_overlap():
    g1 = Gene("g1", "c1", 100, 200, 1, 0, "M" * 33)
    g2 = Gene("g2", "c1", 250, 400, 1, 1, "M" * 50)
    # Gap: 250 - 200 - 1 = 49 bp
    assert _compute_intergenic_distance(g1, g2) == 49

    # Overlap: g2 starts at 195 (overlaps g1 by 6 bp)
    g2_ovl = Gene("g2", "c1", 195, 400, 1, 1, "M" * 68)
    # Overlap is -(200 - 195 + 1) = -6
    assert _compute_intergenic_distance(g1, g2_ovl) == -6


def test_contig_edge_boundary_start():
    # Shell is at index 0 of the contig
    shell_gene = Gene("shell", "ctg1", 100, 500, 1, 0, "M" * 133)
    g1 = Gene("g1", "ctg1", 520, 800, 1, 1, "M" * 93)
    g2 = Gene("g2", "ctg1", 850, 1200, 1, 2, "M" * 116)

    gd = GenomeData("genome_edge", [shell_gene, g1, g2])
    shell_hit = ShellHit(
        "shell", "Linocin", "PF04454", 1e-40, 150.0, 1e-40, 1, shell_gene
    )

    cands = extract_neighborhood(gd, [shell_hit], window_size=5)
    assert len(cands) == 2
    assert [c.candidate_id for c in cands] == ["g1", "g2"]
    assert cands[0].gene_distance == 1
    assert cands[1].gene_distance == 2


def test_contig_edge_boundary_end():
    # Shell is at the last index of the contig
    g1 = Gene("g1", "ctg1", 100, 400, 1, 0, "M" * 100)
    g2 = Gene("g2", "ctg1", 450, 800, 1, 1, "M" * 116)
    shell_gene = Gene("shell", "ctg1", 850, 1200, 1, 2, "M" * 116)

    gd = GenomeData("genome_edge_end", [g1, g2, shell_gene])
    shell_hit = ShellHit(
        "shell", "Linocin", "PF04454", 1e-40, 150.0, 1e-40, 1, shell_gene
    )

    cands = extract_neighborhood(gd, [shell_hit], window_size=5)
    assert len(cands) == 2
    assert [c.candidate_id for c in cands] == ["g1", "g2"]
    assert cands[0].gene_distance == -2
    assert cands[1].gene_distance == -1


def test_single_gene_contig():
    shell_gene = Gene("shell", "ctg1", 100, 500, 1, 0, "M" * 133)
    gd = GenomeData("genome_single", [shell_gene])
    shell_hit = ShellHit(
        "shell", "Linocin", "PF04454", 1e-40, 150.0, 1e-40, 1, shell_gene
    )

    cands = extract_neighborhood(gd, [shell_hit], window_size=5)
    assert cands == []


def test_operon_continuity_broken_by_strand_or_gap():
    # 3 genes on same contig
    g1 = Gene("g1", "ctg1", 100, 300, 1, 0, "M" * 66)
    # g2 opposite strand
    g2 = Gene("g2", "ctg1", 320, 500, -1, 1, "M" * 60)
    g3 = Gene("g3", "ctg1", 520, 700, 1, 2, "M" * 60)

    # Continuity between 0 and 2 broken by opposite strand at 1
    assert not _check_operon_continuity(
        [g1, g2, g3], 0, 2, required_strand=1, max_gap_bp=150
    )

    # Large gap > 150 bp
    g2_same_strand = Gene("g2", "ctg1", 320, 500, 1, 1, "M" * 60)
    g3_distant = Gene("g3", "ctg1", 1000, 1200, 1, 2, "M" * 66)
    assert not _check_operon_continuity(
        [g1, g2_same_strand, g3_distant], 0, 2, required_strand=1, max_gap_bp=150
    )


def test_extract_neighborhood_thermotoga_fixture():
    faa = FIXTURES_DIR / "thermotoga_slice.faa"
    gff = FIXTURES_DIR / "thermotoga_slice.gff3"
    gd = load_protein_and_gff3(faa, gff, genome_id="t_maritima")
    hits = detect_shells(gd, [PROFILE_PATH])
    assert len(hits) == 1

    cands = extract_neighborhood(gd, hits, window_size=5, max_operon_gap_bp=150)
    # Window of 5 around gene 8 (index 7, with 13 total genes)
    # Contig indices [7-5=2 to 7+5=12], total 11 genes minus shell = 10 candidates
    assert len(cands) == 10

    # Specifically verify the adjacent known cargo: thermotoga_slice_0009
    cargo = next(c for c in cands if c.candidate_id == "thermotoga_slice_0009")
    assert cargo.abs_gene_distance == 1
    assert cargo.gene_distance == 1
    # On negative strand, neighbor index 8 (gene 9) is upstream in transcription
    assert cargo.strand_relative_distance == -1
    assert cargo.same_strand is True
    assert cargo.bp_distance == 14
    assert cargo.is_likely_operon is True

    # Check a divergent neighbor: thermotoga_slice_0007 (strand +1)
    cand7 = next(c for c in cands if c.candidate_id == "thermotoga_slice_0007")
    assert cand7.abs_gene_distance == 1
    assert cand7.same_strand is False
    assert cand7.is_likely_operon is False
