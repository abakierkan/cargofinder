"""Unit tests for prophage context screening.

These cover the window logic and the flag rule without needing HMM profiles.
The end-to-end check that real phage markers are found, and that real encapsulin
loci are not flagged, lives in validation/run_validation.py; see the
"Prophage screening control" section of validation/report.md.
"""

from pathlib import Path

from cargofinder.config import PHAGE_MARKER_FLAG_THRESHOLD
from cargofinder.io import Gene, GenomeData
from cargofinder.prophage import (
    ProphageEvidence,
    assess_prophage_context,
    collect_prophage_windows,
    find_phage_marker_profiles,
    scan_genes_for_phage_markers,
)
from cargofinder.shells import ShellHit


def _gene(idx: int, contig: str = "ctg1") -> Gene:
    start = 100 + idx * 1000
    return Gene(f"g{idx}", contig, start, start + 600, 1, idx, "M" * 200)


def _shell_hit(gene: Gene) -> ShellHit:
    return ShellHit(
        gene.gene_id, "Linocin_M18", "PF04454", 1e-40, 150.0, 1e-40, 1, gene
    )


def test_evidence_counts_and_flag_rule():
    ev = ProphageEvidence(shell_gene_id="s1")
    assert ev.marker_count == 0
    assert ev.distinct_family_count == 0
    assert not ev.is_probable_prophage
    assert ev.marker_summary == ""

    # One family across two genes stays below the threshold of 2 families.
    ev.marker_families["g2"] = "Phage_integrase"
    ev.marker_families["g9"] = "Phage_integrase"
    assert ev.marker_count == 2
    assert ev.distinct_family_count == 1
    assert not ev.is_probable_prophage

    # A second distinct family trips the flag.
    ev.marker_families["g5"] = "TerL_ATPase"
    assert ev.distinct_family_count == 2
    assert ev.is_probable_prophage


def test_flag_threshold_matches_config():
    ev = ProphageEvidence(shell_gene_id="s1")
    for i in range(PHAGE_MARKER_FLAG_THRESHOLD - 1):
        ev.marker_families[f"g{i}"] = f"family_{i}"
    assert not ev.is_probable_prophage
    ev.marker_families["extra"] = f"family_{PHAGE_MARKER_FLAG_THRESHOLD}"
    assert ev.is_probable_prophage


def test_marker_summary_is_sorted_and_deterministic():
    ev = ProphageEvidence(shell_gene_id="s1")
    ev.marker_families["gZ"] = "TerL_ATPase"
    ev.marker_families["gA"] = "Phage_integrase"
    assert ev.marker_summary == "gA:Phage_integrase;gZ:TerL_ATPase"


def test_window_is_wider_than_cargo_window_and_excludes_shell():
    genes = [_gene(i) for i in range(41)]
    genome = GenomeData("g", genes)
    shell = _shell_hit(genes[20])

    windows = collect_prophage_windows([genome], [[shell]], prophage_window=20)
    picked = windows["g20"]

    assert len(picked) == 40, "20 genes either side, shell itself excluded"
    assert all(g.gene_id != "g20" for g in picked)
    assert {g.gene_id for g in picked} == {f"g{i}" for i in range(41)} - {"g20"}


def test_window_respects_contig_start_and_end():
    genes = [_gene(i) for i in range(10)]
    genome = GenomeData("g", genes)

    at_start = collect_prophage_windows(
        [genome], [[_shell_hit(genes[0])]], prophage_window=20
    )["g0"]
    assert len(at_start) == 9
    assert "g0" not in {g.gene_id for g in at_start}

    at_end = collect_prophage_windows(
        [genome], [[_shell_hit(genes[9])]], prophage_window=20
    )["g9"]
    assert len(at_end) == 9


def test_window_does_not_cross_contigs():
    """Contigs are linear and separate, matching extract_neighborhood."""
    a = [_gene(i, "ctgA") for i in range(5)]
    b = [
        Gene(f"h{i}", "ctgB", 100 + i * 1000, 700 + i * 1000, 1, i, "M" * 200)
        for i in range(5)
    ]
    genome = GenomeData("g", a + b)

    picked = collect_prophage_windows(
        [genome], [[_shell_hit(a[2])]], prophage_window=20
    )["g2"]
    assert all(g.contig_id == "ctgA" for g in picked)
    assert len(picked) == 4


def test_window_is_small_when_requested():
    genes = [_gene(i) for i in range(41)]
    genome = GenomeData("g", genes)
    picked = collect_prophage_windows(
        [genome], [[_shell_hit(genes[20])]], prophage_window=2
    )["g20"]
    assert {g.gene_id for g in picked} == {"g18", "g19", "g21", "g22"}


def test_no_marker_profiles_yields_empty_evidence_for_every_shell():
    """Without profiles the screen must report nothing, not fail or drop shells."""
    genes = [_gene(i) for i in range(5)]
    evidence = assess_prophage_context({"g0": genes, "g3": genes}, [])
    assert set(evidence) == {"g0", "g3"}
    assert all(not ev.is_probable_prophage for ev in evidence.values())
    assert all(ev.marker_count == 0 for ev in evidence.values())


def test_scan_with_no_profiles_or_no_genes_returns_empty():
    assert scan_genes_for_phage_markers([_gene(0)], []) == {}
    assert scan_genes_for_phage_markers([], ["nonexistent.hmm"]) == {}


def test_missing_marker_directory_is_not_an_error(tmp_path):
    """A profile dir predating this feature must degrade, not crash."""
    assert find_phage_marker_profiles(tmp_path) == []


def test_marker_directory_is_found_when_present(tmp_path):
    sub = tmp_path / "phage_markers"
    sub.mkdir()
    (sub / "PF00589.hmm").write_text("HMMER3/f placeholder\n", encoding="utf-8")
    found = find_phage_marker_profiles(tmp_path)
    assert [Path(p).name for p in found] == ["PF00589.hmm"]


def test_shell_with_no_window_genes_still_gets_an_entry():
    genes = [_gene(0)]
    genome = GenomeData("g", genes)
    windows = collect_prophage_windows(
        [genome], [[_shell_hit(genes[0])]], prophage_window=5
    )
    assert windows == {"g0": []}
    evidence = assess_prophage_context(windows, [])
    assert "g0" in evidence
