"""Unit tests for encapsulin shell detection using pyhmmer."""

from pathlib import Path

import pytest

from cargofinder.io import Gene, GenomeData, load_protein_and_gff3
from cargofinder.shells import detect_shells, load_hmm_profiles

FIXTURES_DIR = Path(__file__).parent / "data"
PROFILE_PATH = FIXTURES_DIR / "PF04454_test.hmm"


def test_load_hmm_profiles_missing():
    with pytest.raises(FileNotFoundError):
        load_hmm_profiles([FIXTURES_DIR / "nonexistent.hmm"])


def test_detect_shells_thermotoga_fixture():
    faa = FIXTURES_DIR / "thermotoga_slice.faa"
    gff = FIXTURES_DIR / "thermotoga_slice.gff3"
    gd = load_protein_and_gff3(faa, gff, genome_id="t_maritima")

    hits = detect_shells(gd, [PROFILE_PATH], evalue_cutoff=1e-5)

    assert len(hits) == 1
    hit = hits[0]
    assert hit.gene_id == "thermotoga_slice_0008"
    assert "Linocin_M18" in hit.profile_name
    assert hit.evalue < 1e-50
    assert hit.bitscore > 200.0
    assert hit.start == 3702
    assert hit.end == 4481
    assert hit.strand == -1


def test_detect_shells_empty_genome():
    gd = GenomeData(genome_id="empty", genes=[])
    hits = detect_shells(gd, [PROFILE_PATH])
    assert hits == []


def test_detect_shells_no_matching_genes():
    # Synthetic genes that do not match the encapsulin profile
    genes = [
        Gene("g1", "ctg1", 100, 700, 1, 0, "M" * 200),
        Gene("g2", "ctg1", 800, 1400, 1, 1, "G" * 200),
    ]
    gd = GenomeData(genome_id="non_shell", genes=genes)
    hits = detect_shells(gd, [PROFILE_PATH])
    assert hits == []


def test_detect_shells_min_length_filter():
    faa = FIXTURES_DIR / "thermotoga_slice.faa"
    gff = FIXTURES_DIR / "thermotoga_slice.gff3"
    gd = load_protein_and_gff3(faa, gff, genome_id="t_maritima")

    # The shell is 259 aa. Filtering with min_length_aa=300 should exclude it.
    hits = detect_shells(gd, [PROFILE_PATH], min_length_aa=300)
    assert len(hits) == 0
