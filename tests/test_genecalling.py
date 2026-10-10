"""Unit tests for Pyrodigal gene calling."""

from pathlib import Path

from cargofinder.genecalling import call_genes_from_genome_fasta, call_genes_pyrodigal

FIXTURES_DIR = Path(__file__).parent / "data"


def test_call_genes_pyrodigal_synthetic():
    # Construct a valid synthetic coding frame (Start ATG, internal codons, stop TAA)
    orf = "ATG" + "GAA" * 50 + "TAA"
    contigs = {"syn_ctg": orf}
    gd = call_genes_pyrodigal(contigs, genome_id="synth", meta=True)
    assert len(gd) >= 1
    g = gd.genes[0]
    assert g.protein_seq.startswith("M")
    assert not g.protein_seq.endswith("*")


def test_call_genes_from_thermotoga_slice():
    slice_fna = FIXTURES_DIR / "thermotoga_slice.fna"
    assert slice_fna.exists()

    gd = call_genes_from_genome_fasta(slice_fna, genome_id="t_maritima_slice")
    assert len(gd) == 13

    # Check the encapsulin shell (TM0785 homolog in slice)
    g_shell = gd.get_gene("thermotoga_slice_0008")
    assert g_shell is not None
    assert g_shell.start == 3702
    assert g_shell.end == 4481
    assert g_shell.strand == -1
    assert g_shell.length_aa == 259

    # Check the ferritin-like cargo (TM0786 homolog in slice)
    g_cargo = gd.get_gene("thermotoga_slice_0009")
    assert g_cargo is not None
    assert g_cargo.start == 4496
    assert g_cargo.end == 4840
    assert g_cargo.strand == -1
    assert g_cargo.length_aa == 114
