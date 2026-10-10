"""Unit tests for the genomic I/O layer (FASTA and GFF3 handling)."""

from pathlib import Path

import pytest

from cargofinder.io import (
    Gene,
    GenomeData,
    load_protein_and_gff3,
    parse_gff3,
    read_fasta_sequences,
    write_fasta_sequences,
)

FIXTURES_DIR = Path(__file__).parent / "data"


def test_gene_model_validation():
    # Valid gene
    g = Gene(
        gene_id="g1",
        contig_id="chr1",
        start=100,
        end=250,
        strand=1,
        gene_index=0,
        protein_seq="MVKTY*",
    )
    assert g.length_nt == 151
    assert g.length_aa == 5  # trailing * stripped
    assert g.protein_seq == "MVKTY"

    # Start > end raises ValueError
    with pytest.raises(ValueError, match="cannot exceed end"):
        Gene(
            gene_id="g2",
            contig_id="chr1",
            start=300,
            end=200,
            strand=1,
            gene_index=1,
            protein_seq="M",
        )

    # Invalid strand
    with pytest.raises(ValueError, match="strand must be"):
        Gene(
            gene_id="g3",
            contig_id="chr1",
            start=100,
            end=200,
            strand=0,
            gene_index=1,
            protein_seq="M",
        )


def test_genome_data_sorting():
    genes = [
        Gene("g2", "ctgA", 500, 700, 1, 0, "MK"),
        Gene("g1", "ctgA", 100, 300, -1, 0, "MA"),
        Gene("g3", "ctgB", 50, 150, 1, 0, "MV"),
    ]
    gd = GenomeData(genome_id="test_gen", genes=genes)
    assert len(gd) == 3
    ctg_a = gd.genes_on_contig("ctgA")
    assert [g.gene_id for g in ctg_a] == ["g1", "g2"]
    assert ctg_a[0].gene_index == 0
    assert ctg_a[1].gene_index == 1
    assert gd.get_gene("g3").contig_id == "ctgB"


def test_read_fasta_sequences(tmp_path):
    fasta_file = tmp_path / "test.fa"
    fasta_file.write_text(
        ">seq1 First sequence description\natgc\nGGTT\n\n>seq2\nCCAA\n"
    )
    records = read_fasta_sequences(fasta_file)
    assert records == {"seq1": "ATGCGGTT", "seq2": "CCAA"}


def test_read_fasta_duplicate_id_raises(tmp_path):
    fasta_file = tmp_path / "dup.fa"
    fasta_file.write_text(">seq1\nAAA\n>seq1\nCCC\n")
    with pytest.raises(ValueError, match="Duplicate FASTA ID"):
        read_fasta_sequences(fasta_file)


def test_write_and_read_fasta(tmp_path):
    out_file = tmp_path / "written.fa"
    data = {"seqA": "M" * 70, "seqB": "K" * 15}
    write_fasta_sequences(data, out_file, line_width=60)
    read_back = read_fasta_sequences(out_file)
    assert read_back == data


def test_parse_gff3_errors(tmp_path):
    # Malformed column count
    bad_gff = tmp_path / "bad.gff3"
    bad_gff.write_text("ctg1\tsource\tCDS\t100\t200\n")
    with pytest.raises(ValueError, match="expected 9 tab-separated"):
        parse_gff3(bad_gff)

    # Non-integer coordinates
    bad_coord = tmp_path / "bad_coord.gff3"
    bad_coord.write_text("ctg1\tS\tCDS\tX\t200\t.\t+\t0\tID=g1\n")
    with pytest.raises(ValueError, match="Invalid integer coordinates"):
        parse_gff3(bad_coord)


def test_load_protein_and_gff3_fixture():
    faa = FIXTURES_DIR / "thermotoga_slice.faa"
    gff = FIXTURES_DIR / "thermotoga_slice.gff3"
    assert faa.exists() and gff.exists()

    gd = load_protein_and_gff3(faa, gff, genome_id="thermotoga_fixture")
    assert len(gd) == 13
    g8 = gd.get_gene("thermotoga_slice_0008")
    assert g8 is not None
    assert g8.start == 3702
    assert g8.end == 4481
    assert g8.strand == -1
    assert g8.length_aa == 259
