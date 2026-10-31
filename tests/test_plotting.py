"""Unit tests for neighborhood visualization."""

from pathlib import Path

from cargofinder.cluster import ClusterResult
from cargofinder.features import load_motif_patterns
from cargofinder.io import load_protein_and_gff3
from cargofinder.neighborhood import extract_neighborhood
from cargofinder.plotting import plot_neighborhood
from cargofinder.scoring import score_candidates
from cargofinder.shells import detect_shells

FIXTURES_DIR = Path(__file__).parent / "data"
PROFILE_PATH = FIXTURES_DIR / "PF04454_test.hmm"
MOTIF_FILE = Path(__file__).parent.parent / "src/cargofinder/data/example_motifs.json"


def test_plot_neighborhood_thermotoga_slice(tmp_path):
    faa = FIXTURES_DIR / "thermotoga_slice.faa"
    gff = FIXTURES_DIR / "thermotoga_slice.gff3"
    gd = load_protein_and_gff3(faa, gff, genome_id="t_maritima")
    hits = detect_shells(gd, [PROFILE_PATH])
    assert len(hits) == 1
    shell = hits[0]

    cands = extract_neighborhood(gd, hits, window_size=5)
    clu = ClusterResult(
        clusters={c.candidate_id: c.candidate_id for c in cands},
        cluster_sizes={c.candidate_id: 1 for c in cands},
        cluster_shell_counts={c.candidate_id: 1 for c in cands},
        cluster_genome_counts={c.candidate_id: 1 for c in cands},
    )
    motifs = load_motif_patterns(MOTIF_FILE)
    scored = score_candidates(cands, clu, motifs=motifs)

    out_png = tmp_path / "neighborhood_TM0785.png"
    saved = plot_neighborhood(
        shell=shell,
        genome=gd,
        candidates=scored,
        out_path=out_png,
        flank_genes=5,
    )

    assert saved.exists()
    assert saved.stat().st_size > 1000  # valid image file generated
    # Check PNG magic bytes
    with saved.open("rb") as fh:
        magic = fh.read(8)
        assert magic == b"\x89PNG\r\n\x1a\n"
