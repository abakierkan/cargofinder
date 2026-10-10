"""Unit tests for terminal peptide feature computation on hand-computed examples."""

import math
from pathlib import Path

from cargofinder.features import (
    compute_gly_ser_fraction,
    compute_hydrophobic_fraction,
    compute_hydrophobicity_kd,
    compute_low_complexity_score,
    compute_net_charge,
    compute_peptide_features,
    compute_shannon_entropy,
    extract_terminal_peptides,
    load_motif_patterns,
    match_motifs,
)

MOTIF_FILE = Path(__file__).parent.parent / "src/cargofinder/data/example_motifs.json"


def test_extract_terminal_peptides_short_and_long():
    # Short peptide (< 20 aa)
    n_term, c_term = extract_terminal_peptides("MVKTY", tail_length=20)
    assert n_term == "MVKTY"
    assert c_term == "MVKTY"

    # Exact length 20
    exact = "A" * 20
    n_term, c_term = extract_terminal_peptides(exact, tail_length=20)
    assert n_term == exact
    assert c_term == exact

    # Long sequence: 10 'A' + 10 'C' + 10 'D' (length 30)
    seq = "A" * 10 + "C" * 10 + "D" * 10
    n_term, c_term = extract_terminal_peptides(seq, tail_length=20)
    assert n_term == "A" * 10 + "C" * 10
    assert c_term == "C" * 10 + "D" * 10


def test_hand_computed_hydrophobicity():
    # IVL: I (4.5) + V (4.2) + L (3.8) = 12.5 / 3 = 4.16666...
    kd_ivl = compute_hydrophobicity_kd("IVL")
    assert math.isclose(kd_ivl, 12.5 / 3, rel_tol=1e-5)
    assert compute_hydrophobic_fraction("IVL") == 1.0

    # KRED: K (-3.9) + R (-4.5) + E (-3.5) + D (-3.5) = -15.4 / 4 = -3.85
    kd_kred = compute_hydrophobicity_kd("KRED")
    assert math.isclose(kd_kred, -3.85, rel_tol=1e-5)
    assert compute_hydrophobic_fraction("KRED") == 0.0


def test_hand_computed_gly_ser_fraction():
    assert compute_gly_ser_fraction("GSGG") == 1.0
    assert compute_gly_ser_fraction("ACDE") == 0.0
    assert compute_gly_ser_fraction("GSAC") == 0.5


def test_hand_computed_shannon_entropy_and_low_complexity():
    # Homopolymer of length 10: Shannon entropy is 0.0, low complexity score is 1.0
    assert compute_shannon_entropy("AAAAAAAAAA") == 0.0
    assert compute_low_complexity_score("AAAAAAAAAA") == 1.0

    # Equimolar 2-letter repeat: "AGAGAGAG" (4 A, 4 G)
    # H = - (0.5 * log2(0.5) + 0.5 * log2(0.5)) = 1.0 bit
    # max_H = log2(8) = 3.0 bits
    # low_complexity = 1.0 - (1.0 / 3.0) = 2/3
    h_ag = compute_shannon_entropy("AGAGAGAG")
    assert math.isclose(h_ag, 1.0, rel_tol=1e-5)
    lc_ag = compute_low_complexity_score("AGAGAGAG")
    assert math.isclose(lc_ag, 2.0 / 3.0, rel_tol=1e-5)


def test_hand_computed_net_charge():
    assert compute_net_charge("KRDE") == 0.0
    assert compute_net_charge("KKRR") == 4.0
    assert compute_net_charge("DDDE") == -4.0


def test_motif_matching_and_example_file():
    patterns = load_motif_patterns(MOTIF_FILE)
    assert "T_maritima_Flp_TP" in patterns

    # Positive match for T. maritima cargo tail
    tm_tail = "EEETSGGSENTGGDLGIRKL"
    matches = match_motifs(tm_tail, patterns)
    assert matches["T_maritima_Flp_TP"] is True
    assert matches["Generic_hydrophobic_c_terminus"] is True

    # Negative match
    neg_tail = "ACDEFGHIKLMNPQRSTVWY"
    matches_neg = match_motifs(neg_tail, patterns)
    assert matches_neg["T_maritima_Flp_TP"] is False


def test_compute_peptide_features_tm0786():
    tm0786_tail = (
        "MADQYHEPVSELTGKDRDFVRALNSLKEEIEAVAWYHQRVVTTKDETVRKILEHNRDEEMEHAAMLLEWLRRNMPGW"
        "DEALRTYLFTDKPITEIEEETSGGSENTGGDLGIRKL"
    )
    motifs = load_motif_patterns(MOTIF_FILE)
    feats = compute_peptide_features(tm0786_tail, tail_length=20, motifs=motifs)

    assert feats.c_term_seq == "EEETSGGSENTGGDLGIRKL"
    assert feats.motif_matches["T_maritima_Flp_TP"] is True
    # In EEETSGGSENTGGDLGIRKL, Gly count is 5, Ser count is 2 -> 7/20 = 0.35
    assert feats.c_gly_ser_fraction == 0.35
    d = feats.to_dict()
    assert "c_kd_hydrophobicity" in d
    assert d["motif_T_maritima_Flp_TP"] is True
