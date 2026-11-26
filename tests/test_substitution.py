"""Tests for BLOSUM62 and the published-TP-consensus similarity feature."""

from __future__ import annotations

import math

import pytest

from cargofinder.features import (
    compute_peptide_features,
    compute_tp_consensus_similarity,
)
from cargofinder.substitution import (
    _BLOSUM62_ORDER,
    BLOSUM62,
    FAMILY1_TP_CONSENSUS,
    best_ungapped_similarity,
    score_pair,
    self_score,
)


def test_blosum62_is_complete_and_symmetric():
    assert len(BLOSUM62) == 20 * 20
    for a in _BLOSUM62_ORDER:
        for b in _BLOSUM62_ORDER:
            assert BLOSUM62[(a, b)] == BLOSUM62[(b, a)], f"asymmetric at {a}{b}"


@pytest.mark.parametrize(
    ("pair", "expected"),
    [
        (("W", "W"), 11),
        (("C", "C"), 9),
        (("H", "H"), 8),
        (("A", "A"), 4),
        (("I", "L"), 2),
        (("K", "R"), 2),
        (("D", "E"), 2),
        (("F", "Y"), 3),
        (("G", "W"), -2),
        (("P", "F"), -4),
        (("D", "W"), -4),
    ],
)
def test_blosum62_known_values(pair, expected):
    """Spot-check against the published NCBI BLOSUM62 matrix."""
    assert score_pair(*pair) == expected


def test_unknown_residues_score_zero():
    assert score_pair("X", "A") == 0
    assert score_pair("A", "*") == 0
    assert self_score("XXX") == 0


def test_self_score_is_sum_of_diagonal():
    # A=4, C=9, W=11 -> 24
    assert self_score("ACW") == 24
    assert self_score(FAMILY1_TP_CONSENSUS) == 62


def test_exact_consensus_match_scores_one():
    assert best_ungapped_similarity(FAMILY1_TP_CONSENSUS) == 1.0


def test_consensus_embedded_in_longer_peptide_is_found():
    """The reference is slid along the peptide, so position must not matter."""
    assert best_ungapped_similarity("AAAA" + FAMILY1_TP_CONSENSUS + "AAAA") == 1.0
    assert best_ungapped_similarity(FAMILY1_TP_CONSENSUS + "PPPPPP") == 1.0


def test_similarity_is_bounded_and_clamped():
    """A peptide that aligns only negatively must clamp to 0.0, never go below."""
    # Trp/Pro run against a Gly/Ser-rich consensus scores negative everywhere.
    assert best_ungapped_similarity("WWWWWWWWWWWW") == 0.0
    for pep in ["", "A", "MKV", "P" * 40, FAMILY1_TP_CONSENSUS]:
        assert 0.0 <= best_ungapped_similarity(pep) <= 1.0


def test_empty_inputs_return_zero():
    assert best_ungapped_similarity("") == 0.0
    assert best_ungapped_similarity("GSLGIGSL", reference="") == 0.0


def test_altenburg_consensus_tp_ranks_above_unrelated_peptide():
    """The T1 consensus TP from Altenburg et al. (2021) must outrank junk.

    DGSLGIGSLKG is the T1 consensus TP characterised in Altenburg et al. (2021)
    Sci Rep 11:4951; it shares the GSLGIGSL core with the Giessen & Silver (2017)
    consensus used as the reference here.
    """
    t1_consensus = best_ungapped_similarity("DGSLGIGSLKG")
    unrelated = best_ungapped_similarity("MKKFLVWWPHQ")
    assert t1_consensus > 0.7
    assert t1_consensus > unrelated


def test_glycine_serine_linker_scores_background_not_zero():
    """Documents a known confounder rather than asserting the ideal behaviour.

    The reference consensus is itself Gly/Ser-rich, so a GGS linker picks up
    appreciable background similarity even though Altenburg et al. (2021) showed
    the GGSGGSGGSGG control peptide does not bind the shell pocket. If this
    assumption ever changes, this test should be revisited along with the
    docstring in features.compute_tp_consensus_similarity.
    """
    ggs = best_ungapped_similarity("GGSGGSGGSGG")
    assert 0.2 < ggs < 0.45


def test_case_and_whitespace_are_normalised():
    assert best_ungapped_similarity(f"  {FAMILY1_TP_CONSENSUS.lower()}  ") == 1.0


def test_feature_wrapper_matches_underlying_function():
    assert compute_tp_consensus_similarity(FAMILY1_TP_CONSENSUS) == 1.0
    assert math.isclose(
        compute_tp_consensus_similarity("DGSLGIGSLKG"),
        best_ungapped_similarity("DGSLGIGSLKG"),
    )


def test_features_report_both_termini_independently():
    """The consensus is at the C-terminus only, so the N-terminal value must be low."""
    protein = "M" + "A" * 40 + FAMILY1_TP_CONSENSUS
    feats = compute_peptide_features(protein, tail_length=20)
    assert feats.c_tp_consensus_similarity == 1.0
    assert feats.n_tp_consensus_similarity < 0.2
    d = feats.to_dict()
    assert d["c_tp_consensus_similarity"] == 1.0
    assert "n_tp_consensus_similarity" in d


def test_tail_length_can_truncate_the_consensus():
    """A slice shorter than the reference cannot reach 1.0; it degrades, not errors."""
    protein = "A" * 40 + FAMILY1_TP_CONSENSUS
    short = compute_peptide_features(protein, tail_length=5)
    full = compute_peptide_features(protein, tail_length=20)
    assert short.c_tp_consensus_similarity < full.c_tp_consensus_similarity
