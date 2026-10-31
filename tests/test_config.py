"""Unit tests for configuration and scoring weights."""

from pathlib import Path

import pytest

import cargofinder
from cargofinder.config import (
    DEFAULT_SHELL_PROFILES,
    EXTENDED_SHELL_PROFILES,
    RunConfig,
    ScoringWeights,
    load_weights,
    select_profiles,
)


def test_default_scoring_weights_valid():
    weights = ScoringWeights()
    weights.validate()
    total = (
        weights.w_distance
        + weights.w_operon
        + weights.w_tail_features
        + weights.w_recurrence
    )
    assert 0.99 <= total <= 1.01


def test_invalid_negative_scoring_weights():
    weights = ScoringWeights(
        w_distance=-0.1, w_operon=0.5, w_tail_features=0.3, w_recurrence=0.3
    )
    with pytest.raises(ValueError, match="non-negative"):
        weights.validate()


def test_invalid_sum_scoring_weights():
    weights = ScoringWeights(
        w_distance=0.5, w_operon=0.5, w_tail_features=0.5, w_recurrence=0.5
    )
    with pytest.raises(ValueError, match="must sum to 1.0"):
        weights.validate()


def test_default_shell_profiles():
    assert len(DEFAULT_SHELL_PROFILES) >= 1
    p0 = DEFAULT_SHELL_PROFILES[0]
    assert p0.pfam_id == "PF04454"
    assert "Linocin_M18" in p0.family


def test_extended_shell_profiles():
    assert "PF19307" in EXTENDED_SHELL_PROFILES
    assert "PF05065" in EXTENDED_SHELL_PROFILES
    assert "PF08967" in EXTENDED_SHELL_PROFILES


def test_run_config_defaults():
    cfg = RunConfig()
    assert cfg.window_size == 5
    assert cfg.tail_length == 20
    assert cfg.max_operon_gap_bp == 150
    assert cfg.evalue_cutoff == 1e-5


def test_load_weights_inline_and_file(tmp_path):
    assert load_weights('{"w_distance": 0.35}') == ScoringWeights()
    shipped = Path(cargofinder.__file__).parent / "data" / "default_weights.json"
    assert load_weights(str(shipped)) == ScoringWeights()


def test_load_weights_rejects_unknown_key_and_bad_json():
    with pytest.raises(ValueError, match="Unknown scoring weight"):
        load_weights('{"distance": 0.9, "w_operon": 0.1}')
    with pytest.raises(ValueError, match="Failed to parse"):
        load_weights("{not json")


def test_select_profiles():
    assert select_profiles() == DEFAULT_SHELL_PROFILES
    assert select_profiles(["PF19307"])[-1].pfam_id == "PF19307"
    with pytest.raises(ValueError, match="Unknown extra profile"):
        select_profiles(["PF00000"])


@pytest.mark.parametrize(
    "kwargs",
    [
        {"window_size": 0},
        {"window_size": -2},
        {"tail_length": 0},
        {"max_operon_gap_bp": -1},
        {"evalue_cutoff": 0.0},
        {"termini": "n"},
    ],
)
def test_run_config_rejects_bad_values(kwargs):
    with pytest.raises(ValueError):
        RunConfig(**kwargs).validate()
