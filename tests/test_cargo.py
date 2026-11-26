"""Unit tests for known cargo class recognition.

The end-to-end check that real cargo gets the right class lives in
validation/run_validation.py; see the "Cargo class recognition" section of
validation/report.md.
"""

from pathlib import Path

import pytest

from cargofinder.cargo import (
    _BY_ACCESSION,
    CargoClassHit,
    _definition_for,
    _outranks,
    classify_cargo_candidates,
    find_cargo_class_profiles,
)
from cargofinder.config import (
    CARGO_CLASS_PROFILES,
    CARGO_CLASS_SCORE_GENERAL,
    CARGO_CLASS_SCORE_SPECIFIC,
)
from cargofinder.io import Gene


def _hit(pfam: str, specific: bool, evalue: float = 1e-20) -> CargoClassHit:
    return CargoClassHit(
        gene_id="g1",
        cargo_class="cls",
        pfam_id=pfam,
        profile_name=pfam,
        evalue=evalue,
        encapsulin_specific=specific,
    )


def test_profile_accessions_are_well_formed_and_unique():
    accs = [d.pfam_id for d in CARGO_CLASS_PROFILES]
    assert len(accs) == len(set(accs)), "duplicate Pfam accession in cargo classes"
    for a in accs:
        assert a.startswith("PF") and a[2:].isdigit() and len(a) == 7, a


def test_exactly_the_two_encapsulin_specific_families_are_marked():
    """PF22277 and PF24309 are the only Pfam families defined as encapsulin cargo."""
    specific = {d.pfam_id for d in CARGO_CLASS_PROFILES if d.encapsulin_specific}
    assert specific == {"PF22277", "PF24309"}


def test_accessory_components_are_not_listed_as_cargo():
    """Andreas & Giessen state accessory operon components are not cargo proteins.

    PF00581 is Rhodanese and PF00583 is the acetyltransferase family; both are
    enriched Family 2A accessory components, and neither belongs here.
    """
    accs = {d.pfam_id for d in CARGO_CLASS_PROFILES}
    assert "PF00581" not in accs
    assert "PF00583" not in accs


def test_score_tiers():
    assert _hit("PF22277", True).score == CARGO_CLASS_SCORE_SPECIFIC
    assert _hit("PF00210", False).score == CARGO_CLASS_SCORE_GENERAL
    assert _hit("PF22277", True).score > _hit("PF00210", False).score


def test_encapsulin_specific_outranks_general_even_with_worse_evalue():
    """A ferritin domain hit must not displace an encapsulated-ferritin call."""
    specific_weak = _hit("PF22277", True, evalue=1e-5)
    general_strong = _hit("PF00210", False, evalue=1e-80)
    assert _outranks(specific_weak, general_strong)
    assert not _outranks(general_strong, specific_weak)


def test_within_a_tier_the_better_evalue_wins():
    assert _outranks(_hit("PF00210", False, 1e-40), _hit("PF02915", False, 1e-10))
    assert not _outranks(_hit("PF00210", False, 1e-10), _hit("PF02915", False, 1e-40))
    assert not _outranks(_hit("PF22277", True, 1e-10), _hit("PF24309", True, 1e-40))


def test_definition_lookup_strips_the_pfam_version():
    """HMM headers carry versioned accessions such as PF00266.26."""
    assert _definition_for("Aminotran_5", "PF00266.26") is _BY_ACCESSION["PF00266"]
    assert _definition_for("Aminotran_5", "PF00266") is _BY_ACCESSION["PF00266"]


def test_unknown_profile_resolves_to_nothing():
    assert _definition_for("Something_else", "PF99999.1") is None


def test_missing_cargo_directory_is_not_an_error(tmp_path):
    """A profile dir predating this feature must degrade, not crash."""
    assert find_cargo_class_profiles(tmp_path) == []


def test_cargo_directory_is_found_when_present(tmp_path):
    sub = tmp_path / "cargo_classes"
    sub.mkdir()
    (sub / "PF22277.hmm").write_text("HMMER3/f placeholder\n", encoding="utf-8")
    assert [Path(p).name for p in find_cargo_class_profiles(tmp_path)] == [
        "PF22277.hmm"
    ]


def test_classify_with_no_profiles_or_no_genes_returns_empty():
    gene = Gene("g1", "c1", 1, 300, 1, 0, "M" * 100)
    assert classify_cargo_candidates([gene], []) == {}
    assert classify_cargo_candidates([], ["nonexistent.hmm"]) == {}


@pytest.mark.parametrize("field", ["cargo_class", "pfam_id"])
def test_hit_carries_identifying_fields(field):
    assert getattr(_hit("PF22277", True), field)
