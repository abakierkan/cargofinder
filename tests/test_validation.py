"""Tests for the validation harness and the integrity of its truth table."""

import importlib.util
import sys
from pathlib import Path

import pandas as pd

VALIDATION_DIR = Path(__file__).parent.parent / "validation"
_spec = importlib.util.spec_from_file_location(
    "run_validation", VALIDATION_DIR / "run_validation.py"
)
run_validation = importlib.util.module_from_spec(_spec)
sys.modules["run_validation"] = run_validation
_spec.loader.exec_module(run_validation)


def test_known_systems_rows_are_complete():
    known = pd.read_csv(VALIDATION_DIR / "known_systems.tsv", sep="\t", dtype=str)
    assert not known.isna().any().any()
    assert set(known["evidence"]) <= {"experimental", "homology_inferred"}
    assert set(known["cargo_location_type"]) <= {
        "adjacent_operon",
        "distal_genomic",
        "natural_fusion",
        "family_2a_operon",
    }


def test_parse_annotation_uses_parent_old_locus_tag(tmp_path):
    gff = tmp_path / "a.gff3"
    gff.write_text(
        "c\tR\tgene\t10\t99\t.\t-\t.\tID=gene-1;locus_tag=X_RS1;old_locus_tag=X0001\n"
        "c\tR\tCDS\t10\t99\t.\t-\t0\tID=cds-1;Parent=gene-1;locus_tag=X_RS1;product=p\n"
    )
    annotation = run_validation.parse_annotation(gff)
    cds = run_validation.find_tag(annotation, "X_0001")
    # On the minus strand the stop codon is at the lower coordinate
    assert cds.stop_key == (-1, 10)
    assert cds.product == "p"


def test_stop_key_ignores_start_site():
    plus = run_validation.stop_key(100, 400, 1)
    assert plus == run_validation.stop_key(160, 400, 1)
    minus = run_validation.stop_key(100, 400, -1)
    assert minus == run_validation.stop_key(100, 340, -1)
