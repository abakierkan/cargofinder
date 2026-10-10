"""CLI smoke tests and end-to-end command line validation."""

from pathlib import Path

import pandas as pd
import pytest

from cargofinder.cli import main

FIXTURES_DIR = Path(__file__).parent / "data"
PROFILES_DIR = FIXTURES_DIR


def test_cli_help_and_version(capsys):
    with pytest.raises(SystemExit) as exc_info:
        main(["--help"])
    assert exc_info.value.code == 0
    captured = capsys.readouterr()
    assert "cargofinder" in captured.out
    assert "run" in captured.out
    assert "fetch-profiles" in captured.out

    with pytest.raises(SystemExit) as exc_run:
        main(["run", "--help"])
    assert exc_run.value.code == 0
    captured_run = capsys.readouterr()
    assert "--gff" in captured_run.out
    assert "--window" in captured_run.out


def test_cli_smoke_nucleotide_fasta(tmp_path):
    fna = FIXTURES_DIR / "thermotoga_slice.fna"
    out_dir = tmp_path / "results_fna"

    ret = main(
        [
            "run",
            str(fna),
            "--out",
            str(out_dir),
            "--profiles",
            str(PROFILES_DIR),
            "--window",
            "5",
            "--tail",
            "20",
        ]
    )
    assert ret == 0

    tsv_file = out_dir / "candidates.tsv"
    json_file = out_dir / "candidates.json"
    summary_file = out_dir / "run_summary.json"
    plot_file = out_dir / "plots" / "neighborhood_thermotoga_slice_0008.png"

    assert tsv_file.exists()
    assert json_file.exists()
    assert summary_file.exists()
    assert plot_file.exists()

    df = pd.read_csv(tsv_file, sep="\t")
    assert len(df) > 0
    # Top ranked candidate should be the known cargo (thermotoga_slice_0009)
    assert df.iloc[0]["candidate_id"] == "thermotoga_slice_0009"
    assert df.iloc[0]["rank"] == 1


def test_cli_smoke_protein_and_gff(tmp_path):
    faa = FIXTURES_DIR / "thermotoga_slice.faa"
    gff = FIXTURES_DIR / "thermotoga_slice.gff3"
    out_dir = tmp_path / "results_gff"

    ret = main(
        [
            "run",
            str(faa),
            "--gff",
            str(gff),
            "--out",
            str(out_dir),
            "--profiles",
            str(PROFILES_DIR),
        ]
    )
    assert ret == 0

    tsv_file = out_dir / "candidates.tsv"
    assert tsv_file.exists()
    df = pd.read_csv(tsv_file, sep="\t")
    assert df.iloc[0]["candidate_id"] == "thermotoga_slice_0009"


def test_cli_missing_input_error(tmp_path, capsys):
    ret = main(["run", str(tmp_path / "nonexistent.fa")])
    assert ret != 0
    captured = capsys.readouterr()
    assert "Execution error" in captured.err
