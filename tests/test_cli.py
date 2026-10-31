"""CLI smoke tests and end-to-end command line validation."""

from pathlib import Path

import pandas as pd
import pytest

from cargofinder.cli import main
from cargofinder.pipeline import run_pipeline

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
    assert (out_dir / "shells.tsv").exists()

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


def test_cli_protein_fasta_without_gff_is_an_error(tmp_path, capsys):
    faa = FIXTURES_DIR / "thermotoga_slice.faa"
    ret = main(
        ["run", str(faa), "--out", str(tmp_path), "--profiles", str(PROFILES_DIR)]
    )
    assert ret != 0
    assert "--gff" in capsys.readouterr().err


def test_cli_explicit_empty_profile_dir_is_an_error(tmp_path, capsys):
    fna = FIXTURES_DIR / "thermotoga_slice.fna"
    empty = tmp_path / "empty"
    empty.mkdir()
    ret = main(["run", str(fna), "--out", str(tmp_path), "--profiles", str(empty)])
    assert ret != 0
    assert "No .hmm" in capsys.readouterr().err


@pytest.mark.parametrize(
    "flag", [("--window", "0"), ("--tail", "0"), ("--window", "-2")]
)
def test_cli_rejects_bad_parameters(tmp_path, capsys, flag):
    fna = FIXTURES_DIR / "thermotoga_slice.fna"
    ret = main(
        [
            "run",
            str(fna),
            "--out",
            str(tmp_path),
            "--profiles",
            str(PROFILES_DIR),
            *flag,
        ]
    )
    assert ret != 0
    assert "at least 1" in capsys.readouterr().err


def test_cli_unknown_weight_key_is_an_error(tmp_path, capsys):
    fna = FIXTURES_DIR / "thermotoga_slice.fna"
    ret = main(["run", str(fna), "--out", str(tmp_path), "--weights", '{"x": 1}'])
    assert ret != 0
    assert "Unknown scoring weight" in capsys.readouterr().err


def test_cli_several_inputs_with_gff_is_an_error(tmp_path, capsys):
    faa = FIXTURES_DIR / "thermotoga_slice.faa"
    gff = FIXTURES_DIR / "thermotoga_slice.gff3"
    ret = main(["run", str(faa), str(faa), "--gff", str(gff), "--out", str(tmp_path)])
    assert ret != 0
    assert "one protein FASTA" in capsys.readouterr().err


def test_single_shell_run_needs_no_mmseqs(tmp_path):
    fna = FIXTURES_DIR / "thermotoga_slice.fna"
    res = run_pipeline(
        fna,
        out_dir=tmp_path,
        profile_dir=PROFILES_DIR,
        generate_plots=False,
        mmseqs_bin="/nonexistent/mmseqs",
    )
    assert len(res.shells) == 1
    assert res.summary["clustering"].startswith("skipped")
    assert res.summary["environment"]["mmseqs2"] == "not used"
    assert res.summary["profiles"][0]["sha256"]


@pytest.mark.external
def test_multi_genome_run_measures_recurrence(tmp_path):
    fna = FIXTURES_DIR / "thermotoga_slice.fna"
    copy = tmp_path / "copy_of_slice.fna"
    copy.write_bytes(fna.read_bytes())
    ret = main(
        [
            "run",
            str(fna),
            str(copy),
            "--out",
            str(tmp_path / "out"),
            "--profiles",
            str(PROFILES_DIR),
            "--no-plot",
        ]
    )
    assert ret == 0
    df = pd.read_csv(tmp_path / "out" / "candidates.tsv", sep="\t")
    assert set(df["genome_id"]) == {"thermotoga_slice", "copy_of_slice"}
    assert df["candidate_id"].str.contains("|", regex=False).all()
    assert (df["cluster_shell_count"] == 2).all()
    assert set(df["rank_in_shell"]) >= {1}
