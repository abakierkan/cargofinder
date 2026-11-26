"""Command-line interface (CLI) for cargofinder."""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path
from typing import Optional

import cargofinder
from cargofinder.config import EXTENDED_SHELL_PROFILES, load_weights, select_profiles
from cargofinder.pipeline import run_pipeline
from cargofinder.profiles import (
    fetch_cargo_class_profiles,
    fetch_default_profiles,
    fetch_phage_marker_profiles,
)


def build_parser() -> argparse.ArgumentParser:
    """Build command-line parser."""
    parser = argparse.ArgumentParser(
        prog="cargofinder",
        description=(
            "cargofinder: Discovery and ranking of encapsulin cargo candidates "
            "from prokaryotic genomes."
        ),
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {cargofinder.__version__}",
    )

    subparsers = parser.add_subparsers(dest="command", help="Subcommand to execute")

    # Subcommand: run
    run_parser = subparsers.add_parser(
        "run",
        help="Scan a genome for encapsulin shells and rank candidate cargo enzymes.",
    )
    run_parser.add_argument(
        "input",
        type=Path,
        nargs="+",
        help=(
            "Nucleotide genome FASTA file(s), or one protein FASTA if --gff is given. "
            "With several genomes, cluster recurrence is measured across them."
        ),
    )
    run_parser.add_argument(
        "--gff",
        "-g",
        type=Path,
        default=None,
        help="GFF3 annotation file. If provided, input is treated as a protein FASTA.",
    )
    run_parser.add_argument(
        "--out",
        "-o",
        type=Path,
        default=Path("results"),
        help="Output directory (default: results).",
    )
    run_parser.add_argument(
        "--window",
        "-w",
        type=int,
        default=5,
        help="Number of genes upstream/downstream of shell to evaluate (default: 5).",
    )
    run_parser.add_argument(
        "--tail",
        "-t",
        type=int,
        default=20,
        help="Length of terminal peptides to extract in residues (default: 20).",
    )
    run_parser.add_argument(
        "--max-gap",
        type=int,
        default=150,
        help="Maximum gap in bp for likely operon co-transcription (default: 150).",
    )
    run_parser.add_argument(
        "--evalue",
        type=float,
        default=1e-5,
        help="Maximum E-value cutoff for encapsulin shell detection (default: 1e-5).",
    )
    run_parser.add_argument(
        "--profiles",
        "-p",
        type=Path,
        default=None,
        help="Directory containing shell HMM profile models (.hmm).",
    )
    run_parser.add_argument(
        "--motifs",
        "-m",
        type=Path,
        default=None,
        help="JSON file containing user-defined targeting peptide regex patterns.",
    )
    run_parser.add_argument(
        "--weights",
        type=str,
        default=None,
        help="Custom scoring weights as JSON string or path to JSON file.",
    )
    run_parser.add_argument(
        "--termini",
        choices=["c", "both"],
        default="c",
        help=(
            "Terminus used in the tail score: C-terminal only (default) or the better "
            "of both, for N-terminally targeted cargo."
        ),
    )
    run_parser.add_argument(
        "--min-shell-length",
        type=int,
        default=150,
        help="Minimum shell protein length in residues (default: 150).",
    )
    run_parser.add_argument(
        "--mmseqs-bin",
        type=Path,
        default=None,
        help="Explicit path to MMseqs2 executable.",
    )
    run_parser.add_argument(
        "--drop-prophage-loci",
        action="store_true",
        help=(
            "Exclude shell loci whose neighborhood contains two or more distinct "
            "phage marker families. Encapsulin shells share the HK97 fold with phage "
            "capsids, so some hits are prophage capsids. Off by default: flagged loci "
            "are reported either way, in shells.tsv and on every candidate row."
        ),
    )
    run_parser.add_argument(
        "--no-plot",
        action="store_true",
        help="Skip generating neighborhood visual diagrams.",
    )
    run_parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Enable verbose / debug logging.",
    )

    # Subcommand: fetch-profiles
    fetch_parser = subparsers.add_parser(
        "fetch-profiles",
        help="Download encapsulin shell HMM profiles (default: PF04454) from InterPro.",
    )
    fetch_parser.add_argument(
        "--extra",
        nargs="*",
        default=[],
        metavar="PFAM_ID",
        help=(
            "Additional profiles to fetch; choose from "
            f"{', '.join(sorted(EXTENDED_SHELL_PROFILES))}."
        ),
    )
    fetch_parser.add_argument(
        "--out",
        "-o",
        type=Path,
        default=Path("profiles"),
        help="Output directory for HMM profiles (default: profiles).",
    )
    fetch_parser.add_argument(
        "--timeout",
        type=int,
        default=30,
        help="Network timeout in seconds (default: 30).",
    )
    fetch_parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Enable verbose / debug logging.",
    )

    return parser


def main(argv: Optional[list[str]] = None) -> int:
    """CLI entrypoint."""
    parser = build_parser()
    args = parser.parse_args(argv)

    if not args.command:
        parser.print_help()
        return 0

    log_level = logging.DEBUG if getattr(args, "verbose", False) else logging.INFO
    logging.basicConfig(
        level=log_level,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%H:%M:%S",
    )

    if args.command == "fetch-profiles":
        print(f"Downloading encapsulin shell profiles to '{args.out}'...")
        try:
            files = fetch_default_profiles(
                args.out, select_profiles(args.extra), timeout=args.timeout
            )
            print(f"Successfully downloaded {len(files)} profile(s):")
            for f in files:
                print(f"  -> {f}")
            marker_files = fetch_phage_marker_profiles(args.out, timeout=args.timeout)
            print(
                f"Downloaded {len(marker_files)} phage marker profile(s) for "
                f"prophage screening into '{args.out / 'phage_markers'}'."
            )
            cargo_files = fetch_cargo_class_profiles(args.out, timeout=args.timeout)
            print(
                f"Downloaded {len(cargo_files)} cargo class profile(s) into "
                f"'{args.out / 'cargo_classes'}'."
            )
            print(f"Manifest written to: {args.out / 'manifest.json'}")
            return 0
        except Exception as exc:
            print(f"Error fetching profiles: {exc}", file=sys.stderr)
            return 1

    if args.command == "run":
        try:
            weights = load_weights(args.weights) if args.weights else None
        except Exception as exc:
            print(f"Configuration error: {exc}", file=sys.stderr)
            return 1

        print(f"Running cargofinder on: {', '.join(str(p) for p in args.input)}")
        try:
            res = run_pipeline(
                input_path=args.input,
                gff_path=args.gff,
                out_dir=args.out,
                profile_dir=args.profiles,
                motif_file=args.motifs,
                window_size=args.window,
                tail_length=args.tail,
                max_operon_gap_bp=args.max_gap,
                evalue_cutoff=args.evalue,
                weights=weights,
                min_shell_length_aa=args.min_shell_length,
                termini=args.termini,
                mmseqs_bin=args.mmseqs_bin,
                generate_plots=not args.no_plot,
                drop_prophage_loci=args.drop_prophage_loci,
            )

            print("\nScan completed:")
            print(f"  - Shells detected:       {len(res.shells)}")
            print(f"  - Candidates ranked:     {len(res.candidates)}")
            classed = sum(1 for c in res.candidates if c.cargo_class)
            if classed:
                print(f"  - Known cargo classes:   {classed} candidate(s)")
            flagged = sorted(
                shell_id
                for shell_id, ev in res.prophage.items()
                if ev.is_probable_prophage
            )
            if flagged:
                print(
                    f"  - Probable prophage loci: {len(flagged)} ({', '.join(flagged)})"
                )
                if not args.drop_prophage_loci:
                    print(
                        "    These loci carry two or more distinct phage marker "
                        "families nearby and may be phage capsids rather than "
                        "encapsulins. Re-run with --drop-prophage-loci to exclude "
                        "them."
                    )
            print(f"  - Results directory:     {args.out}")
            if "tsv" in res.output_files:
                print(f"  - Ranked TSV:            {res.output_files['tsv']}")
            if not res.shells:
                print(
                    "\nNo shell found with the profiles used. This does not mean the "
                    "genome has no encapsulin: the default profile covers only the "
                    "Linocin_M18 family (see `fetch-profiles --extra`)."
                )
            if res.candidates:
                top = res.candidates[0]
                print(
                    f"\nTop Candidate: {top.candidate_id} "
                    f"(Score: {top.score_total:.3f}, Shell: {top.shell_gene_id}, "
                    f"Dist: {top.gene_distance:+d})"
                )
            return 0
        except Exception as exc:
            print(f"\nExecution error: {exc}", file=sys.stderr)
            if getattr(args, "verbose", False):
                import traceback

                traceback.print_exc()
            return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
