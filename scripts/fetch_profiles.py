#!/usr/bin/env python3
"""Helper script to download encapsulin shell HMM profiles from InterPro."""

import argparse
import logging
import sys
from pathlib import Path

# Add src to path so script can be run without pip install
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from cargofinder.config import DEFAULT_SHELL_PROFILES  # noqa: E402
from cargofinder.profiles import fetch_default_profiles  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Download encapsulin shell HMM profiles from InterPro/Pfam."
    )
    parser.add_argument(
        "--out",
        "-o",
        type=Path,
        default=Path("profiles"),
        help="Output directory to save HMM profiles and manifest (default: profiles)",
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=30,
        help="Network timeout in seconds (default: 30)",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Enable debug logging",
    )
    args = parser.parse_args()

    log_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(level=log_level, format="%(levelname)s: %(message)s")

    n_prof = len(DEFAULT_SHELL_PROFILES)
    print(f"Fetching {n_prof} default shell profiles into '{args.out}'...")
    for p in DEFAULT_SHELL_PROFILES:
        print(f"  - {p.pfam_id}: {p.family} ({p.description})")

    try:
        files = fetch_default_profiles(args.out, timeout=args.timeout)
        print(f"Successfully downloaded {len(files)} profile(s):")
        for f in files:
            print(f"  -> {f}")
        print(f"Manifest written to: {args.out / 'manifest.json'}")
        return 0
    except Exception as exc:
        print(f"Error fetching profiles: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
