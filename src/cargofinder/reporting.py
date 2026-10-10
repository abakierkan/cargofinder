"""Export and serialization utilities for ranked cargo candidates and run summaries."""

from __future__ import annotations

import json
import logging
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional

import pandas as pd
import pyhmmer
import pyrodigal

import cargofinder
from cargofinder.config import RunConfig
from cargofinder.io import GenomeData
from cargofinder.scoring import ScoredCandidate
from cargofinder.shells import ShellHit

logger = logging.getLogger(__name__)


def export_candidates_tsv(
    candidates: List[ScoredCandidate],
    out_path: Path | str,
) -> Path:
    """Export ranked candidates table to a TSV file.

    Parameters
    ----------
    candidates : list of ScoredCandidate
        Ranked candidates.
    out_path : Path or str
        Destination path.

    Returns
    -------
    Path
        Path to the written TSV.
    """
    out_file = Path(out_path)
    out_file.parent.mkdir(parents=True, exist_ok=True)

    if not candidates:
        logger.warning("No candidates to export; writing empty TSV table")
        out_file.write_text("rank\tcandidate_id\tshell_gene_id\tscore_total\n")
        return out_file

    rows = [c.to_dict() for c in candidates]
    df = pd.DataFrame(rows)
    df.to_csv(out_file, sep="\t", index=False)
    logger.info("Exported %d ranked candidate(s) to %s", len(df), out_file)
    return out_file


def export_candidates_json(
    candidates: List[ScoredCandidate],
    out_path: Path | str,
    metadata: Optional[Dict[str, Any]] = None,
) -> Path:
    """Export ranked candidates table to a structured JSON file."""
    out_file = Path(out_path)
    out_file.parent.mkdir(parents=True, exist_ok=True)

    data: Dict[str, Any] = {
        "metadata": metadata or {},
        "candidates": [c.to_dict() for c in candidates],
    }

    out_file.write_text(json.dumps(data, indent=2), encoding="utf-8")
    logger.info("Exported candidates JSON to %s", out_file)
    return out_file


def generate_run_summary(
    genome: GenomeData,
    shells: List[ShellHit],
    candidates: List[ScoredCandidate],
    config: RunConfig,
) -> Dict[str, Any]:
    """Compile counts, parameters, and environment versions for the run."""
    unique_clusters = len(set(c.cluster_id for c in candidates)) if candidates else 0

    return {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "cargofinder_version": cargofinder.__version__,
        "environment": {
            "python": sys.version.split()[0],
            "platform": platform.platform(),
            "pyhmmer": getattr(pyhmmer, "__version__", "unknown"),
            "pyrodigal": getattr(pyrodigal, "__version__", "unknown"),
        },
        "parameters": {
            "window_size": config.window_size,
            "tail_length": config.tail_length,
            "max_operon_gap_bp": config.max_operon_gap_bp,
            "evalue_cutoff": config.evalue_cutoff,
            "min_shell_length_aa": config.min_shell_length_aa,
            "scoring_weights": {
                "w_distance": config.weights.w_distance,
                "w_operon": config.weights.w_operon,
                "w_tail_features": config.weights.w_tail_features,
                "w_recurrence": config.weights.w_recurrence,
            },
        },
        "counts": {
            "genome_id": genome.genome_id,
            "total_genes": len(genome.genes),
            "shells_detected": len(shells),
            "candidates_extracted": len(candidates),
            "distinct_clusters": unique_clusters,
        },
    }


def save_run_summary(summary: Dict[str, Any], out_path: Path | str) -> Path:
    """Save the run summary to a JSON file and a human-readable text file."""
    out_file = Path(out_path)
    out_file.parent.mkdir(parents=True, exist_ok=True)

    out_file.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    # Also save a human-readable summary text file alongside JSON
    txt_file = out_file.with_suffix(".txt")
    lines = [
        "===========================================================",
        "cargofinder: Encapsulin Cargo Candidate Discovery Summary",
        "===========================================================",
        f"Genome ID:           {summary['counts']['genome_id']}",
        f"Total Genes:         {summary['counts']['total_genes']}",
        f"Shells Detected:     {summary['counts']['shells_detected']}",
        f"Candidates Ranked:   {summary['counts']['candidates_extracted']}",
        f"Distinct Clusters:   {summary['counts']['distinct_clusters']}",
        "-----------------------------------------------------------",
        f"Window Size:         +/- {summary['parameters']['window_size']} genes",
        f"Tail Length:         {summary['parameters']['tail_length']} residues",
        f"Max Operon Gap:      {summary['parameters']['max_operon_gap_bp']} bp",
        f"E-value Cutoff:      {summary['parameters']['evalue_cutoff']}",
        "-----------------------------------------------------------",
        f"cargofinder:         v{summary['cargofinder_version']}",
        f"pyhmmer:             v{summary['environment']['pyhmmer']}",
        f"pyrodigal:           v{summary['environment']['pyrodigal']}",
        f"Python:              {summary['environment']['python']}",
        "===========================================================",
    ]
    txt_file.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return out_file
