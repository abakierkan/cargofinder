"""Clustering candidate cargo sequences using MMseqs2 to detect cross-system recurrence.

Biological rationale:
Encapsulin shell genes frequently co-occur with homologous cargo enzymes
(e.g. ferritin-like proteins, DyP peroxidases) across diverse bacterial and
archaeal phyla. When analyzing multiple encapsulin loci or multiple genomes,
recurrent appearance of homologous neighbor families next to independent shell
genes provides strong evolutionary evidence for bona fide functional cargo.
"""

from __future__ import annotations

import logging
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional

from cargofinder.neighborhood import NeighborCandidate

logger = logging.getLogger(__name__)


@dataclass
class ClusterResult:
    """Clustering output and recurrence metrics across candidate cargo genes."""

    # candidate_id -> cluster representative ID
    clusters: Dict[str, str]
    # cluster_id -> total number of candidate sequences in cluster
    cluster_sizes: Dict[str, int]
    # cluster_id -> number of distinct encapsulin shell genes where this cluster appears
    cluster_shell_counts: Dict[str, int]
    # cluster_id -> number of distinct genomes where this cluster appears
    cluster_genome_counts: Dict[str, int]

    def get_cluster_id(self, candidate_id: str) -> str:
        return self.clusters.get(candidate_id, candidate_id)

    def get_shell_recurrence(self, candidate_id: str) -> int:
        cid = self.get_cluster_id(candidate_id)
        return self.cluster_shell_counts.get(cid, 1)

    def get_cluster_size(self, candidate_id: str) -> int:
        cid = self.get_cluster_id(candidate_id)
        return self.cluster_sizes.get(cid, 1)


def find_mmseqs_binary(custom_path: Optional[str | Path] = None) -> Path:
    """Locate the MMseqs2 executable or raise an informative error.

    Checks:
    1. User-supplied custom path (if provided)
    2. System PATH (shutil.which)
    3. Python virtual environment / conda / pixi bin directory (sys.prefix/bin/mmseqs)
    """
    if custom_path:
        cp = Path(custom_path)
        if cp.is_file():
            return cp
        raise FileNotFoundError(
            f"Specified MMseqs2 executable does not exist: {custom_path}"
        )

    on_path = shutil.which("mmseqs")
    if on_path:
        return Path(on_path)

    env_bin = Path(sys.prefix) / "bin" / "mmseqs"
    if env_bin.is_file():
        return env_bin

    raise RuntimeError(
        "MMseqs2 executable ('mmseqs') was not found on PATH or in current env.\n"
        "MMseqs2 is required to cluster neighbor sequences across encapsulin loci.\n"
        "To install MMseqs2:\n"
        "  - Pixi: pixi add mmseqs2\n"
        "  - Conda/Mamba: conda install -c conda-forge -c bioconda mmseqs2\n"
        "  - Ubuntu/Debian: apt-get install mmseqs2\n"
        "Alternatively, specify the binary path using --mmseqs-bin."
    )


def cluster_candidates(
    candidates: List[NeighborCandidate],
    min_seq_id: float = 0.3,
    coverage: float = 0.8,
    mmseqs_bin: Optional[str | Path] = None,
    work_dir: Optional[Path | str] = None,
) -> ClusterResult:
    """Cluster candidate cargo proteins using MMseqs2.

    Parameters
    ----------
    candidates : list of NeighborCandidate
        Extracted candidate genes across all shells in the run.
    min_seq_id : float
        Minimum sequence identity threshold (0.0 to 1.0, default 0.3).
    coverage : float
        Minimum coverage threshold (0.0 to 1.0, default 0.8).
    mmseqs_bin : Path or str, optional
        Explicit path to mmseqs binary.
    work_dir : Path or str, optional
        Directory for intermediate files. If None, uses a temporary directory.

    Returns
    -------
    ClusterResult
        Clustering assignments and recurrence metrics.
    """
    if not candidates:
        logger.debug("No candidates provided for clustering; returning empty result")
        return ClusterResult(
            clusters={},
            cluster_sizes={},
            cluster_shell_counts={},
            cluster_genome_counts={},
        )

    # If only 1 candidate, clustering is trivial
    if len(candidates) == 1:
        c = candidates[0]
        return ClusterResult(
            clusters={c.candidate_id: c.candidate_id},
            cluster_sizes={c.candidate_id: 1},
            cluster_shell_counts={c.candidate_id: 1},
            cluster_genome_counts={c.candidate_id: 1},
        )

    binary = find_mmseqs_binary(custom_path=mmseqs_bin)

    with tempfile.TemporaryDirectory(dir=work_dir, prefix="cargofinder_mmseqs_") as tmp:
        tmp_path = Path(tmp)
        in_fasta = tmp_path / "candidates.faa"
        out_prefix = tmp_path / "clu"
        mmseqs_tmp = tmp_path / "tmp"
        mmseqs_tmp.mkdir(parents=True, exist_ok=True)

        # Write candidate protein sequences (deduplicating identical candidate_ids)
        written_ids = set()
        with in_fasta.open("w", encoding="utf-8") as fh:
            for cand in candidates:
                cid = cand.candidate_id
                if cid not in written_ids and cand.candidate_protein_seq:
                    fh.write(f">{cid}\n{cand.candidate_protein_seq}\n")
                    written_ids.add(cid)

        cmd = [
            str(binary),
            "easy-cluster",
            str(in_fasta),
            str(out_prefix),
            str(mmseqs_tmp),
            "--min-seq-id",
            str(min_seq_id),
            "-c",
            str(coverage),
            "-v",
            "1",  # quiet logging
        ]

        logger.debug("Running MMseqs2 command: %s", " ".join(cmd))
        try:
            res = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=False,
            )
        except OSError as exc:
            raise RuntimeError(f"Failed to execute MMseqs2 ({binary}): {exc}") from exc

        tsv_path = tmp_path / "clu_cluster.tsv"
        if res.returncode != 0 or not tsv_path.exists():
            stderr_msg = res.stderr.strip() if res.stderr else "Unknown error"
            raise RuntimeError(
                f"MMseqs2 easy-cluster failed (exit {res.returncode}):\n{stderr_msg}"
            )

        # Parse TSV: col 0 = representative_id, col 1 = member_id
        clusters: Dict[str, str] = {}
        with tsv_path.open("r", encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                parts = line.split("\t")
                if len(parts) >= 2:
                    rep_id, member_id = parts[0], parts[1]
                    clusters[member_id] = rep_id

        # Fallback for any candidate that might not have appeared in TSV
        for cand in candidates:
            clusters.setdefault(cand.candidate_id, cand.candidate_id)

    # Compute recurrence statistics
    cluster_members: Dict[str, List[NeighborCandidate]] = {}
    for cand in candidates:
        cid = clusters.get(cand.candidate_id, cand.candidate_id)
        cluster_members.setdefault(cid, []).append(cand)

    cluster_sizes = {cid: len(mem) for cid, mem in cluster_members.items()}
    cluster_shell_counts = {
        cid: len(set(m.shell_gene_id for m in mem))
        for cid, mem in cluster_members.items()
    }
    cluster_genome_counts = {
        cid: len(set(m.contig_id for m in mem)) for cid, mem in cluster_members.items()
    }

    logger.info(
        "MMseqs2 clustered %d candidate(s) into %d distinct cluster(s)",
        len(candidates),
        len(cluster_sizes),
    )
    return ClusterResult(
        clusters=clusters,
        cluster_sizes=cluster_sizes,
        cluster_shell_counts=cluster_shell_counts,
        cluster_genome_counts=cluster_genome_counts,
    )
