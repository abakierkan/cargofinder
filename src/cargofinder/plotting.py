"""Publication-ready gene neighborhood plotting with matplotlib."""

from __future__ import annotations

import logging
import os
import tempfile
from pathlib import Path
from typing import List, Optional

# Configure headless backend before importing pyplot
if "MPLCONFIGDIR" not in os.environ:
    os.environ["MPLCONFIGDIR"] = tempfile.gettempdir()

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.patches as patches  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

from cargofinder.io import Gene, GenomeData  # noqa: E402
from cargofinder.scoring import ScoredCandidate  # noqa: E402
from cargofinder.shells import ShellHit  # noqa: E402

logger = logging.getLogger(__name__)

# Color palette: colorblind-safe, clean academic styling
COLOR_SHELL = "#1f77b4"  # Deep blue for encapsulin shell
COLOR_TOP_CARGO = "#d95f02"  # Vibrant vermilion/orange for top candidate
COLOR_CARGO_OTHER = "#7570b3"  # Muted purple for other candidate cargo
COLOR_NON_CANDIDATE = "#d0d7de"  # Soft neutral gray for background genes

EDGE_SHELL = "#0b3c61"
EDGE_TOP_CARGO = "#993f00"
EDGE_CARGO_OTHER = "#4c4878"
EDGE_NON_CANDIDATE = "#6e7781"


def _clean_gene_label(gene_id: str, contig_id: str) -> str:
    """Format a clean, readable gene label to prevent visual crowding."""
    if gene_id.startswith(f"{contig_id}_"):
        suffix = gene_id[len(contig_id) + 1 :]
        return f"gene_{suffix}"
    if len(gene_id) > 18 and "_" in gene_id:
        return f"..._{gene_id.split('_')[-1]}"
    return gene_id


def _draw_gene_arrow(
    ax: plt.Axes,
    start: int,
    end: int,
    strand: int,
    y: float = 0.45,
    height: float = 0.34,
    color: str = "#d0d7de",
    edgecolor: str = "#6e7781",
) -> None:
    """Draw a directional chevron/pentagon arrow representing a gene."""
    length = max(1, end - start)
    arrow_len = min(length * 0.35, 300)

    y_bottom = y - height / 2.0
    y_top = y + height / 2.0

    if length <= arrow_len:
        if strand == 1:
            verts = [(start, y_bottom), (end, y), (start, y_top)]
        else:
            verts = [(end, y_bottom), (start, y), (end, y_top)]
    else:
        if strand == 1:
            body_end = end - arrow_len
            verts = [
                (start, y_bottom),
                (body_end, y_bottom),
                (end, y),
                (body_end, y_top),
                (start, y_top),
            ]
        else:
            body_start = start + arrow_len
            verts = [
                (end, y_bottom),
                (body_start, y_bottom),
                (start, y),
                (body_start, y_top),
                (end, y_top),
            ]

    poly = patches.Polygon(
        verts,
        closed=True,
        facecolor=color,
        edgecolor=edgecolor,
        linewidth=1.1,
        zorder=3,
    )
    ax.add_patch(poly)


def plot_neighborhood(
    shell: ShellHit,
    genome: GenomeData,
    candidates: List[ScoredCandidate],
    out_path: Path | str,
    flank_genes: int = 6,
    title: Optional[str] = None,
) -> Path:
    """Generate a clean genomic neighborhood diagram centered on an encapsulin shell.

    Parameters
    ----------
    shell : ShellHit
        The detected encapsulin shell gene.
    genome : GenomeData
        Annotated genome containing full contig gene order.
    candidates : list of ScoredCandidate
        Ranked candidates around this shell.
    out_path : Path or str
        Destination image path (PNG or SVG).
    flank_genes : int
        Number of genes to display on each side of the shell (default 6).
    title : str, optional
        Custom title. Defaults to shell ID and contig coordinates.

    Returns
    -------
    Path
        Path to the saved figure.
    """
    out_file = Path(out_path)
    out_file.parent.mkdir(parents=True, exist_ok=True)

    contig_genes = genome.genes_on_contig(shell.contig_id)
    shell_idx = -1
    for idx, g in enumerate(contig_genes):
        if g.gene_id == shell.gene_id:
            shell_idx = idx
            break

    if shell_idx == -1:
        raise ValueError(f"Shell {shell.gene_id} not found on contig {shell.contig_id}")

    min_idx = max(0, shell_idx - flank_genes)
    max_idx = min(len(contig_genes) - 1, shell_idx + flank_genes)
    visible_genes: List[Gene] = contig_genes[min_idx : max_idx + 1]

    # Map candidate IDs to their ScoredCandidate info
    cand_by_id = {
        c.candidate_id: c for c in candidates if c.shell_gene_id == shell.gene_id
    }
    top_cand_id = (
        min(cand_by_id.values(), key=lambda c: c.rank).candidate_id
        if cand_by_id
        else None
    )

    # Calculate coordinate bounds
    min_coord = min(g.start for g in visible_genes)
    max_coord = max(g.end for g in visible_genes)
    coord_padding = max(500, int((max_coord - min_coord) * 0.05))
    x_min = min_coord - coord_padding
    x_max = max_coord + coord_padding

    fig, ax = plt.subplots(figsize=(11.5, 4.4), dpi=200)

    # Contig backbone line
    ax.axhline(y=0.45, color="#8c959f", linestyle="-", linewidth=1.5, zorder=1)

    for idx, g in enumerate(visible_genes):
        is_shell = g.gene_id == shell.gene_id
        cand_info = cand_by_id.get(g.gene_id)

        if is_shell:
            color = COLOR_SHELL
            edge_color = EDGE_SHELL
            badge = "[Shell]"
            fontweight = "bold"
        elif cand_info is not None:
            if g.gene_id == top_cand_id:
                color = COLOR_TOP_CARGO
                edge_color = EDGE_TOP_CARGO
                badge = f"Rank #1 ({cand_info.score_total:.2f})"
                fontweight = "bold"
            else:
                color = COLOR_CARGO_OTHER
                edge_color = EDGE_CARGO_OTHER
                badge = f"Rank #{cand_info.rank}"
                fontweight = "normal"
        else:
            color = COLOR_NON_CANDIDATE
            edge_color = EDGE_NON_CANDIDATE
            badge = ""
            fontweight = "normal"

        _draw_gene_arrow(
            ax=ax,
            start=g.start,
            end=g.end,
            strand=g.strand,
            y=0.45,
            height=0.34,
            color=color,
            edgecolor=edge_color,
        )

        clean_id = _clean_gene_label(g.gene_id, shell.contig_id)
        full_text = f"{clean_id}\n{badge}" if badge else clean_id
        mid_x = (g.start + g.end) / 2.0

        # 2-tier alternating tick height to ensure adjacent genes never collide
        tier = idx % 2
        y_tick_top = 0.68 if tier == 0 else 0.82
        y_text = y_tick_top + 0.02

        ax.plot(
            [mid_x, mid_x],
            [0.63, y_tick_top],
            color="#8c959f",
            linewidth=0.75,
            linestyle="-",
            zorder=2,
        )

        ax.text(
            mid_x,
            y_text,
            full_text,
            ha="left",
            va="bottom",
            fontsize=8,
            fontweight=fontweight,
            rotation=36,
            rotation_mode="anchor",
            zorder=4,
        )

    # Configure axes
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(-0.05, 1.55)
    ax.set_yticks([])
    ax.set_xlabel(
        f"Genomic position on {shell.contig_id} (bp)",
        fontsize=9.5,
        fontweight="bold",
        labelpad=8,
    )

    # Format ticks with commas
    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f"{int(x):,}"))

    if title is None:
        strand_s = "(+)" if shell.strand == 1 else "(-)"
        clean_shell_id = _clean_gene_label(shell.gene_id, shell.contig_id)
        title = (
            f"Encapsulin Locus: {clean_shell_id} {strand_s} | "
            f"{shell.contig_id}:{shell.start:,}..{shell.end:,}"
        )
    ax.set_title(title, fontsize=11, pad=14, fontweight="bold")

    # Legend elements
    legend_elements = [
        patches.Patch(
            facecolor=COLOR_SHELL, edgecolor=EDGE_SHELL, label="Encapsulin shell"
        ),
        patches.Patch(
            facecolor=COLOR_TOP_CARGO,
            edgecolor=EDGE_TOP_CARGO,
            label="Top candidate cargo (Rank #1)",
        ),
        patches.Patch(
            facecolor=COLOR_CARGO_OTHER,
            edgecolor=EDGE_CARGO_OTHER,
            label="Other candidate cargo",
        ),
        patches.Patch(
            facecolor=COLOR_NON_CANDIDATE,
            edgecolor=EDGE_NON_CANDIDATE,
            label="Flanking gene",
        ),
    ]
    ax.legend(
        handles=legend_elements,
        loc="upper right",
        bbox_to_anchor=(1.0, -0.22),
        ncol=4,
        fontsize=8,
        frameon=False,
    )

    # Style borders
    for spine in ["top", "left", "right"]:
        ax.spines[spine].set_visible(False)
    ax.spines["bottom"].set_color("#666666")

    plt.tight_layout()
    fig.savefig(out_file, dpi=300, bbox_inches="tight")
    plt.close(fig)

    logger.info("Saved neighborhood diagram to %s", out_file)
    return out_file
