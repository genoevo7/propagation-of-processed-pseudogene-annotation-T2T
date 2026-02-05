#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
gapmasked_alignment_fel_meme_plot.py

Purpose
-------
Visualize codon-level selection signals (FEL and MEME) on a codon-based CDS alignment
by drawing alignment-aware (gap-masked) bars per species and overlaying:
- FEL purifying sites (vertical black ticks)
- MEME episodic diversifying sites (red dots)

IMPORTANT
---------
This script does NOT perform selection inference. FEL and MEME are assumed to have been
run externally (e.g., Datamonkey/HyPhy), and their site tables are provided as Excel files.

Inputs (default filenames)
-------------------------
- Alignment FASTA (codon-based, stop codon removed):
    septin14_cds_alignment_no_stop.fas

- FEL results (Excel export):
    fel_results.xlsx
  Required columns:
    - codon (1-based)
    - p-value
    - class  (e.g., Purifying)

- MEME results (Excel export):
    meme_results.xlsx
  Required columns:
    - Codon (1-based)
    - p-value
    - Class (e.g., Diversifying)

Outputs
-------
- Figure_SEPTIN14_gapmasked_FEL_MEME.png
- Figure_SEPTIN14_gapmasked_FEL_MEME.pdf
"""

import re
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import matplotlib.lines as mlines
import numpy as np
from Bio import SeqIO

# ---------------------------
# Global style (paper-ready)
# ---------------------------
plt.rcParams["font.family"] = "Arial"
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42


def extract_species_name(record):
    """
    Extract species name from FASTA header if possible.

    Supported patterns:
    - organism=Species name
    - [Species name]
    Fallback: record.id
    """
    desc = record.description
    m1 = re.search(r"organism=([A-Za-z_ ]+)", desc)
    if m1:
        return m1.group(1).strip()
    m2 = re.search(r"\[([^\]]+)\]", desc)
    if m2:
        return m2.group(1).strip()
    return record.id


def read_fel_sites(xlsx, p):
    """
    Read FEL sites from an Excel file and return significant purifying codons.

    Required columns:
      - codon (1-based)
      - p-value
      - class
    """
    df = pd.read_excel(xlsx)
    return sorted(
        df[(df["class"].astype(str).str.lower() == "purifying") & (df["p-value"] <= p)]
        ["codon"].astype(int).unique()
    )


def read_meme_sites(xlsx, p):
    """
    Read MEME sites from an Excel file and return significant episodic diversifying codons.

    Required columns:
      - Codon (1-based)
      - p-value
      - Class
    """
    df = pd.read_excel(xlsx)
    return sorted(
        df[(df["Class"].astype(str).str.lower() == "diversifying") & (df["p-value"] <= p)]
        ["Codon"].astype(int).unique()
    )


def parse_alignment(seq):
    """
    Validate codon-based alignment and return:
      - total length (bp)
      - list of gap codon blocks as (start_bp, width_bp) for codons containing '-'
    """
    if len(seq) % 3 != 0:
        raise ValueError("Alignment length is not divisible by 3 (not codon-based).")
    gaps = []
    for i in range(0, len(seq), 3):
        if "-" in seq[i:i + 3]:
            gaps.append((i, 3))
    return len(seq), gaps


def main():
    # ---------------------------
    # Inputs / parameters
    # ---------------------------
    aln = "septin14_cds_alignment_no_stop.fas"
    fel = "fel_results.xlsx"
    meme = "meme_results.xlsx"

    out_base = "Figure_SEPTIN14_gapmasked_FEL_MEME"
    pval = 0.05

    # Highlight specific bp positions (alignment-based coordinates)
    highlight_bp = [1123, 1299]

    # ---------------------------
    # Load data
    # ---------------------------
    records = list(SeqIO.parse(aln, "fasta"))
    species = [extract_species_name(r) for r in records]

    fel_sites = read_fel_sites(fel, pval)
    meme_sites = read_meme_sites(meme, pval)

    parsed = [parse_alignment(str(r.seq)) for r in records]
    max_len = max(t for t, _ in parsed)

    # ---------------------------
    # Plot
    # ---------------------------
    fig_h = max(6, 0.38 * len(records))
    fig, ax = plt.subplots(figsize=(18, fig_h), dpi=600)

    cmap = plt.get_cmap("tab20")
    bar_h = 0.78
    y_pos = np.arange(len(records))[::-1]

    for i, (y, (total_len, gaps)) in enumerate(zip(y_pos, parsed)):
        # Background bar
        ax.add_patch(
            Rectangle(
                (0, y - bar_h / 2),
                total_len,
                bar_h,
                facecolor=cmap(i % 20),
                edgecolor="black",
                alpha=0.22,
                linewidth=0.6
            )
        )

        # Gap-masked codons
        for start, width in gaps:
            ax.add_patch(
                Rectangle(
                    (start, y - bar_h / 2),
                    width,
                    bar_h,
                    facecolor="white",
                    edgecolor=None,
                    zorder=3
                )
            )

        # FEL purifying sites as vertical ticks (codon -> bp start)
        xs_fel = [(c - 1) * 3 for c in fel_sites if (c - 1) * 3 < total_len]
        ax.vlines(xs_fel, y - bar_h / 2, y + bar_h / 2, color="black", linewidth=0.7, alpha=0.9)

        # MEME episodic diversifying sites as red dots (codon -> bp start)
        xs_meme = [(c - 1) * 3 for c in meme_sites if (c - 1) * 3 < total_len]
        ax.scatter(xs_meme, [y] * len(xs_meme), color="red", s=20, zorder=4)

    # ---------------------------
    # Highlight specific bp positions (labels below the bars)
    # ---------------------------
    label_y = y_pos[-1] - bar_h * 2

    for bp in highlight_bp:
        if bp < max_len:
            ax.axvline(bp, color="blue", linestyle="--", linewidth=1.0, alpha=0.6, zorder=2)
            ax.text(
                bp, label_y,
                f"{bp} bp",
                ha="center", va="top",
                fontsize=12, color="blue", fontweight="bold"
            )

    ax.set_yticks(y_pos)
    ax.set_yticklabels(species, fontsize=14, fontweight="bold")
    ax.set_xlim(0, max_len * 1.01)

    ax.set_xlabel("CDS position (bp; alignment-based, gaps masked)", fontsize=14, fontweight="bold")
    ax.set_title("SEPTIN14: alignment-aware CDS with gap-masked bars", fontsize=14, fontweight="bold")

    # Legend
    l1 = mlines.Line2D([], [], color="black", linewidth=2, label="FEL purifying sites")
    l2 = mlines.Line2D([], [], color="red", marker="o", linestyle="None", markersize=6,
                       label="MEME episodic diversifying sites")

    ax.legend(
        handles=[l1, l2],
        frameon=False,
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        fontsize=12,
    )

    ax.grid(axis="x", alpha=0.12)
    fig.subplots_adjust(left=0.18, right=0.82, top=0.90, bottom=0.12)

    # ---------------------------
    # Save (PNG + PDF)
    # ---------------------------
    fig.savefig(f"{out_base}.png", dpi=600)
    fig.savefig(f"{out_base}.pdf")
    plt.close()

    print("[OK] saved:")
    print(f" - {out_base}.png (dpi 600)")
    print(f" - {out_base}.pdf (vector)")


if __name__ == "__main__":
    main()
