"""Generate the ConfPert overview / architecture figure (fig_overview.pdf).

A pure schematic (no data dependency): the end-to-end pipeline read left to
right -- predictors x datasets x splits -> six per-population distributional
discrepancies -> four off-the-shelf conformal head levels (finite-sample
coverage threshold) -> the three pre-registered knockout claims K1/K2/K3.

One takeaway: ConfPert wraps off-the-shelf conformal procedures around six
distributional discrepancies to give a pre-registered, finite-sample
coverage evaluation of single-cell perturbation predictors.

Output is vector PDF written straight into paper_ai4science/ so the LaTeX
\\includegraphics{fig_overview.pdf} resolves from the compile directory.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "paper_ai4science" / "fig_overview.pdf"
OUT_PAPER = ROOT / "paper" / "fig_overview.pdf"

# Palette (colorblind-safe, muted; matches the green/orange/red family used in
# the K1/K2 figures so the figure set reads as one paper).
C_INPUT = "#dbe7f3"   # light blue
C_DISC = "#e7f0d9"    # light green
C_HEAD = "#fbe6cf"    # light orange
C_OUT = "#ecddef"     # light purple
C_EDGE = "#34495e"
C_BANNER = "#2c3e50"


def _box(ax, x, y, w, h, text, fc, fontsize=9, weight="normal", ec=C_EDGE,
         lw=1.0, ha="center", text_color="black"):
    ax.add_patch(FancyBboxPatch(
        (x, y), w, h, boxstyle="round,pad=0.15,rounding_size=0.6",
        linewidth=lw, edgecolor=ec, facecolor=fc, zorder=2))
    tx = x + w / 2 if ha == "center" else x + 0.6
    ax.text(tx, y + h / 2, text, ha=ha, va="center",
            fontsize=fontsize, weight=weight, color=text_color, zorder=3)


def _stage_title(ax, x, y, text, color):
    ax.text(x, y, text, ha="center", va="center", fontsize=11.5,
            weight="bold", color=color, zorder=3)


def _arrow(ax, x0, y0, x1, y1, color=C_EDGE, lw=1.6):
    ax.add_patch(FancyArrowPatch(
        (x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=14,
        linewidth=lw, color=color, zorder=1,
        shrinkA=2, shrinkB=2))


def main() -> int:
    fig, ax = plt.subplots(figsize=(14, 6.4))
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.axis("off")

    # ---- Pre-registration banner (top, spans the whole pipeline) ----------
    ax.add_patch(FancyBboxPatch(
        (2, 90.5), 96, 7.5, boxstyle="round,pad=0.2,rounding_size=0.8",
        linewidth=1.2, edgecolor=C_BANNER, facecolor="#f4f6f7", zorder=2))
    ax.text(50, 94.2,
            "Pre-registered before any first run (git-stamped c7046e4b, 2026-05-03): "
            "hypotheses, permutation-test thresholds, and BH-FDR corrections committed in advance",
            ha="center", va="center", fontsize=10.2, weight="bold", color=C_BANNER,
            zorder=3)

    # Column centers for the four stages.
    cx = [13, 38, 63, 87.5]
    stage_y = 85.5
    _stage_title(ax, cx[0], stage_y, "Predictors x datasets x splits", "#2166ac")
    _stage_title(ax, cx[1], stage_y, "Six distributional\ndiscrepancies", "#4d7c0f")
    _stage_title(ax, cx[2], stage_y, "Four conformal\nhead levels", "#b15928")
    _stage_title(ax, cx[3], stage_y, "Three pre-registered\nknockout claims", "#762a83")

    # ---- Stage A: inputs --------------------------------------------------
    _box(ax, 2, 60, 22, 18,
         "8 predictors, ~0 to 6x10^8 params\n"
         "Mean - Bilinear ridge - scGen - CPA\n"
         "sVAE+ - biolord - GEARS - STATE-600M",
         C_INPUT, fontsize=8.6, weight="bold")
    _box(ax, 2, 40, 22, 16,
         "5 datasets\n"
         "Norman - Replogle K562 / RPE1\n"
         "Adamson - Tahoe-100M",
         C_INPUT, fontsize=8.6)
    _box(ax, 2, 23, 22, 13,
         "3 split types\n"
         "within-pert - held-out-pert\n"
         "cross-cell-line (K562 -> RPE1)",
         C_INPUT, fontsize=8.6)
    # what flows out of stage A
    ax.text(13, 17.5,
            r"each cell: predicted population $X_{\mathrm{pred}}$ vs observed $X_{\mathrm{obs}}$",
            ha="center", va="center", fontsize=8.2, style="italic", color="#2166ac")

    # ---- Stage B: six discrepancies --------------------------------------
    disc = ["KS (per gene)", "Wasserstein-1", "Energy distance",
            "MMD-RBF", "Bimodality match", "Variance ratio"]
    by0, bh, bgap = 24, 8.2, 1.8
    for i, d in enumerate(reversed(disc)):
        _box(ax, 29, by0 + i * (bh + bgap), 18, bh, d, C_DISC, fontsize=8.8)

    # ---- Stage C: four heads ---------------------------------------------
    heads = [
        ("Per-gene", "CD-split / HPD-split + CQR"),
        ("Per-perturbation", "split-conformal on scores"),
        ("Per-population", "OT-CP Brenier-merge"),
        ("Subgroup-conditional", "Mondrian per-subgroup"),
    ]
    hy0, hh, hgap = 25, 11.0, 2.2
    for i, (name, sub) in enumerate(reversed(heads)):
        _box(ax, 54, hy0 + i * (hh + hgap), 18, hh,
             f"{name}\n({sub})", C_HEAD, fontsize=8.4,
             weight="bold")
    ax.text(63, 19.5,
            r"finite-sample threshold $\hat\tau$, coverage $1-\alpha$",
            ha="center", va="center", fontsize=8.4, style="italic", color="#b15928")

    # ---- Stage D: K1/K2/K3 ------------------------------------------------
    _box(ax, 78, 60, 19.5, 16,
         "K1  Calibration benchmark\n"
         "8 predictors x 5 datasets\n"
         "sVAE+ exact 0.800 nominal\n"
         "on bimodality + variance",
         C_OUT, fontsize=8.2, weight="bold")
    _box(ax, 78, 41, 19.5, 16,
         "K2  Capacity hypothesis\n"
         "pre-registered Spearman test\n"
         "honest null: 2 of 5 pass\n"
         "(RPE1, Tahoe); H1 fails overall",
         C_OUT, fontsize=8.2, weight="bold")
    _box(ax, 78, 23, 19.5, 15,
         "K3  PRISM drug resistance\n"
         "BH-FDR over compound x pathway\n"
         "445 calibrated vs 363 uncal.\n"
         "robust across signature sizes",
         C_OUT, fontsize=8.2, weight="bold")

    # ---- Inter-stage arrows ----------------------------------------------
    _arrow(ax, 24.3, 50, 28.7, 50)          # A -> B
    _arrow(ax, 47.3, 50, 53.7, 50)          # B -> C
    _arrow(ax, 72.3, 50, 77.7, 50)          # C -> D

    # release footer
    ax.text(50, 6.5,
            "Released as a pip-installable confpert library and a Cell-Eval plug-in "
            "(5 of 6 discrepancies absent from the Arc Institute registry)",
            ha="center", va="center", fontsize=8.8, color=C_BANNER)

    fig.tight_layout(pad=0.4)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT_PAPER.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, bbox_inches="tight")
    fig.savefig(OUT_PAPER, bbox_inches="tight")
    plt.close(fig)
    print(f"[fig] -> {OUT}")
    print(f"[fig] -> {OUT_PAPER}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
