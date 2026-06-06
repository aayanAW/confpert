"""Generate the ConfPert overview / architecture figure (fig_overview.pdf).

A pure schematic of the end-to-end pipeline, read left to right:
predictors x datasets x splits -> six distributional discrepancies -> four
off-the-shelf conformal head levels -> three pre-registered claims K1/K2/K3.

Boxes carry only short labels so the figure stays legible at two-column
\\textwidth; all specifics (predictor names, parameter range, dataset names,
head methods, claim numbers) live in the LaTeX caption.

Output is vector PDF written into paper_ai4science/ and paper/.
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

C_INPUT = "#dbe7f3"
C_DISC = "#e7f0d9"
C_HEAD = "#fbe6cf"
C_OUT = "#ecddef"
C_EDGE = "#34495e"
C_BANNER = "#2c3e50"

FS_TITLE = 11
FS_BOX = 9
FS_BANNER = 8.5
FS_NOTE = 8


def _box(ax, x, y, w, h, text, fc, fontsize=FS_BOX, weight="normal"):
    ax.add_patch(FancyBboxPatch(
        (x, y), w, h, boxstyle="round,pad=0.2,rounding_size=0.7",
        linewidth=1.1, edgecolor=C_EDGE, facecolor=fc, zorder=2))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center",
            fontsize=fontsize, weight=weight, color="black", zorder=3)


def _title(ax, x, y, text, color):
    ax.text(x, y, text, ha="center", va="center", fontsize=FS_TITLE,
            weight="bold", color=color, zorder=3)


def _arrow(ax, x0, y0, x1, y1):
    ax.add_patch(FancyArrowPatch(
        (x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=16,
        linewidth=1.8, color=C_EDGE, zorder=1, shrinkA=1, shrinkB=1))


def main() -> int:
    fig, ax = plt.subplots(figsize=(7.4, 4.2))
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.axis("off")

    # four well-separated columns (centres), boxes narrow + short labels
    cx = [12, 37, 62, 87]
    HW = 21  # half-... actually full box width
    def x0(c):  # left edge for a centred box of width HW
        return c - HW / 2

    # ---- banner ----------------------------------------------------------
    ax.add_patch(FancyBboxPatch(
        (1, 90.5), 98, 8, boxstyle="round,pad=0.2,rounding_size=0.8",
        linewidth=1.1, edgecolor=C_BANNER, facecolor="#f4f6f7", zorder=2))
    ax.text(50, 94.5,
            "Pre-registered before any first run "
            "(hypotheses, permutation thresholds, BH-FDR corrections)",
            ha="center", va="center", fontsize=FS_BANNER, weight="bold",
            color=C_BANNER, zorder=3)

    # ---- stage titles ----------------------------------------------------
    ty = 87
    _title(ax, cx[0], ty, "Inputs", "#2166ac")
    _title(ax, cx[1], ty, "Discrepancies", "#4d7c0f")
    _title(ax, cx[2], ty, "Conformal heads", "#b15928")
    _title(ax, cx[3], ty, "Claims", "#762a83")

    # ---- stage A: inputs (short) ----------------------------------------
    _box(ax, x0(cx[0]), 60, HW, 12, "8 predictors", C_INPUT, weight="bold")
    _box(ax, x0(cx[0]), 43, HW, 12, "5 datasets", C_INPUT, weight="bold")
    _box(ax, x0(cx[0]), 26, HW, 12, "3 split types", C_INPUT, weight="bold")
    ax.text(cx[0], 19,
            r"$X_{\mathrm{pred}}$ vs $X_{\mathrm{obs}}$",
            ha="center", va="center", fontsize=FS_NOTE, style="italic",
            color="#2166ac")

    # ---- stage B: six discrepancies -------------------------------------
    disc = ["KS", "Wasserstein-1", "Energy", "MMD-RBF",
            "Bimodality", "Variance ratio"]
    bh, bgap, by0 = 7.9, 1.7, 23
    for i, d in enumerate(reversed(disc)):
        _box(ax, x0(cx[1]), by0 + i * (bh + bgap), HW, bh, d, C_DISC)

    # ---- stage C: four heads (short) ------------------------------------
    heads = ["Per-gene", "Per-perturbation", "Per-population",
             "Subgroup-conditional"]
    hh, hgap, hy0 = 11.5, 2.4, 26
    for i, h in enumerate(reversed(heads)):
        _box(ax, x0(cx[2]), hy0 + i * (hh + hgap), HW, hh, h, C_HEAD,
             weight="bold")
    ax.text(cx[2], 19,
            r"threshold $\hat\tau$, coverage $1-\alpha$",
            ha="center", va="center", fontsize=FS_NOTE, style="italic",
            color="#b15928")

    # ---- stage D: K1/K2/K3 (short) --------------------------------------
    _box(ax, x0(cx[3]), 59, HW, 13,
         "K1\nCalibration\nbenchmark", C_OUT, weight="bold")
    _box(ax, x0(cx[3]), 42, HW, 13,
         "K2\nCapacity\n(honest null)", C_OUT, weight="bold")
    _box(ax, x0(cx[3]), 25, HW, 13,
         "K3\nPRISM\nenrichment", C_OUT, weight="bold")

    # ---- arrows ----------------------------------------------------------
    for a, b in zip(cx[:-1], cx[1:]):
        _arrow(ax, a + HW / 2 + 0.5, 48, b - HW / 2 - 0.5, 48)

    # ---- footer ----------------------------------------------------------
    ax.text(50, 8,
            "pip-installable confpert library + Cell-Eval plug-in",
            ha="center", va="center", fontsize=FS_NOTE, color=C_BANNER)

    fig.tight_layout(pad=0.3)
    for out in (OUT, OUT_PAPER):
        out.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out, bbox_inches="tight")
        print(f"[fig] -> {out}")
    plt.close(fig)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
