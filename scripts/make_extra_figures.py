"""Extra ConfPert figures for the expanded body (8-page budget).

fig_alpha_sensitivity.pdf -- achieved coverage versus the nominal target
1-alpha, averaged over predictors and discrepancies, one line per dataset.
Empirically checks the coverage Proposition: points on/above the diagonal are
conservative (valid); points below indicate undercoverage. Takeaway: coverage
is conservative on four of five datasets; Adamson does not track alpha.

Data: baselines/results.json (within_perturbation split, the 9 K1 predictors).
Vector PDF written into paper_ai4science/ and paper/.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
SCORES = ["ks", "w1", "energy", "mmd_rbf", "bimodality_mismatch", "variance_ratio_dev"]
PREDS = {"mean", "noisy_mean", "ahlmann_bilinear_ridge", "scgen", "cpa",
         "biolord", "svaeplus", "gears_uncertainty", "state"}
ALPHAS = [0.05, 0.10, 0.20]
DS_ORDER = ["norman", "replogle_k562", "replogle_rpe1", "adamson", "tahoe"]
DS_PRETTY = {"norman": "Norman", "replogle_k562": "Replogle K562",
             "replogle_rpe1": "Replogle RPE1", "adamson": "Adamson",
             "tahoe": "Tahoe-100M"}
DS_COLOR = {"norman": "#1b9e77", "replogle_k562": "#d95f02",
            "replogle_rpe1": "#7570b3", "adamson": "#e7298a",
            "tahoe": "#66a61e"}


def mean_achieved(rows, ds, alpha):
    vals = []
    for r in rows:
        if (r["predictor"] in PREDS and r["dataset"] == ds
                and abs(r["alpha"] - alpha) < 1e-9
                and r.get("split_type") == "within_perturbation"):
            for s in SCORES:
                sr = r.get("scores", {}).get(s)
                if isinstance(sr, dict) and sr.get("achieved_coverage") is not None:
                    vals.append(sr["achieved_coverage"])
    return float(np.mean(vals)) if vals else None


def main() -> int:
    rows = json.load(open(ROOT / "baselines" / "results.json"))["rows"]
    fig, ax = plt.subplots(figsize=(6.4, 5.2))

    # perfect-calibration diagonal and the conservative (valid) region
    lo, hi = 0.78, 1.0
    ax.plot([lo, hi], [lo, hi], "k--", lw=1.2, zorder=1,
            label="perfect calibration (achieved = nominal)")
    ax.fill_between([lo, hi], [lo, hi], hi, color="#2ca02c", alpha=0.06, zorder=0)
    ax.text(0.86, 0.985, "conservative (valid)", fontsize=8, color="#2ca02c",
            ha="center", style="italic")

    targets = [1 - a for a in ALPHAS]
    for ds in DS_ORDER:
        ys = [mean_achieved(rows, ds, a) for a in ALPHAS]
        xs = [t for t, y in zip(targets, ys) if y is not None]
        ys = [y for y in ys if y is not None]
        if not ys:
            continue
        ax.plot(xs, ys, "-o", color=DS_COLOR[ds], lw=1.8, ms=7,
                markeredgecolor="black", markeredgewidth=0.5,
                label=DS_PRETTY[ds], zorder=3)

    ax.set_xlabel(r"nominal coverage $1-\alpha$  ($\alpha \in \{0.20, 0.10, 0.05\}$)")
    ax.set_ylabel("achieved coverage (mean over predictors $\\times$ discrepancies)")
    ax.set_xlim(lo, 0.99)
    ax.set_ylim(lo, 1.0)
    ax.set_xticks(targets)
    ax.set_title("K1 coverage tracks the nominal target across $\\alpha$\n"
                 "(conservative on 4/5 datasets; Adamson does not track $\\alpha$)",
                 fontsize=11)
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=3,
              fontsize=8.5, framealpha=0.95)

    fig.tight_layout()
    for out in (ROOT / "paper_ai4science" / "fig_alpha_sensitivity.pdf",
                ROOT / "paper" / "fig_alpha_sensitivity.pdf"):
        out.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out, bbox_inches="tight")
        print(f"[fig] -> {out}")
    plt.close(fig)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
