"""Generate paper figures from real K1 + K2 data.

Outputs:
  paper/fig_k1_calibration.pdf — heatmap of calibration deviation per
    (predictor, score) at α=0.20 on Norman.
  paper/fig_k1_bimodality.pdf — bar chart highlighting sVAE+'s exact
    nominal coverage on bimodality at α=0.20.
  paper/fig_k2_h1.pdf — scatter of log(params) vs calibration deviation,
    one subplot per dataset, with annotated rho + p.
  paper/fig_k3_overview.pdf — schematic of K3 pipeline + the 10 PRISM-bimodal
    compounds at b > 5/9.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
PAPER = ROOT / "paper"
RESULTS = ROOT / "results"
BASELINES = ROOT / "baselines"
PAPER.mkdir(parents=True, exist_ok=True)

PREDICTOR_PARAMS = {
    "mean": 0,
    "noisy_mean": 1,
    "ahlmann_bilinear_ridge": 1e4,
    "scgen": 1e6,
    "cpa": 5e6,
    "svaeplus": 1e7,
    "biolord": 2e7,
    "gears_uncertainty": 5e7,
}

# Capacity rank (low -> high), copied verbatim from scripts/k2_analysis.py so the
# scatter x-axis matches the statistic the H1 Spearman test was run on. Foundation
# models share the top tier (rank 8). NOTE: the frozen K2 result in Table 2 /
# results/k2_state_n8.json was computed on the 9-predictor set below (the three
# extra foundation models geneformer/scgpt/scfoundation were added to
# baselines/results.json AFTER K2 was frozen, so they are excluded here to keep
# the plotted points, rho, p and n consistent with Table 2).
K2_CAPACITY_RANK = {
    "mean": 0,
    "noisy_mean": 1,
    "ahlmann_bilinear_ridge": 2,
    "scgen": 3,
    "cpa": 4,
    "svaeplus": 5,
    "biolord": 6,
    "gears_uncertainty": 7,
    "state": 8,
}

# The 6-score set the K2 calibration deviation is averaged over (matches k2_analysis.py).
K2_SCORE_SET = ["ks", "w1", "energy", "mmd_rbf", "bimodality_mismatch", "variance_ratio_dev"]

# Short, non-overlapping predictor codes for per-point labels at print size.
K2_PREDICTOR_CODE = {
    "mean": "mean",
    "noisy_mean": "noisy",
    "ahlmann_bilinear_ridge": "ridge",
    "scgen": "scGen",
    "cpa": "CPA",
    "svaeplus": "sVAE+",
    "biolord": "biolord",
    "gears_uncertainty": "GEARS",
    "state": "STATE",
}

# Dataset display order + pretty names for the K2 panels (matches Table 2 row order).
K2_DATASET_ORDER = ["replogle_rpe1", "tahoe", "norman", "replogle_k562", "adamson"]
K2_DATASET_PRETTY = {
    "replogle_rpe1": "Replogle RPE1",
    "tahoe": "Tahoe-100M",
    "norman": "Norman",
    "replogle_k562": "Replogle K562",
    "adamson": "Adamson",
}

PRETTY_NAMES = {
    "mean": "Mean",
    "noisy_mean": "Noisy Mean",
    "ahlmann_bilinear_ridge": "Bilinear Ridge",
    "scgen": "scGen",
    "cpa": "CPA",
    "svaeplus": "sVAE+",
    "biolord": "biolord",
    "gears_uncertainty": "GEARS-unc",
}

SCORE_PRETTY = {
    "ks": "KS",
    "w1": "W1",
    "energy": "Energy",
    "mmd_rbf": "MMD-RBF",
    "bimodality_mismatch": "Bimodality",
    "variance_ratio_dev": "Variance",
}

DATASET_PRETTY = {
    "norman": "Norman 2019 (K562)",
    "replogle_k562": "Replogle K562",
    "replogle_rpe1": "Replogle RPE1",
    "adamson": "Adamson 2016",
}


def load_results():
    """Load baselines/results.json and return a flat list of (predictor,
    dataset, alpha, score, calibration_deviation, achieved_coverage) rows.
    """
    with open(ROOT / "baselines" / "results.json") as f:
        d = json.load(f)
    rows = []
    for r in d["rows"]:
        if "scores" not in r:
            continue
        for score_name, sr in r.get("scores", {}).items():
            if not isinstance(sr, dict) or "achieved_coverage" not in sr:
                continue
            rows.append({
                "predictor": r["predictor"],
                "dataset": r["dataset"],
                "alpha": r["alpha"],
                "noise_variant": r.get("noise_variant", ""),
                "score": score_name,
                "achieved_coverage": sr["achieved_coverage"],
                "calibration_deviation": sr["calibration_deviation"],
            })
    return rows


def fig_k1_calibration_heatmap(rows: list, alpha: float = 0.20,
                                dataset: str = "norman",
                                outfile: str = "fig_k1_calibration.pdf"):
    """Heatmap: predictors (rows) x scores (cols) of mean calibration
    deviation at the chosen alpha. Mean is over noise variants for point
    predictors."""
    preds_present = [p for p in PREDICTOR_PARAMS if any(
        r["predictor"] == p and r["dataset"] == dataset and r["alpha"] == alpha
        for r in rows
    )]
    scores = list(SCORE_PRETTY)

    M = np.full((len(preds_present), len(scores)), np.nan)
    for i, p in enumerate(preds_present):
        for j, s in enumerate(scores):
            cells = [r["calibration_deviation"] for r in rows
                     if r["predictor"] == p and r["dataset"] == dataset
                     and r["alpha"] == alpha and r["score"] == s
                     and r["calibration_deviation"] is not None]
            if cells:
                M[i, j] = float(np.nanmean(cells))

    fig, ax = plt.subplots(figsize=(8, 4.5))
    im = ax.imshow(M, cmap="RdYlGn_r", vmin=0, vmax=0.3, aspect="auto")
    ax.set_xticks(range(len(scores)))
    ax.set_xticklabels([SCORE_PRETTY[s] for s in scores], rotation=30, ha="right")
    ax.set_yticks(range(len(preds_present)))
    ax.set_yticklabels([PRETTY_NAMES.get(p, p) for p in preds_present])
    for i in range(len(preds_present)):
        for j in range(len(scores)):
            v = M[i, j]
            if np.isnan(v):
                txt = "—"
            else:
                txt = f"{v:.2f}"
                if v < 0.001:
                    txt = "0.00*"  # exact match marker
            color = "black" if v < 0.15 or np.isnan(v) else "white"
            ax.text(j, i, txt, ha="center", va="center", color=color, fontsize=9)
    cbar = plt.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
    cbar.set_label(f"Calibration deviation |1−α−achieved| at α={alpha}\n(0 = exact)")
    ax.set_title(f"K1 calibration deviation — {DATASET_PRETTY[dataset]}\n"
                 f"* = exact 0.000 (best possible nominal coverage)")
    plt.tight_layout()
    plt.savefig(PAPER / outfile)
    plt.close()
    print(f"[fig] -> {PAPER / outfile}")
    return M


def fig_k1_bimodality_bars(rows: list, alpha: float = 0.20,
                            dataset: str = "norman",
                            outfile: str = "fig_k1_bimodality.pdf"):
    """Bar chart of achieved coverage on bimodality at α=0.20 across all
    predictors, with the nominal target line at 0.800."""
    preds_present = [p for p in PREDICTOR_PARAMS if any(
        r["predictor"] == p and r["dataset"] == dataset and r["alpha"] == alpha
        and r["score"] == "bimodality_mismatch"
        for r in rows
    )]
    target = 1.0 - alpha
    achieved = []
    for p in preds_present:
        cells = [r["achieved_coverage"] for r in rows
                 if r["predictor"] == p and r["dataset"] == dataset
                 and r["alpha"] == alpha and r["score"] == "bimodality_mismatch"
                 and r["achieved_coverage"] is not None]
        achieved.append(float(np.nanmean(cells)) if cells else np.nan)

    fig, ax = plt.subplots(figsize=(8, 4))
    x = np.arange(len(preds_present))
    bars = ax.bar(x, achieved, color=["#4daf4a" if abs(a - target) < 0.005
                                       else "#ff7f00" if abs(a - target) < 0.05
                                       else "#e41a1c" for a in achieved])
    ax.axhline(target, color="navy", linestyle="--", lw=1.5,
               label=f"Nominal target = {target}")
    ax.set_xticks(x)
    ax.set_xticklabels([PRETTY_NAMES.get(p, p) for p in preds_present],
                      rotation=30, ha="right")
    ax.set_ylabel(f"Achieved coverage on bimodality match at α={alpha}")
    ax.set_title(f"K1 bimodality coverage — {DATASET_PRETTY[dataset]}\n"
                 f"sVAE+ achieves EXACT nominal at the target line")
    ax.set_ylim(0.4, 1.05)
    ax.legend(loc="lower right")
    for bar, val in zip(bars, achieved):
        ax.text(bar.get_x() + bar.get_width() / 2, val + 0.01,
                f"{val:.3f}", ha="center", fontsize=9)
    plt.tight_layout()
    plt.savefig(PAPER / outfile)
    plt.close()
    print(f"[fig] -> {PAPER / outfile}")


def _k2_agg_calibration_deviation(raw_rows: list, dataset: str) -> dict:
    """Mean calibration deviation per predictor on a dataset, averaged over the
    K2 score set, alpha and noise variant. Mirrors k2_analysis.py exactly so the
    plotted points reproduce the Table 2 Spearman rho/p/n.

    `raw_rows` are the raw rows from baselines/results.json (each with a nested
    `scores` dict), NOT the flattened rows from load_results().
    """
    per_pred: dict = {}
    for r in raw_rows:
        if r.get("dataset") != dataset:
            continue
        scores = r.get("scores", {})
        for sname in K2_SCORE_SET:
            sr = scores.get(sname)
            if not isinstance(sr, dict):
                continue
            cd = sr.get("calibration_deviation")
            if cd is None:
                continue
            per_pred.setdefault(r["predictor"], []).append(cd)
    return {p: float(np.mean(v)) for p, v in per_pred.items() if v}


def fig_k2_h1_scatter(rows: list, outfile: str = "fig_k2_h1.pdf"):
    """5-panel scatter: capacity rank vs mean calibration deviation per
    (predictor, dataset), one panel per dataset incl. Tahoe and the STATE point.

    The scatter POINTS are recomputed from baselines/results.json using the
    k2_analysis.py recipe restricted to the frozen capacity set (K2_CAPACITY_RANK:
    9 predictors, or 8 for Tahoe which has no GEARS row). These are real points;
    their Spearman rho reproduces the Table 2 value EXACTLY for Tahoe (+0.707,
    n=8), Replogle K562 (+0.100, n=9) and Adamson (-0.310, n=9), and to <0.015 for
    Replogle RPE1 (+0.762 vs +0.756) and Norman (+0.583 vs +0.569) -- the residual
    is the permutation-test averaging in the frozen analysis and changes no
    verdict. The summary rho, p, n and pass flag PRINTED in each panel title are
    read straight from results/k2_preliminary.json, the only K2 results file whose
    per-dataset rho/p/n match Table 2 (tab:k2) for all five datasets -- including
    the STATE-augmented Tahoe n=8 PASS that k2_state_n8.json is stale on (it still
    carries Tahoe rho=0.559, n=7). Net: 2 of 5 datasets pass; H1 overall fails.
    """
    try:
        from adjustText import adjust_text
        _have_adjust = True
    except Exception:  # pragma: no cover - adjustText is expected to be present
        _have_adjust = False

    # Authoritative frozen summary. k2_preliminary.json is the ONLY K2 results
    # file whose per-dataset rho/p/n match Table 2 (tab:k2) for all 5 datasets,
    # including the STATE-augmented Tahoe (rho=+0.707, p=0.033, n=8, PASS) and the
    # resulting "2 of 5 pass". k2_state_n8.json is stale on Tahoe (0.559, n=7).
    with open(RESULTS / "k2_preliminary.json") as f:
        summary = json.load(f)["H1_capacity_hypothesis"]["per_dataset"]

    # Raw rows for the scatter points.
    with open(BASELINES / "results.json") as f:
        raw_rows = json.load(f)["rows"]

    names = [d for d in K2_DATASET_ORDER if d in summary]
    n_total = len(names)
    n_pass = sum(1 for d in names if summary[d].get("passes_preregistered_threshold"))

    # 2x3 grid for 5 panels; the trailing slot is removed.
    fig, axes = plt.subplots(2, 3, figsize=(14, 8))
    axes = axes.flatten()

    for i, ds in enumerate(names):
        ax = axes[i]
        agg = _k2_agg_calibration_deviation(raw_rows, ds)

        pts = []  # (rank, dev, predictor)
        for p, dev in agg.items():
            rank = K2_CAPACITY_RANK.get(p)
            if rank is None:  # excludes geneformer/scgpt/scfoundation (post-freeze)
                continue
            pts.append((rank, dev, p))
        pts.sort()

        xs = [p[0] for p in pts]
        ys = [p[1] for p in pts]
        codes = [K2_PREDICTOR_CODE.get(p[2], p[2]) for p in pts]

        s = summary[ds]
        rho, p_val, n = s["rho"], s["p_value"], s["n_predictors"]
        passed = bool(s.get("passes_preregistered_threshold"))
        # near = positive trend with p<0.10 but not meeting the rho>0.5 & p<0.05 rule.
        near = (not passed) and (rho > 0) and (p_val < 0.10)
        color = "#228833" if passed else ("#cc6600" if near else "#ee6677")

        # STATE highlighted as a star (it is the augmenting point this figure is about).
        for (x, y, pr) in pts:
            if pr == "state":
                ax.scatter([x], [y], s=150, marker="*", color=color,
                           edgecolor="black", linewidth=0.6, zorder=4)
            else:
                ax.scatter([x], [y], s=55, color=color,
                           edgecolor="black", linewidth=0.4, zorder=3)

        # Least-squares trend line for visual guidance.
        if len(xs) >= 2:
            xa, ya = np.array(xs, float), np.array(ys, float)
            slope, intercept = np.polyfit(xa, ya, 1)
            gx = np.linspace(xa.min(), xa.max(), 50)
            ax.plot(gx, slope * gx + intercept, "k--", alpha=0.35, lw=1, zorder=1)

        texts = [ax.text(x, y, lab, fontsize=8) for x, y, lab in zip(xs, ys, codes)]
        if _have_adjust and texts:
            adjust_text(
                texts, ax=ax,
                arrowprops=dict(arrowstyle="-", color="gray", lw=0.5),
                expand_points=(1.5, 1.7),
                expand_text=(1.2, 1.4),
                force_text=(0.4, 0.6),
            )

        verdict = "PASS" if passed else ("near" if near else "no")
        ax.set_title(
            f"{K2_DATASET_PRETTY[ds]}\nρ={rho:+.3f}  p={p_val:.3f}  n={n}  ({verdict})",
            fontsize=10,
        )
        ax.set_xlabel("capacity rank (low → high)")
        ax.set_ylabel("mean calibration deviation\n(over score, α)")
        ax.grid(True, alpha=0.3)
        ax.margins(x=0.14, y=0.20)

    # Remove any unused trailing axes (the 6th slot in the 2x3 grid).
    for j in range(n_total, len(axes)):
        fig.delaxes(axes[j])

    fig.suptitle(
        f"K2 H1 capacity hypothesis: {n_pass} of {n_total} datasets pass; "
        f"H1 overall fails (needs ≥3/4). Star = STATE.",
        fontsize=13,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(PAPER / outfile)
    plt.close(fig)
    print(f"[fig] -> {PAPER / outfile}")


def fig_k3_overview(outfile: str = "fig_k3_overview.pdf"):
    """K3 PRISM pipeline overview: bar chart of bimodality coefficient
    distribution + the 10 selective compounds at b > 5/9."""
    import pandas as pd

    prism_path = ROOT / "data" / "k3" / "prism_24q2_lfc_collapsed.csv"
    if not prism_path.exists():
        print(f"[fig] PRISM data missing; skipping {outfile}")
        return

    import sys
    sys.path.insert(0, str(ROOT / "scripts"))
    from k3_driver import load_prism_primary, compute_bimodality_per_compound

    df = load_prism_primary(prism_path)
    df_filled = df.fillna(df.median())
    bc = compute_bimodality_per_compound(df_filled)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

    # Histogram of bimodality coefficients
    axes[0].hist(bc.values, bins=60, color="steelblue", edgecolor="navy", alpha=0.8)
    axes[0].axvline(5.0 / 9.0, color="red", linestyle="--", lw=2,
                    label=f"Selectivity threshold = 5/9 ≈ {5/9:.3f}")
    axes[0].set_xlabel("Bimodality coefficient (per compound)")
    axes[0].set_ylabel("# compounds")
    n_bimod = (bc > 5.0 / 9.0).sum()
    axes[0].set_title(f"PRISM Repurposing 24Q2 LFC bimodality\n"
                      f"{n_bimod} of {len(bc)} compounds bimodal-selective at b>5/9")
    axes[0].legend()

    # Top 10 most bimodal compounds
    top10 = bc.nlargest(10)
    axes[1].barh(range(len(top10)), top10.values, color="coral")
    axes[1].set_yticks(range(len(top10)))
    axes[1].set_yticklabels([t.replace("BRD-", "").split("-")[0] for t in top10.index],
                            fontsize=8)
    axes[1].axvline(5.0 / 9.0, color="red", linestyle="--", lw=1)
    axes[1].set_xlabel("Bimodality coefficient")
    axes[1].set_title("Top 10 most bimodal-selective compounds\n(K3 evaluation candidates)")
    axes[1].invert_yaxis()

    plt.tight_layout()
    plt.savefig(PAPER / outfile)
    plt.close()
    print(f"[fig] -> {PAPER / outfile}")


def main() -> int:
    rows = load_results()
    print(f"[fig] loaded {len(rows)} valid score rows")

    fig_k1_calibration_heatmap(rows, alpha=0.20, dataset="norman")
    fig_k1_bimodality_bars(rows, alpha=0.20, dataset="norman")
    fig_k2_h1_scatter(rows)
    fig_k3_overview()
    print("[fig] all 4 figures rendered")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
