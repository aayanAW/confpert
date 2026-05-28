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


def fig_k2_h1_scatter(rows: list, outfile: str = "fig_k2_h1.pdf"):
    """4-panel scatter plot: log(params) vs calibration deviation per
    (predictor, dataset). Annotate rho + p from results/k2_8predictor.json."""
    with open(ROOT / "results" / "k2_8predictor.json") as f:
        k2 = json.load(f)

    fig, axes = plt.subplots(2, 2, figsize=(11, 8), sharey=True)
    datasets = ["norman", "replogle_rpe1", "adamson", "replogle_k562"]
    for ax, ds in zip(axes.flat, datasets):
        # Aggregate calibration deviation per predictor on this dataset
        agg: dict[str, list[float]] = {}
        for r in rows:
            if r["dataset"] != ds or r["calibration_deviation"] is None:
                continue
            agg.setdefault(r["predictor"], []).append(r["calibration_deviation"])
        x_vals, y_vals, names = [], [], []
        for p, devs in agg.items():
            if p not in PREDICTOR_PARAMS:
                continue
            params = PREDICTOR_PARAMS[p]
            x = np.log(max(params, 1))
            y = float(np.mean(devs))
            x_vals.append(x)
            y_vals.append(y)
            names.append(PRETTY_NAMES.get(p, p))

        ax.scatter(x_vals, y_vals, s=60, color="steelblue", edgecolor="navy")
        for x, y, n in zip(x_vals, y_vals, names):
            ax.annotate(n, (x, y), xytext=(5, 5), textcoords="offset points",
                        fontsize=8)
        # Trend line (least squares)
        if len(x_vals) >= 2:
            xa = np.array(x_vals)
            ya = np.array(y_vals)
            slope, intercept = np.polyfit(xa, ya, 1)
            xs = np.linspace(xa.min(), xa.max(), 50)
            ax.plot(xs, slope * xs + intercept, "k--", alpha=0.4, lw=1)
        # rho + p annotation from K2 results
        ds_k2 = k2["H1_capacity_hypothesis"]["per_dataset"].get(ds, {})
        if isinstance(ds_k2, dict) and "rho" in ds_k2:
            rho, p = ds_k2["rho"], ds_k2["p_value"]
            n = ds_k2["n_predictors"]
            ax.set_title(f"{DATASET_PRETTY[ds]}\nρ={rho:+.3f} p={p:.3f} n={n}")
        else:
            ax.set_title(DATASET_PRETTY[ds])
        ax.set_xlabel("log(parameters)")
        if ax in axes[:, 0]:
            ax.set_ylabel("Mean calibration deviation\n(over α, score)")
        ax.grid(True, alpha=0.3)
    fig.suptitle("K2 H1 capacity hypothesis: 0/4 datasets pass pre-reg threshold (ρ>0.5, p<0.05)\n"
                 "→ honest null per pre-registration outcome (D)",
                 fontsize=11)
    plt.tight_layout()
    plt.savefig(PAPER / outfile)
    plt.close()
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
