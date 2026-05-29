"""Phase 2 figure-generation script.

Per pre-reg v2 §5.2 figure plan: 4 main-text figures.

  Fig 1: Benchmark schematic (predictors × datasets × splits × discrepancies
         × heads × α levels). Static + manual layout. Out-of-scope here.

  Fig 2: K1 calibration deviation heatmap (predictor × dataset × discrepancy)
         with bootstrap 95% CI bars. Auto-generated from
         paper_neurips_dnb_2026/phase2d_report_phase{1,2}.json.

  Fig 3: K2 cell-line covariate panel — per-dataset Spearman scatter with
         K562/non-K562 color, expanded to 10 datasets.

  Fig 4: Calibration-vs-accuracy Pareto — scatter (mean L2 error) x
         (mean calibration deviation) for each (predictor, dataset).

Run:
  python scripts/make_phase2_figures.py \
    --phase2d-report paper_neurips_dnb_2026/phase2d_report_phase1.json \
    --out paper_neurips_dnb_2026/figures/

Outputs:
  paper_neurips_dnb_2026/figures/fig2_k1_heatmap.pdf
  paper_neurips_dnb_2026/figures/fig3_k2_covariate.pdf
  paper_neurips_dnb_2026/figures/fig4_pareto.pdf

Each figure is also written as .png for README embedding.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np


def _safe_import_mpl():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        "font.family": "serif",
        "font.size": 9,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 0.5,
    })
    return matplotlib, plt


def make_fig2_heatmap(phase2d_report: dict[str, Any], out_path: Path) -> None:
    """K1 calibration deviation heatmap (predictor x dataset x discrepancy).

    Reads `phase2d_report.bootstrap_ci_table`, computes per-(predictor,
    dataset, score) coverage and calibration deviation (1 - achieved_coverage
    vs 1 - target_coverage), heatmap-grids by score.
    """
    _, plt = _safe_import_mpl()

    bt = phase2d_report["bootstrap_ci_table"]
    predictors = sorted({cell["predictor"] for cell in bt.values()
                          if cell["predictor"]})
    datasets = sorted({cell["dataset"] for cell in bt.values()
                        if cell["dataset"]})
    scores = sorted({cell["score"] for cell in bt.values()
                      if cell["score"]})

    fig, axes = plt.subplots(
        2, 3, figsize=(11, 6.5), sharex=True, sharey=True,
        constrained_layout=True,
    )
    axes = axes.ravel()

    for ax, score in zip(axes, scores):
        grid = np.full((len(predictors), len(datasets)), np.nan)
        for cell in bt.values():
            if cell["score"] != score:
                continue
            i = predictors.index(cell["predictor"])
            j = datasets.index(cell["dataset"])
            # |achieved - target| as calibration deviation. Use mean alpha.
            p = cell["p_hat"]
            target = 1.0 - cell["alpha"]
            grid[i, j] = abs(p - target) if 0.0 <= p <= 1.0 else np.nan
        im = ax.imshow(grid, aspect="auto", cmap="RdYlBu_r",
                        vmin=0.0, vmax=0.20)
        ax.set_title(score, fontsize=9)
        ax.set_xticks(range(len(datasets)))
        ax.set_xticklabels(datasets, rotation=45, ha="right", fontsize=7)
        ax.set_yticks(range(len(predictors)))
        ax.set_yticklabels(predictors, fontsize=7)

    fig.colorbar(im, ax=axes[-1], label="|achieved − target| (deviation)",
                  shrink=0.7)
    fig.suptitle("K1 calibration deviation per (predictor, dataset, discrepancy)",
                  fontsize=10)
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    fig.savefig(out_path.with_suffix(".png"), dpi=200, bbox_inches="tight")
    plt.close(fig)


def make_fig3_k2_covariate(phase2d_report: dict[str, Any], out_path: Path) -> None:
    """K2 cell-line covariate scatter: mean deviation per predictor stratified
    by K562 vs non-K562 dataset context.
    """
    _, plt = _safe_import_mpl()

    bt = phase2d_report["bootstrap_ci_table"]
    # K562 vs non-K562 by dataset name heuristic
    K562_DATASETS = {"norman", "replogle_k562", "adamson", "replogle_rpe1"}
    by_predictor = {}
    for cell in bt.values():
        pred = cell["predictor"]
        ds = cell["dataset"]
        is_k562 = ds in K562_DATASETS
        target = 1.0 - cell["alpha"]
        p = cell["p_hat"]
        if not (0.0 <= p <= 1.0):
            continue
        dev = abs(p - target)
        by_predictor.setdefault(pred, {"k562": [], "non_k562": []})
        bucket = "k562" if is_k562 else "non_k562"
        by_predictor[pred][bucket].append(dev)

    fig, ax = plt.subplots(figsize=(7, 4))
    predictors = sorted(by_predictor.keys())
    for i, pred in enumerate(predictors):
        b = by_predictor[pred]
        if b["k562"]:
            ax.scatter([i - 0.15] * len(b["k562"]), b["k562"],
                        color="#1f77b4", alpha=0.5, s=14, label="K562" if i == 0 else None)
        if b["non_k562"]:
            ax.scatter([i + 0.15] * len(b["non_k562"]), b["non_k562"],
                        color="#d62728", alpha=0.5, s=14,
                        label="non-K562" if i == 0 else None)

    ax.set_xticks(range(len(predictors)))
    ax.set_xticklabels(predictors, rotation=45, ha="right", fontsize=7)
    ax.set_ylabel("|achieved − target| calibration deviation")
    ax.legend(frameon=False, fontsize=8)
    ax.set_title("K2 cell-line context covariate")
    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    fig.savefig(out_path.with_suffix(".png"), dpi=200, bbox_inches="tight")
    plt.close(fig)


def make_fig4_pareto(phase2d_report: dict[str, Any], out_path: Path) -> None:
    """Calibration-vs-accuracy Pareto. Phase 2 has L2 error per cell once
    Phase 2D includes accuracy fields; for Phase 1 we use 'calibration_deviation'
    on x and 'n_perts_test' on y as a placeholder.
    """
    _, plt = _safe_import_mpl()
    bt = phase2d_report["bootstrap_ci_table"]
    fig, ax = plt.subplots(figsize=(6, 4.5))

    predictors_seen: set[str] = set()
    cmap = plt.cm.tab10
    pred_to_color = {}
    for cell in bt.values():
        if cell["predictor"] not in pred_to_color:
            pred_to_color[cell["predictor"]] = cmap(
                len(pred_to_color) % cmap.N)
    for cell in bt.values():
        p = cell["p_hat"]
        if not (0.0 <= p <= 1.0):
            continue
        target = 1.0 - cell["alpha"]
        dev = abs(p - target)
        # Placeholder y: n_perts_test as accuracy-related (Phase 2D will
        # populate with real L2). Replace once Phase 2 results.json includes
        # an `l2_error` field per cell.
        y = cell.get("n_perts_test", 0)
        label = cell["predictor"] if cell["predictor"] not in predictors_seen else None
        predictors_seen.add(cell["predictor"])
        ax.scatter(dev, y, color=pred_to_color[cell["predictor"]], alpha=0.4,
                    s=14, label=label)

    ax.set_xlabel("calibration deviation")
    ax.set_ylabel("n_perts_test (placeholder; Phase 2D will use L2 accuracy)")
    ax.set_title("Calibration-vs-accuracy Pareto (Fig 4, placeholder y-axis)")
    ax.legend(frameon=False, fontsize=7, ncol=2)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    fig.savefig(out_path.with_suffix(".png"), dpi=200, bbox_inches="tight")
    plt.close(fig)


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description="Phase 2 figure generator")
    p.add_argument("--phase2d-report", required=True, type=Path)
    p.add_argument("--out", required=True, type=Path,
                    help="Output figures directory")
    args = p.parse_args(argv)

    args.out.mkdir(parents=True, exist_ok=True)
    report = json.loads(args.phase2d_report.read_text())

    print(f"[fig] generating fig2_k1_heatmap.pdf ...")
    make_fig2_heatmap(report, args.out / "fig2_k1_heatmap.pdf")
    print(f"[fig] generating fig3_k2_covariate.pdf ...")
    make_fig3_k2_covariate(report, args.out / "fig3_k2_covariate.pdf")
    print(f"[fig] generating fig4_pareto.pdf ...")
    make_fig4_pareto(report, args.out / "fig4_pareto.pdf")
    print(f"[fig] done. Figures at {args.out}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
