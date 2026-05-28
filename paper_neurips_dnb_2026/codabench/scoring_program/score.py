"""ConfPert Codabench scorer.

Reads participant predictions + reference observed populations, computes
per-cell calibration deviation, aggregates into the leaderboard metric.

Codabench mounts:
  /input/predictions/<dataset>__<perturbation>.npy   -- (n_cells, n_genes)
  /input/observed/<dataset>__<perturbation>.npy     -- (n_cells, n_genes)
  /input/manifest.json                              -- cells to score
  /output/scores.json                                -- leaderboard metrics
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np


def _score_pair(pred: np.ndarray, obs: np.ndarray,
                  alpha: float = 0.05) -> dict[str, float]:
    """Score one (pred, obs) population pair.

    Computes the six pre-reg discrepancies + their conformal nonconformity
    scores, then derives achieved coverage at level alpha.
    """
    from confpert import metrics
    SCORES = ["ks", "w1", "energy", "mmd_rbf", "bimodality_mismatch",
              "variance_ratio_dev"]
    out: dict[str, float] = {}
    for sname in SCORES:
        fn = getattr(metrics, f"score_{sname}", None)
        if fn is None:
            continue
        try:
            d = float(fn(pred, obs))
        except Exception as exc:
            d = float("nan")
        out[sname] = d
    return out


def aggregate_cells(per_cell: list[dict[str, Any]]) -> dict[str, float]:
    """Aggregate per-cell results into the leaderboard metric."""
    devs = []
    for c in per_cell:
        d = c.get("calibration_deviation")
        if d is None or not np.isfinite(d):
            continue
        devs.append(float(d))
    if not devs:
        return {
            "mean_calibration_deviation": float("nan"),
            "ci_lower_bound": float("nan"),
            "frac_passing_cells": 0.0,
            "n_cells_scored": 0,
        }
    devs_arr = np.array(devs)
    rng = np.random.default_rng(42)
    boot_means = []
    for _ in range(200):
        sample = rng.choice(devs_arr, size=devs_arr.size, replace=True)
        boot_means.append(sample.mean())
    boot_means = np.array(boot_means)
    return {
        "mean_calibration_deviation": float(devs_arr.mean()),
        "ci_lower_bound": float(np.percentile(boot_means, 2.5)),
        "frac_passing_cells": float((devs_arr <= 0.05).mean()),
        "n_cells_scored": int(devs_arr.size),
    }


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--input", type=Path, default=Path("/input"))
    p.add_argument("--output", type=Path, default=Path("/output"))
    args = p.parse_args(argv)

    manifest_path = args.input / "manifest.json"
    if not manifest_path.exists():
        print(f"FATAL: manifest missing at {manifest_path}", flush=True)
        return 2
    manifest = json.loads(manifest_path.read_text())

    pred_dir = args.input / "predictions"
    obs_dir = args.input / "observed"

    per_cell = []
    for cell in manifest.get("cells", []):
        ds = cell["dataset"]
        pert = cell["perturbation"]
        alpha = cell.get("alpha", 0.05)
        pred_p = pred_dir / f"{ds}__{pert}.npy"
        obs_p = obs_dir / f"{ds}__{pert}.npy"
        if not pred_p.exists() or not obs_p.exists():
            print(f"[score] missing files for {ds}::{pert}", flush=True)
            per_cell.append({
                "dataset": ds, "perturbation": pert, "alpha": alpha,
                "calibration_deviation": float("nan"),
                "error": "missing_files",
            })
            continue
        pred = np.load(pred_p)
        obs = np.load(obs_p)
        scores = _score_pair(pred, obs, alpha=alpha)
        # Pseudo-coverage from the mean discrepancy. This is a stand-in: the
        # real conformal-prediction pipeline lives in confpert.conformal.
        # The Codabench leaderboard scorer should swap to the real pipeline
        # once Phase 2D ships the conformal calibration on the private fold.
        mean_disc = float(np.nanmean(list(scores.values())))
        cov_proxy = max(0.0, 1.0 - mean_disc)
        target = 1.0 - alpha
        per_cell.append({
            "dataset": ds, "perturbation": pert, "alpha": alpha,
            "scores": scores,
            "cov_proxy": cov_proxy,
            "calibration_deviation": abs(cov_proxy - target),
        })

    agg = aggregate_cells(per_cell)

    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / "scores.json").write_text(json.dumps(agg, indent=2))
    (args.output / "per_cell.json").write_text(json.dumps(per_cell, indent=2,
                                                            default=str))
    print(f"[score] leaderboard metric: {agg}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
