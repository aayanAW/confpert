"""Merge a Modal predictor result JSON into baselines/results.json.

Modal predictor outputs (e.g. confpert_gears_norman.json) have a different shape
than the K1 sweep rows. This script flattens the per-(alpha, score) results into
the K1 row format so they can be appended to baselines/results.json with the
correct metadata.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def merge(modal_json: Path, results_json: Path) -> str:
    with open(modal_json) as f:
        modal = json.load(f)
    cfg = modal["config"]
    predictor = cfg["predictor"]
    # Cross-cell-line entrypoints set "dataset" to a pair name like
    # "replogle_k562_to_rpe1" and add split_type="cross_cell_line".
    dataset = cfg["dataset"]
    split_type = cfg.get("split_type", "within_perturbation")
    fit_sec = modal.get("fit_sec", 0.0)
    n_pred = modal.get("n_pred_pops", 0)

    if results_json.exists():
        with open(results_json) as f:
            results = json.load(f)
    else:
        results = {"rows": [], "sha256": None}

    new_rows = 0
    for alpha_key, alpha_data in modal.get("calibration_results", {}).items():
        # alpha_key is like "alpha_0.05"
        alpha = float(alpha_key.split("_", 1)[1])
        scores = {}
        for score_name, sr in alpha_data.items():
            scores[score_name] = sr
        # For cross_cell_line, n_pred is the number of common perts (each
        # used as both calib AND test, with separate ctrl-input arms). For
        # within_perturbation, n_pred is split 50/50 random.
        if split_type == "cross_cell_line":
            n_calib_perts = n_pred
            n_test_perts = n_pred
        else:
            n_calib_perts = n_pred // 2
            n_test_perts = n_pred - n_pred // 2
        row = {
            "predictor": predictor,
            "dataset": dataset,
            "split_type": split_type,
            "alpha": alpha,
            "noise_variant": "native_uncertainty",
            "seed": 42,
            "n_predict_cells": n_pred,
            "n_calib_perts": n_calib_perts,
            "n_test_perts": n_test_perts,
            "fit_sec": fit_sec,
            "eval_sec": 0.0,
            "scores": scores,
            "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
            "source": str(modal_json.name),
        }
        results["rows"].append(row)
        new_rows += 1

    with open(results_json, "w") as f:
        json.dump(results, f, indent=2, default=float)
    h = hashlib.sha256()
    h.update(results_json.read_bytes())
    sha = h.hexdigest()
    results["sha256"] = sha
    with open(results_json, "w") as f:
        json.dump(results, f, indent=2, default=float)
    return sha


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--modal_json", required=True)
    p.add_argument("--results_json", default=str(ROOT / "baselines" / "results.json"))
    args = p.parse_args()
    sha = merge(Path(args.modal_json), Path(args.results_json))
    print(f"merged {args.modal_json}; sha256 = {sha[:16]}...")


if __name__ == "__main__":
    main()
