"""K2 analysis: pre-registered Spearman tests on H1 (capacity) and H1b (data).

Per preregistration.md (commit c7046e4b):

  H1 (capacity hypothesis): calibration error increases with parameter count.
    Spearman rho between log(params) and |1 - alpha - achieved_coverage| > 0.5
    with p < 0.05 (permutation test) on at least 3 of 4 datasets.

  H1b (data hypothesis): calibration error decreases with training-set-size.
    Spearman rho between log(N_train) and |dev| < -0.5 with p < 0.05 on at
    least 3 of 4 datasets.

This script reads `baselines/results.json` (locked at Phase 1 lock-in) and
computes both Spearman rhos with permutation-test p-values. Output to
`results/k2_preliminary.json` (or `_final.json` once full sweep completes).
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy import stats


# Pre-registered predictor parameter counts and approximate training-set sizes.
# These are the K2 instrument; lock them with the preregistration.
PREDICTOR_PARAMS = {
    "mean": 0,
    "additive": 1,
    "noisy_mean": 1,
    "ahlmann_bilinear_ridge": int(1e4),
    "scgen": int(1e6),
    "cpa": int(5e6),
    "biolord": int(2e7),
    "svaeplus": int(1e7),
    "gears_uncertainty": int(5e7),
    "state": int(6e8),
}

# Approximate training cells per perturbation per dataset on average.
# Real numbers will be filled in from the dataset loaders at K1 lock-in time.
TRAIN_CELLS_BY_DATASET = {
    "norman": 800,
    "replogle_k562": 250,
    "replogle_rpe1": 150,
    "adamson": 600,
    "tahoe": 400,  # subset cells_per_condition; pre-reg label "Tahoe-100M cross-cell-line subset"
    "tahoe_subset": 400,  # legacy alias if any results.json rows use this name
}


def permutation_test_spearman(x: np.ndarray, y: np.ndarray, alternative: str = "greater",
                               n_perm: int = 10_000, seed: int = 0) -> tuple[float, float]:
    """Permutation test for Spearman correlation. Returns (rho, p_value).

    alternative:
      'greater' for H1 (rho > 0)
      'less' for H1b (rho < 0)
      'two-sided' for omnibus
    """
    rho, _ = stats.spearmanr(x, y)
    rng = np.random.RandomState(seed)
    n_extreme = 0
    for _ in range(n_perm):
        y_perm = rng.permutation(y)
        rho_perm, _ = stats.spearmanr(x, y_perm)
        if alternative == "greater":
            if rho_perm >= rho:
                n_extreme += 1
        elif alternative == "less":
            if rho_perm <= rho:
                n_extreme += 1
        else:
            if abs(rho_perm) >= abs(rho):
                n_extreme += 1
    p = (n_extreme + 1) / (n_perm + 1)
    return float(rho), float(p)


def load_calibration_table(results_json: Path,
                            include_split_types: tuple[str, ...] | None = None
                            ) -> list[dict]:
    """Flatten baselines/results.json rows into a per-(predictor, dataset, alpha,
    variant, score) calibration-deviation table.

    By default excludes ``cross_cell_line`` rows (which are split-type 3 and
    not a 5th K2 dataset). Pass ``include_split_types=()`` to include all.
    """
    if include_split_types is None:
        # K2 default: only within-perturbation cells (split #1)
        include_split_types = ("within_perturbation", "")  # legacy rows have no split_type
    with open(results_json) as f:
        data = json.load(f)
    rows = []
    for r in data.get("rows", []):
        if "error" in r and "scores" not in r:
            continue
        if r.get("split_type", "") not in include_split_types:
            continue
        for score_name, sr in r.get("scores", {}).items():
            if isinstance(sr, dict) and "error" not in sr:
                rows.append({
                    "predictor": r["predictor"],
                    "dataset": r["dataset"],
                    "alpha": r["alpha"],
                    "noise_variant": r.get("noise_variant", ""),
                    "score": score_name,
                    "achieved_coverage": sr.get("achieved_coverage"),
                    "calibration_deviation": sr.get("calibration_deviation"),
                })
    return rows


def aggregate_per_predictor_per_dataset(rows: list[dict]) -> dict:
    """For each (predictor, dataset), average calibration deviation over alpha and
    variant and score. Returns {(predictor, dataset): mean_deviation}.
    """
    by_key: dict[tuple, list[float]] = {}
    for r in rows:
        if r["calibration_deviation"] is None:
            continue
        key = (r["predictor"], r["dataset"])
        by_key.setdefault(key, []).append(r["calibration_deviation"])
    return {k: float(np.mean(v)) for k, v in by_key.items() if v}


def k2_test(rows: list[dict], hypothesis: str = "H1",
            n_perm: int = 10_000, seed: int = 0) -> dict:
    """Run pre-registered K2 test on the calibration table.

    hypothesis: 'H1' (capacity, expects rho > 0) or 'H1b' (data, expects rho < 0).
    """
    agg = aggregate_per_predictor_per_dataset(rows)
    if not agg:
        return {"hypothesis": hypothesis, "status": "no_data", "datasets_tested": 0}

    # Group by dataset
    datasets = sorted(set(d for (_, d) in agg))
    per_dataset = {}
    for ds in datasets:
        preds = [(p, dev) for (p, d), dev in agg.items() if d == ds]
        if len(preds) < 3:
            per_dataset[ds] = {"status": "too_few_predictors", "n": len(preds)}
            continue
        x = []
        y = []
        for (pred_name, dev) in preds:
            if hypothesis == "H1":
                params = PREDICTOR_PARAMS.get(pred_name)
                if params is None or params <= 0:
                    # log(0) undefined; use log1p trick
                    x_val = np.log1p(max(params or 0, 0))
                else:
                    x_val = float(np.log(params))
            else:  # H1b
                n_train = TRAIN_CELLS_BY_DATASET.get(ds)
                if n_train is None:
                    continue
                x_val = float(np.log(n_train))
                # H1b uses one X per dataset (constant); the variance is across
                # predictors so the test on a single dataset is not informative.
                # Better: run H1b across DATASETS (one per dataset) with
                # average-deviation-across-predictors as Y. Done in cross_dataset_h1b().
            x.append(x_val)
            y.append(dev)
        x = np.array(x)
        y = np.array(y)
        if len(x) < 3:
            per_dataset[ds] = {"status": "too_few_predictors", "n": len(x)}
            continue
        alternative = "greater" if hypothesis == "H1" else "less"
        rho, p = permutation_test_spearman(x, y, alternative=alternative,
                                            n_perm=n_perm, seed=seed)
        threshold = (rho > 0.5 if hypothesis == "H1" else rho < -0.5) and p < 0.05
        per_dataset[ds] = {
            "rho": rho,
            "p_value": p,
            "n_predictors": int(len(x)),
            "passes_preregistered_threshold": bool(threshold),
        }

    n_datasets_passing = sum(1 for v in per_dataset.values()
                              if isinstance(v, dict) and v.get("passes_preregistered_threshold"))
    return {
        "hypothesis": hypothesis,
        "alternative": "greater" if hypothesis == "H1" else "less",
        "n_datasets_tested": len(per_dataset),
        "n_datasets_passing_threshold": n_datasets_passing,
        "preregistered_success_at_least_3_of_4": n_datasets_passing >= 3,
        "per_dataset": per_dataset,
    }


def cross_dataset_h1b(rows: list[dict],
                      n_perm: int = 10_000, seed: int = 0) -> dict:
    """H1b is more naturally tested across datasets (one X per dataset, Y = mean
    deviation across predictors)."""
    agg = aggregate_per_predictor_per_dataset(rows)
    by_dataset = {}
    for (p, d), dev in agg.items():
        by_dataset.setdefault(d, []).append(dev)
    if len(by_dataset) < 3:
        return {"status": "too_few_datasets"}
    x = []
    y = []
    for ds, devs in by_dataset.items():
        n_train = TRAIN_CELLS_BY_DATASET.get(ds)
        if n_train is None:
            continue
        x.append(np.log(n_train))
        y.append(float(np.mean(devs)))
    x = np.array(x)
    y = np.array(y)
    if len(x) < 3:
        return {"status": "too_few_datasets"}
    rho, p = permutation_test_spearman(x, y, alternative="less",
                                        n_perm=n_perm, seed=seed)
    return {
        "hypothesis": "H1b_cross_dataset",
        "rho": rho,
        "p_value": p,
        "n_datasets": int(len(x)),
        "passes_preregistered_threshold": bool(rho < -0.5 and p < 0.05),
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--results", default="baselines/results.json")
    p.add_argument("--out", default="results/k2_preliminary.json")
    p.add_argument("--n_perm", type=int, default=10_000)
    p.add_argument("--seed", type=int, default=0)
    args = p.parse_args()

    rows = load_calibration_table(Path(args.results))
    print(f"[k2] loaded {len(rows)} rows from {args.results}")

    # Count predictors and datasets present
    preds = sorted(set(r["predictor"] for r in rows))
    datasets = sorted(set(r["dataset"] for r in rows))
    print(f"[k2] predictors present: {preds}")
    print(f"[k2] datasets present: {datasets}")

    h1 = k2_test(rows, hypothesis="H1", n_perm=args.n_perm, seed=args.seed)
    h1b_per = k2_test(rows, hypothesis="H1b", n_perm=args.n_perm, seed=args.seed)
    h1b_cross = cross_dataset_h1b(rows, n_perm=args.n_perm, seed=args.seed)

    out = {
        "n_rows": len(rows),
        "predictors_present": preds,
        "datasets_present": datasets,
        "H1_capacity_hypothesis": h1,
        "H1b_data_hypothesis_per_dataset": h1b_per,
        "H1b_data_hypothesis_cross_dataset": h1b_cross,
        "preregistration_commit": "c7046e4b4c08106acd52c37e7ea8410e126d2ae4",
    }
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2, default=float)
    print(f"[k2] -> {args.out}")
    print(json.dumps({"H1_n_passing": h1.get("n_datasets_passing_threshold"),
                      "H1_overall_pass": h1.get("preregistered_success_at_least_3_of_4"),
                      "H1b_cross_pass": h1b_cross.get("passes_preregistered_threshold")},
                     indent=2))


if __name__ == "__main__":
    main()
