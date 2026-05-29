"""K1 baseline runner.

Per preregistration.md (commit c7046e4b):
  - Fixed dataset, predictor, split-type, alpha grid.
  - Six discrepancies (KS, W1, energy, MMD-RBF, bimodality match, variance ratio).
  - For point-estimate predictors (Mean, Bilinear ridge), sweep four noise variants.
  - Lock baselines/results.json with SHA-256 hash.

Usage:
    python scripts/run_k1_baseline.py --dataset norman --predictor mean
    python scripts/run_k1_baseline.py --dataset norman --predictor ahlmann_bilinear_ridge
    python scripts/run_k1_baseline.py --dataset norman --predictor noisy_mean

Output:
    Appends one row per (predictor, dataset, split, alpha, noise_variant) to
    baselines/results.json. SHA-256 of results.json computed and recorded.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

import numpy as np  # noqa: E402

from confpert.conformal import (  # noqa: E402
    EnergyConformal,
    PerturbationConformal,
)
from confpert.data import (  # noqa: E402
    load_adamson,
    load_norman,
    load_replogle,
    load_replogle_rpe1,
    load_tahoe,
)
# Phase 2 dataset loaders (session 6 additions, see CHANGELOG)
from confpert.data.frangieh import load_frangieh  # noqa: E402
from confpert.data.datlinger import load_datlinger  # noqa: E402
from confpert.data.schmidt import load_schmidt  # noqa: E402
from confpert.data.mcfaline_figueroa import load_mcfaline_figueroa  # noqa: E402
from confpert.data.lara_astiaso import load_lara_astiaso  # noqa: E402
from confpert.metrics import SCORES  # noqa: E402
from confpert.predictors import (  # noqa: E402
    AhlmannEltzeBilinearRidge,
    MeanPredictor,
    NoisyMeanPredictor,
)


DATASET_LOADERS = {
    "norman": load_norman,
    "replogle_k562": load_replogle,
    "replogle_rpe1": load_replogle_rpe1,
    "adamson": load_adamson,
    "tahoe": load_tahoe,
    "frangieh": load_frangieh,
    "datlinger": load_datlinger,
    "schmidt": load_schmidt,
    "mcfaline_figueroa": load_mcfaline_figueroa,
    "lara_astiaso": load_lara_astiaso,
}

# Predictors that need a noise-variant sweep (point-estimate baselines).
NEEDS_NOISE_SWEEP = {"mean", "ahlmann_bilinear_ridge"}

NOISE_VARIANTS = ["A_no_noise", "B_isotropic", "C_per_gene_marginal", "D_full_covariance"]

ALPHAS = [0.05, 0.10, 0.20]
SCORE_NAMES = list(SCORES.keys())  # 6 scores


def build_within_pert_split(ds, train_frac: float = 0.5, calib_frac: float = 0.25,
                            seed: int = 42):
    """50% train, 25% calib, 25% test per perturbation. Returns dicts keyed by pert."""
    rng = np.random.RandomState(seed)
    n_ctrl = ds.X_ctrl.shape[0]
    perm = rng.permutation(n_ctrl)
    n_t = int(n_ctrl * train_frac)
    X_ctrl_train = ds.X_ctrl[perm[:n_t]]
    X_ctrl_other = ds.X_ctrl[perm[n_t:]]

    X_pert_train_list, perts_train_list = [], []
    X_pert_calib, X_pert_test = {}, {}
    for p, X_p in ds.X_pert.items():
        n_p = X_p.shape[0]
        perm_p = rng.permutation(n_p)
        n_t_p = int(n_p * train_frac)
        n_c_p = int(n_p * calib_frac)
        X_pert_train_list.append(X_p[perm_p[:n_t_p]])
        perts_train_list.append(np.array([p] * n_t_p))
        X_pert_calib[p] = X_p[perm_p[n_t_p:n_t_p + n_c_p]]
        X_pert_test[p] = X_p[perm_p[n_t_p + n_c_p:]]
    X_pert_train = np.concatenate(X_pert_train_list).astype(np.float32)
    perts_train = np.concatenate(perts_train_list)
    return X_ctrl_train, X_ctrl_other, X_pert_train, perts_train, X_pert_calib, X_pert_test


def fit_predictor(name: str, X_ctrl_train, perts_train, X_pert_train,
                  noise_variant: str = "A_no_noise", seed: int = 42):
    if name == "mean":
        return MeanPredictor(noise_variant=noise_variant, seed=seed).fit(
            X_ctrl_train, perts_train, X_pert_train
        )
    if name == "ahlmann_bilinear_ridge":
        return AhlmannEltzeBilinearRidge(K=10, lam=0.1, noise_variant=noise_variant,
                                         seed=seed).fit(X_ctrl_train, perts_train, X_pert_train)
    if name == "noisy_mean":
        return NoisyMeanPredictor(seed=seed).fit(X_ctrl_train, perts_train, X_pert_train)
    raise KeyError(f"Unknown predictor {name}")


def run_one_cell(predictor_name: str, dataset_name: str, alpha: float,
                 noise_variant: str, seed: int = 42, n_predict_cells: int = 200):
    """Single (predictor, dataset, alpha, noise_variant) calibration cell.

    Returns a dict suitable for one row of baselines/results.json.
    """
    t0 = time.time()
    loader = DATASET_LOADERS[dataset_name]
    h5ad_name_map = {
        "norman": "norman_2019.h5ad",
        "replogle_k562": "replogle_k562_essential.h5ad",
        "replogle_rpe1": "replogle_rpe1.h5ad",
        "adamson": "adamson_2016.h5ad",
        "tahoe": "tahoe_subset.h5ad",
        "frangieh": "frangieh.h5ad",
        "datlinger": "datlinger.h5ad",
        "schmidt": "schmidt.h5ad",
        "mcfaline_figueroa": "mcfaline_figueroa.h5ad",
        "lara_astiaso": "lara_astiaso_exvivo.h5ad",
    }
    if dataset_name not in h5ad_name_map:
        raise KeyError(f"No h5ad-name mapping for {dataset_name}")
    h5ad = ROOT / "data" / h5ad_name_map[dataset_name]
    print(f"[k1] loading {dataset_name} ...", flush=True)
    # Replogle K562 / RPE1 OOM-killed at full-load on tight-RAM hosts; the
    # chunked path uses backed-mode obs filtering before materializing rows.
    if dataset_name in {"replogle_k562", "replogle_rpe1"}:
        ds = loader(h5ad_path=str(h5ad), chunked=True)
    else:
        ds = loader(h5ad_path=str(h5ad))
    print(f"[k1] loaded ds: n_ctrl={ds.X_ctrl.shape[0]}, n_perts={len(ds.X_pert)}, "
          f"d_genes={ds.X_ctrl.shape[1]}", flush=True)

    X_ctrl_train, X_ctrl_other, X_pert_train, perts_train, X_pert_calib, X_pert_test = \
        build_within_pert_split(ds, seed=seed)

    # Pick n_predict input control cells from X_ctrl_other
    rng = np.random.RandomState(seed)
    if X_ctrl_other.shape[0] > n_predict_cells:
        idx = rng.choice(X_ctrl_other.shape[0], n_predict_cells, replace=False)
        X_input = X_ctrl_other[idx]
    else:
        X_input = X_ctrl_other

    # Fit predictor
    pred = fit_predictor(predictor_name, X_ctrl_train, perts_train, X_pert_train,
                         noise_variant=noise_variant, seed=seed)
    fit_t = time.time() - t0

    # Build calibration and test (pred, obs) pairs per perturbation
    calib_perts = sorted(X_pert_calib.keys())
    rng2 = np.random.RandomState(seed + 1)
    rng2.shuffle(calib_perts)
    n_split = len(calib_perts) // 2
    calib_set = calib_perts[:n_split]
    test_set = calib_perts[n_split:]

    pred_calib_pops = [pred.predict_samples(X_input, p, n_cells=X_pert_calib[p].shape[0])
                       for p in calib_set]
    obs_calib_pops = [X_pert_calib[p] for p in calib_set]

    pred_test_pops = [pred.predict_samples(X_input, p, n_cells=X_pert_test[p].shape[0])
                      for p in test_set]
    obs_test_pops = [X_pert_test[p] for p in test_set]

    # Per-perturbation conformal head, one threshold per discrepancy
    score_results = {}
    for score_name in SCORE_NAMES:
        score_fn = SCORES[score_name]
        try:
            pc = PerturbationConformal(score_fn=score_fn, alpha=alpha)
            pc.calibrate(pred_calib_pops, obs_calib_pops)
            cov = pc.coverage(pred_test_pops, obs_test_pops)
            score_results[score_name] = cov
        except Exception as e:
            score_results[score_name] = {"error": f"{type(e).__name__}: {e}"}

    eval_t = time.time() - t0 - fit_t
    return {
        "predictor": predictor_name,
        "dataset": dataset_name,
        "split_type": "within_perturbation",
        "alpha": alpha,
        "noise_variant": noise_variant,
        "seed": seed,
        "n_predict_cells": n_predict_cells,
        "n_calib_perts": len(calib_set),
        "n_test_perts": len(test_set),
        "fit_sec": fit_t,
        "eval_sec": eval_t,
        "scores": score_results,
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
    }


def update_results_json(out_path: Path, row: dict) -> str:
    """Append row to results.json and recompute SHA-256."""
    if out_path.exists():
        with open(out_path) as f:
            results = json.load(f)
    else:
        results = {"rows": [], "sha256": None}
    results["rows"].append(row)
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2, default=float)
    # SHA-256 of the file (not just the rows)
    h = hashlib.sha256()
    h.update(out_path.read_bytes())
    sha = h.hexdigest()
    # Re-write with the SHA at the top so it's always current
    results["sha256"] = sha
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2, default=float)
    return sha


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--predictor", required=True,
                   choices=list({"mean", "ahlmann_bilinear_ridge", "noisy_mean"}))
    p.add_argument("--dataset", required=True, choices=list(DATASET_LOADERS))
    p.add_argument("--alphas", type=float, nargs="+", default=ALPHAS)
    p.add_argument("--noise_variants", nargs="+", default=NOISE_VARIANTS)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--n_predict_cells", type=int, default=200)
    p.add_argument("--out", default=str(ROOT / "baselines" / "results.json"))
    args = p.parse_args()

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    # Sample sweep ranges
    if args.predictor in NEEDS_NOISE_SWEEP:
        variant_iter = args.noise_variants
    else:
        variant_iter = ["A_no_noise"]  # noise_variant column is required but unused

    for alpha in args.alphas:
        for variant in variant_iter:
            print(f"\n[k1] === {args.predictor} | {args.dataset} | alpha={alpha} | "
                  f"variant={variant} ===", flush=True)
            try:
                row = run_one_cell(
                    predictor_name=args.predictor,
                    dataset_name=args.dataset,
                    alpha=alpha,
                    noise_variant=variant,
                    seed=args.seed,
                    n_predict_cells=args.n_predict_cells,
                )
                sha = update_results_json(out_path, row)
                print(f"[k1] sha256 = {sha[:12]}...", flush=True)
                # Print summary
                for sn, sr in row["scores"].items():
                    if "error" in sr:
                        print(f"  {sn}: ERROR {sr['error']}", flush=True)
                    else:
                        print(f"  {sn}: target={sr['target_coverage']:.2f} "
                              f"achieved={sr['achieved_coverage']:.3f} "
                              f"dev={sr['calibration_deviation']:.3f} tau={sr['tau']:.3f}",
                              flush=True)
            except Exception as e:
                err_row = {
                    "predictor": args.predictor,
                    "dataset": args.dataset,
                    "alpha": alpha,
                    "noise_variant": variant,
                    "seed": args.seed,
                    "error": f"{type(e).__name__}: {e}",
                    "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
                }
                update_results_json(out_path, err_row)
                print(f"[k1] FAILED with {type(e).__name__}: {e}", flush=True)


if __name__ == "__main__":
    main()
