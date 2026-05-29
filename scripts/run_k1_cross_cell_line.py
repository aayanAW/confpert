"""Cross-cell-line K1 baseline runner (pre-reg L32-35, split #3).

Pre-registration locks three K1 split types per dataset:
  1. Within-perturbation 70/30 cell-level split  (run_k1_baseline.py)
  2. Held-out perturbation 5-fold CV             (effective via conformal calib)
  3. Cross-cell-line: K562 calibration, RPE1 test (this script)

For Replogle: train on K562 essential, calibrate on K562 held-out, test on
RPE1. We require perturbations to overlap between K562 and RPE1 datasets
(both are essential-gene CRISPRi screens, so the intersection is large).

For Tahoe-100M (when the subset lands): train on N-2 cell lines, calibrate
on 1, test on 1 held out.

Usage:
    python scripts/run_k1_cross_cell_line.py --predictor mean
    python scripts/run_k1_cross_cell_line.py --predictor ahlmann_bilinear_ridge
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

from confpert.conformal import PerturbationConformal  # noqa: E402
from confpert.data import (  # noqa: E402
    load_replogle,
    load_replogle_rpe1,
    load_tahoe,
)
from confpert.metrics import SCORES  # noqa: E402
from confpert.predictors import (  # noqa: E402
    AhlmannEltzeBilinearRidge,
    MeanPredictor,
    NoisyMeanPredictor,
)


ALPHAS = [0.05, 0.10, 0.20]
NOISE_VARIANTS = ["A_no_noise", "B_isotropic", "C_per_gene_marginal", "D_full_covariance"]
NEEDS_NOISE_SWEEP = {"mean", "ahlmann_bilinear_ridge"}


def _gene_intersection(ds_train, ds_test):
    """Return the gene intersection (in stable order from ds_train) and the
    column indices to apply to each dataset's matrices.
    """
    g_train = list(ds_train.gene_names)
    g_test_set = set(ds_test.gene_names)
    keep = [i for i, g in enumerate(g_train) if g in g_test_set]
    if len(keep) < 50:
        raise RuntimeError(
            f"Only {len(keep)} genes overlap between train and test; "
            f"cross-cell-line K1 requires >=50 shared genes for stable conformal."
        )
    train_idx = np.asarray(keep, dtype=np.int64)
    train_genes = [g_train[i] for i in keep]
    test_name_to_idx = {g: i for i, g in enumerate(ds_test.gene_names)}
    test_idx = np.asarray([test_name_to_idx[g] for g in train_genes], dtype=np.int64)
    return train_genes, train_idx, test_idx


def build_replogle_cross_split(seed: int = 42, train_frac: float = 0.5,
                                calib_frac: float = 0.25,
                                n_top_per_dataset: int = 200):
    """Train on Replogle K562 perturbations (50% cells); calib on K562 (25%);
    test on Replogle RPE1 perturbations (overlap subset).

    K562 essential and RPE1 essential are distinct gene sets by biology, so
    the perturbation overlap at the default n_top=50 each is only ~3 genes.
    We pull n_top_per_dataset=200 from each to enlarge the intersection;
    for context, the published Replogle 2022 K562/RPE1 essentialome files
    share ~700 perturbations at the genome-scale level.
    """
    print("[k1-xcl] loading Replogle K562 (train + calib) ...", flush=True)
    ds_k562 = load_replogle(
        h5ad_path=str(ROOT / "data" / "replogle_k562_essential.h5ad"),
        n_top_perturbations=n_top_per_dataset,
        chunked=True,
    )
    print(f"[k1-xcl]   K562: n_ctrl={ds_k562.X_ctrl.shape[0]} "
          f"n_perts={len(ds_k562.X_pert)} d_genes={ds_k562.X_ctrl.shape[1]}",
          flush=True)
    print("[k1-xcl] loading Replogle RPE1 (test) ...", flush=True)
    ds_rpe1 = load_replogle_rpe1(
        h5ad_path=str(ROOT / "data" / "replogle_rpe1.h5ad"),
        n_top_perturbations=n_top_per_dataset,
        chunked=True,
    )
    print(f"[k1-xcl]   RPE1: n_ctrl={ds_rpe1.X_ctrl.shape[0]} "
          f"n_perts={len(ds_rpe1.X_pert)} d_genes={ds_rpe1.X_ctrl.shape[1]}",
          flush=True)

    # Align gene space
    train_genes, train_idx, test_idx = _gene_intersection(ds_k562, ds_rpe1)
    n_genes = len(train_genes)
    print(f"[k1-xcl]   gene intersection: {n_genes}", flush=True)

    X_ctrl_k562 = ds_k562.X_ctrl[:, train_idx]
    X_ctrl_rpe1 = ds_rpe1.X_ctrl[:, test_idx]

    # Cell-level train/calib split on K562 controls
    rng = np.random.RandomState(seed)
    n_ctrl = X_ctrl_k562.shape[0]
    perm = rng.permutation(n_ctrl)
    n_train = int(n_ctrl * train_frac)
    X_ctrl_train = X_ctrl_k562[perm[:n_train]]
    X_ctrl_calib = X_ctrl_k562[perm[n_train:]]

    # Perturbation-level: train predictor on (50% cells of) K562 perts
    X_pert_train_list, perts_train_list = [], []
    X_pert_calib: dict[str, np.ndarray] = {}
    common_perts = sorted(set(ds_k562.X_pert.keys()) & set(ds_rpe1.X_pert.keys()))
    print(f"[k1-xcl]   K562 ∩ RPE1 perts: {len(common_perts)} "
          f"(K562={len(ds_k562.X_pert)}, RPE1={len(ds_rpe1.X_pert)})", flush=True)
    if len(common_perts) < 5:
        raise RuntimeError(
            f"Only {len(common_perts)} perturbations overlap K562 and RPE1; "
            "cross-cell-line K1 needs >=5 for conformal calibration."
        )

    # K562 train + calib pool: use ALL K562 perts the predictor will see for
    # training; calibration coverage is computed only on the common subset.
    for p, X_p in ds_k562.X_pert.items():
        X_p = X_p[:, train_idx]
        n_p = X_p.shape[0]
        perm_p = rng.permutation(n_p)
        n_train_p = int(n_p * train_frac)
        n_calib_p = int(n_p * calib_frac)
        X_pert_train_list.append(X_p[perm_p[:n_train_p]])
        perts_train_list.append(np.array([p] * n_train_p))
        if p in common_perts:
            # Held-out K562 cells per overlap pert -> calibration arm
            X_pert_calib[p] = X_p[perm_p[n_train_p:n_train_p + n_calib_p]]

    X_pert_train = np.concatenate(X_pert_train_list).astype(np.float32)
    perts_train = np.concatenate(perts_train_list)

    # Test arm: RPE1 perts in the overlap, full RPE1 cell counts (cross-context)
    X_pert_test: dict[str, np.ndarray] = {}
    for p in common_perts:
        X_p = ds_rpe1.X_pert[p][:, test_idx]
        if X_p.shape[0] >= 30:
            X_pert_test[p] = X_p

    return (X_ctrl_train, X_ctrl_calib, X_ctrl_rpe1,
            X_pert_train, perts_train, X_pert_calib, X_pert_test, train_genes)


def build_tahoe_cross_split(cell_lines_train: list[str] | None = None,
                              cell_lines_calib: list[str] | None = None,
                              cell_lines_test: list[str] | None = None,
                              seed: int = 42, train_frac: float = 0.5):
    """Tahoe-100M cross-cell-line split. Trains on N-2 lines, calibrates on
    one held-out line, tests on another. By default (when args=None) takes
    the cell line list from the subset h5ad's adata.uns["chosen_cell_lines"]
    and partitions 70/15/15.
    """
    print("[k1-xcl] loading Tahoe subset (multi-line) ...", flush=True)
    import scanpy as sc
    h5ad = ROOT / "data" / "tahoe_subset.h5ad"
    if not h5ad.exists():
        raise FileNotFoundError(
            f"Tahoe subset not at {h5ad}; build it via "
            "`modal run --detach scripts/modal_launch.py::tahoe_subset_build`."
        )
    adata = sc.read_h5ad(h5ad)
    if "cell_line_id" not in adata.obs.columns:
        raise RuntimeError(
            f"Tahoe subset has no cell_line_id col; obs={list(adata.obs.columns)}"
        )
    all_lines = sorted(adata.obs["cell_line_id"].astype(str).unique().tolist())
    print(f"[k1-xcl]   Tahoe lines available: {len(all_lines)}: {all_lines}",
          flush=True)
    if cell_lines_train is None:
        rng = np.random.RandomState(seed)
        perm = rng.permutation(len(all_lines))
        n_test = max(1, len(all_lines) // 7)
        n_calib = max(1, len(all_lines) // 7)
        cell_lines_test = [all_lines[i] for i in perm[:n_test]]
        cell_lines_calib = [all_lines[i] for i in perm[n_test:n_test + n_calib]]
        cell_lines_train = [all_lines[i] for i in perm[n_test + n_calib:]]
    print(f"[k1-xcl]   train_lines={cell_lines_train}", flush=True)
    print(f"[k1-xcl]   calib_lines={cell_lines_calib}", flush=True)
    print(f"[k1-xcl]   test_lines={cell_lines_test}", flush=True)

    # Use load_tahoe per-line view to keep preprocessing identical to K1 main
    ds_train = load_tahoe(h5ad_path=str(h5ad), cell_line=cell_lines_train[0])
    train_genes, train_idx, _ = _gene_intersection(ds_train, ds_train)  # identity
    # ... For brevity, we expose only Replogle cross-split here. Tahoe x-cl is
    # exercised at a sweep level via the per-line K1 runner once the subset
    # lands; the line-stratified split lives in the tahoe.py loader's
    # `cell_line=` arg for now.
    return None  # placeholder for upstream call site


def fit_predictor(name: str, X_ctrl_train, perts_train, X_pert_train,
                  noise_variant: str = "A_no_noise", seed: int = 42):
    if name == "mean":
        return MeanPredictor(noise_variant=noise_variant, seed=seed).fit(
            X_ctrl_train, perts_train, X_pert_train)
    if name == "ahlmann_bilinear_ridge":
        return AhlmannEltzeBilinearRidge(K=10, lam=0.1,
                                          noise_variant=noise_variant,
                                          seed=seed).fit(
            X_ctrl_train, perts_train, X_pert_train)
    if name == "noisy_mean":
        return NoisyMeanPredictor(seed=seed).fit(
            X_ctrl_train, perts_train, X_pert_train)
    raise KeyError(f"Unknown predictor {name}")


def run_one_cross_cell_line(predictor_name: str, alpha: float,
                             noise_variant: str, dataset_pair: str = "replogle_k562_to_rpe1",
                             seed: int = 42, n_predict_cells: int = 200):
    t0 = time.time()
    if dataset_pair != "replogle_k562_to_rpe1":
        raise NotImplementedError(
            f"dataset_pair={dataset_pair} not yet supported. "
            "Currently only `replogle_k562_to_rpe1` is implemented; "
            "Tahoe cross-cell-line is exercised via per-line load_tahoe(cell_line=...)."
        )
    (X_ctrl_train, X_ctrl_calib, X_ctrl_rpe1, X_pert_train, perts_train,
     X_pert_calib, X_pert_test, train_genes) = build_replogle_cross_split(seed=seed)

    rng = np.random.RandomState(seed)
    if X_ctrl_calib.shape[0] > n_predict_cells:
        idx = rng.choice(X_ctrl_calib.shape[0], n_predict_cells, replace=False)
        X_input_calib = X_ctrl_calib[idx]
    else:
        X_input_calib = X_ctrl_calib
    if X_ctrl_rpe1.shape[0] > n_predict_cells:
        idx = rng.choice(X_ctrl_rpe1.shape[0], n_predict_cells, replace=False)
        X_input_test = X_ctrl_rpe1[idx]
    else:
        X_input_test = X_ctrl_rpe1

    pred = fit_predictor(predictor_name, X_ctrl_train, perts_train, X_pert_train,
                         noise_variant=noise_variant, seed=seed)
    fit_t = time.time() - t0

    calib_perts = sorted(X_pert_calib.keys())
    test_perts = sorted(X_pert_test.keys())
    common = sorted(set(calib_perts) & set(test_perts))
    if len(common) < 4:
        raise RuntimeError(f"Only {len(common)} common perts for calib+test")
    pred_calib_pops = [pred.predict_samples(X_input_calib, p,
                                              n_cells=X_pert_calib[p].shape[0])
                       for p in common]
    obs_calib_pops = [X_pert_calib[p] for p in common]
    pred_test_pops = [pred.predict_samples(X_input_test, p,
                                              n_cells=X_pert_test[p].shape[0])
                       for p in common]
    obs_test_pops = [X_pert_test[p] for p in common]

    score_results = {}
    for sn in SCORES:
        try:
            pc = PerturbationConformal(score_fn=SCORES[sn], alpha=alpha)
            pc.calibrate(pred_calib_pops, obs_calib_pops)
            cov = pc.coverage(pred_test_pops, obs_test_pops)
            score_results[sn] = cov
        except Exception as e:
            score_results[sn] = {"error": f"{type(e).__name__}: {e}"}

    eval_t = time.time() - t0 - fit_t
    return {
        "predictor": predictor_name,
        "dataset": "replogle_cross_cell_line",
        "split_type": "cross_cell_line",
        "calib_context": "replogle_k562",
        "test_context": "replogle_rpe1",
        "alpha": alpha,
        "noise_variant": noise_variant,
        "seed": seed,
        "n_predict_cells": n_predict_cells,
        "n_common_perts": len(common),
        "n_genes": len(train_genes),
        "fit_sec": fit_t,
        "eval_sec": eval_t,
        "scores": score_results,
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
    }


def update_results_json(out_path: Path, row: dict) -> str:
    if out_path.exists():
        with open(out_path) as f:
            results = json.load(f)
    else:
        results = {"rows": [], "sha256": None}
    results["rows"].append(row)
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2, default=float)
    h = hashlib.sha256()
    h.update(out_path.read_bytes())
    sha = h.hexdigest()
    results["sha256"] = sha
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2, default=float)
    return sha


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--predictor", required=True,
                   choices=list({"mean", "ahlmann_bilinear_ridge", "noisy_mean"}))
    p.add_argument("--alphas", type=float, nargs="+", default=ALPHAS)
    p.add_argument("--noise_variants", nargs="+", default=NOISE_VARIANTS)
    p.add_argument("--dataset_pair", default="replogle_k562_to_rpe1")
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--n_predict_cells", type=int, default=200)
    p.add_argument("--out", default=str(ROOT / "baselines" / "results.json"))
    args = p.parse_args()

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    if args.predictor in NEEDS_NOISE_SWEEP:
        variant_iter = args.noise_variants
    else:
        variant_iter = ["A_no_noise"]

    for alpha in args.alphas:
        for variant in variant_iter:
            print(f"\n[k1-xcl] {args.predictor} | {args.dataset_pair} | "
                  f"alpha={alpha} | variant={variant}", flush=True)
            try:
                row = run_one_cross_cell_line(
                    predictor_name=args.predictor,
                    alpha=alpha,
                    noise_variant=variant,
                    dataset_pair=args.dataset_pair,
                    seed=args.seed,
                    n_predict_cells=args.n_predict_cells,
                )
                sha = update_results_json(out_path, row)
                print(f"[k1-xcl] sha256 = {sha[:12]}...", flush=True)
                for sn, sr in row["scores"].items():
                    if "error" in sr:
                        print(f"  {sn}: ERROR {sr['error']}", flush=True)
                    else:
                        print(f"  {sn}: target={sr['target_coverage']:.2f} "
                              f"achieved={sr['achieved_coverage']:.3f} "
                              f"dev={sr['calibration_deviation']:.3f}",
                              flush=True)
            except Exception as e:
                err_row = {
                    "predictor": args.predictor,
                    "dataset": "replogle_cross_cell_line",
                    "split_type": "cross_cell_line",
                    "alpha": alpha,
                    "noise_variant": variant,
                    "error": f"{type(e).__name__}: {e}",
                    "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
                }
                update_results_json(out_path, err_row)
                print(f"[k1-xcl] FAILED: {type(e).__name__}: {e}", flush=True)


if __name__ == "__main__":
    main()
