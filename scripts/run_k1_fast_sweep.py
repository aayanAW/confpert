"""Fast local lightweight sweep: loads each dataset ONCE then iterates
predictors x alphas x noise variants. Bypasses the per-call h5ad load
that makes scripts/run_k1_baseline.py slow on 1.5GB Frangieh.

Writes rows to baselines/results.json using update_results_json from
run_k1_baseline.
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

import numpy as np
import run_k1_baseline as rkb  # noqa: E402


def run_dataset(dataset_name: str, predictors: list[str], alphas: list[float],
                noise_variants: list[str], out_path: Path, seed: int = 42,
                n_predict_cells: int = 200):
    """Load dataset once; iterate predictors x alphas x variants."""
    print(f"\n[fast] === DATASET: {dataset_name} ===", flush=True)
    loader = rkb.DATASET_LOADERS[dataset_name]
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
    h5ad = ROOT / "data" / h5ad_name_map[dataset_name]
    print(f"[fast] loading {dataset_name} from {h5ad} ...", flush=True)
    t0 = time.time()
    if dataset_name in {"replogle_k562", "replogle_rpe1"}:
        ds = loader(h5ad_path=str(h5ad), chunked=True)
    else:
        ds = loader(h5ad_path=str(h5ad))
    print(f"[fast] loaded ds: n_ctrl={ds.X_ctrl.shape[0]} "
          f"n_perts={len(ds.X_pert)} d_genes={ds.X_ctrl.shape[1]} "
          f"({time.time()-t0:.1f}s)", flush=True)

    (X_ctrl_train, X_ctrl_other, X_pert_train, perts_train,
     X_pert_calib, X_pert_test) = rkb.build_within_pert_split(ds, seed=seed)
    rng = np.random.RandomState(seed)
    if X_ctrl_other.shape[0] > n_predict_cells:
        idx = rng.choice(X_ctrl_other.shape[0], n_predict_cells, replace=False)
        X_predict_input = X_ctrl_other[idx]
    else:
        X_predict_input = X_ctrl_other

    for predictor in predictors:
        variants = (noise_variants if predictor in rkb.NEEDS_NOISE_SWEEP
                    else ["A_no_noise"])
        for alpha in alphas:
            for variant in variants:
                print(f"[fast] {predictor} | {dataset_name} | "
                      f"a={alpha} v={variant}", flush=True)
                try:
                    t1 = time.time()
                    pred = rkb.fit_predictor(predictor, X_ctrl_train,
                                              perts_train, X_pert_train,
                                              noise_variant=variant,
                                              seed=seed)
                    # Replicate run_one_cell post-fit logic
                    pred_calib, obs_calib = [], []
                    pred_test, obs_test = [], []
                    for p in X_pert_calib:
                        n_calib = X_pert_calib[p].shape[0]
                        n_test = X_pert_test[p].shape[0]
                        if n_calib < 5 or n_test < 5:
                            continue
                        pc = pred.predict_samples(X_predict_input, p,
                                                   n_cells=n_calib)
                        pt = pred.predict_samples(X_predict_input, p,
                                                   n_cells=n_test)
                        pred_calib.append(pc)
                        obs_calib.append(X_pert_calib[p])
                        pred_test.append(pt)
                        obs_test.append(X_pert_test[p])
                    from confpert.conformal import PerturbationConformal
                    from confpert.metrics import SCORES
                    scores_out = {}
                    for sn, sfn in SCORES.items():
                        try:
                            pc_obj = PerturbationConformal(score_fn=sfn,
                                                            alpha=alpha)
                            pc_obj.calibrate(pred_calib, obs_calib)
                            scores_out[sn] = pc_obj.coverage(pred_test, obs_test)
                        except Exception as e:
                            scores_out[sn] = {"error": f"{type(e).__name__}: {e}"}
                    row = {
                        "predictor": predictor,
                        "dataset": dataset_name,
                        "alpha": alpha,
                        "noise_variant": variant,
                        "seed": seed,
                        "scores": scores_out,
                        "fit_sec": float(time.time() - t1),
                        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
                    }
                    sha = rkb.update_results_json(out_path, row)
                    print(f"[fast]   sha={sha[:12]}  fit={row['fit_sec']:.1f}s",
                          flush=True)
                except Exception as e:
                    err_row = {
                        "predictor": predictor, "dataset": dataset_name,
                        "alpha": alpha, "noise_variant": variant, "seed": seed,
                        "error": f"{type(e).__name__}: {e}",
                        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
                    }
                    rkb.update_results_json(out_path, err_row)
                    print(f"[fast]   FAILED: {type(e).__name__}: {e}",
                          flush=True)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--datasets", nargs="+", required=True)
    p.add_argument("--predictors", nargs="+",
                   default=["mean", "ahlmann_bilinear_ridge", "noisy_mean"])
    p.add_argument("--alphas", nargs="+", type=float, default=rkb.ALPHAS)
    p.add_argument("--noise_variants", nargs="+", default=rkb.NOISE_VARIANTS)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--out", default=str(ROOT / "baselines" / "results.json"))
    args = p.parse_args()

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    for ds in args.datasets:
        run_dataset(ds, args.predictors, args.alphas, args.noise_variants,
                    out_path, seed=args.seed)


if __name__ == "__main__":
    main()
