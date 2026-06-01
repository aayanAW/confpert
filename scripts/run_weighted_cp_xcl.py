"""Weighted-CP cross-cell-line coverage recovery (Track B Task 5 Step 8).

Demonstrates that Tibshirani et al. (2019) weighted split-conformal restores
finite-sample coverage on the Replogle K562 -> RPE1 transfer where the unweighted
split-conformal head (the paper's reported result) collapses toward 0.

Runs on CPU. The two Replogle h5ads (K562 ~6.7 GB, RPE1) live on the Modal volume
`causeflow-artifacts:/data/`; pull them, or (cleaner, avoids local disk) run on a
Modal CPU worker via `scripts/modal_weighted_cp_xcl.py`, which mounts the volume and
calls `run()` here verbatim.

EXECUTED 2026-05-31 on a Modal CPU worker (ahlmann_bilinear_ridge predictor,
C_per_gene_marginal noise, seed 42); result committed at
`results/weighted_cp_xcl.json`. Finding (35 common perts, 92 genes, ESS 25.2/35):
weighted CP recovers cross-cell-line coverage on bimodality (alpha=0.20:
0.057 -> 0.857) but NOT on the four shape metrics KS/W1/energy/MMD (0.000 ->
0.000) -- for those the RPE1 test scores lie entirely above the K562 calibration
support, so reweighting the calibration set cannot cover them (not a weight-
degeneracy / low-ESS failure; the ESS is healthy). Partial recovery; motivates
Mondrian-by-cell-line stratification. Archival-paper material; the accepted poster
is unchanged.

Pipeline:
  1. Build the K562 (calib) -> RPE1 (test) split via the existing runner's
     `build_replogle_cross_split` (gene-intersected, >=35 common perturbations).
  2. For each common perturbation, score = per-population discrepancy between the
     predictor's samples and the observed cells, computed in the calib (K562) and
     test (RPE1) contexts.
  3. Unweighted CP: threshold = split_conformal_quantile(calib_scores, alpha);
     achieved coverage = mean(test_scores <= tau). Reproduces the collapse.
  4. Weighted CP: weight each calibration perturbation by an estimated covariate-
     shift likelihood ratio w_i = dP_test/dP_calib(x_i), where x_i is the
     perturbation's pseudobulk profile and the LR is read off a logistic-regression
     domain classifier trained on K562 (label 0) vs RPE1 (label 1) pseudobulks.
     WeightedPerturbationConformal then takes the weighted (1-alpha) quantile.
  5. Report unweighted vs weighted achieved coverage per (score, alpha).

Honest by construction: the weight estimate is a real domain classifier, not an
oracle; whatever recovery it yields is reported as-is, including failures.

Usage:
    python scripts/run_weighted_cp_xcl.py --predictor ahlmann_bilinear_ridge \
        --noise C_per_gene_marginal --out results/weighted_cp_xcl.json
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

import numpy as np  # noqa: E402

from confpert.conformal import (  # noqa: E402
    PerturbationConformal,
    WeightedPerturbationConformal,
    split_conformal_quantile,
)
from confpert.metrics import SCORES  # noqa: E402

# Reuse the exact, reviewed cross-split logic.
sys.path.insert(0, str(ROOT / "scripts"))
from run_k1_cross_cell_line import build_replogle_cross_split, fit_predictor  # noqa: E402


def pseudobulk(pop: np.ndarray) -> np.ndarray:
    """Mean expression vector over a cell population (the per-perturbation covariate)."""
    return np.asarray(pop, dtype=np.float64).mean(axis=0)


def domain_classifier_weights(calib_covariates: np.ndarray,
                              k562_bulks: np.ndarray,
                              rpe1_bulks: np.ndarray,
                              clip: float = 20.0) -> np.ndarray:
    """Tibshirani 2019 importance weights w_i = dP_test/dP_calib(x_i), estimated by
    a logistic-regression domain classifier on K562 (source, label 0) vs RPE1
    (target, label 1) perturbation pseudobulks.

    For a calibrated classifier p(target|x), the density ratio is
        dP_target/dP_source(x) = [p/(1-p)] * [(1-pi)/pi]
    where pi is the target prior; the constant (1-pi)/pi cancels in the weighted
    quantile normalisation, so we use the odds p/(1-p) = exp(decision_function).
    The log-odds is clipped to +/-clip for numerical stability.
    """
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler
    from sklearn.pipeline import make_pipeline

    X = np.vstack([k562_bulks, rpe1_bulks])
    y = np.concatenate([np.zeros(len(k562_bulks)), np.ones(len(rpe1_bulks))])
    clf = make_pipeline(
        StandardScaler(),
        LogisticRegression(max_iter=2000, C=1.0, class_weight="balanced"),
    )
    clf.fit(X, y)
    logits = clf.decision_function(calib_covariates)
    logits = np.clip(logits, -clip, clip)
    return np.exp(logits)


def run(predictor_name: str, noise_variant: str, seed: int = 42,
        alphas=(0.05, 0.10, 0.20), n_predict_cells: int = 200) -> dict:
    t0 = time.time()
    print(f"[wcp] building Replogle K562->RPE1 split (predictor={predictor_name}, "
          f"noise={noise_variant}) ...", flush=True)
    (X_ctrl_train, X_ctrl_calib, X_ctrl_rpe1, X_pert_train, perts_train,
     X_pert_calib, X_pert_test, train_genes) = build_replogle_cross_split(seed=seed)

    rng = np.random.RandomState(seed)
    X_input_calib = (X_ctrl_calib[rng.choice(X_ctrl_calib.shape[0], n_predict_cells, replace=False)]
                     if X_ctrl_calib.shape[0] > n_predict_cells else X_ctrl_calib)
    X_input_test = (X_ctrl_rpe1[rng.choice(X_ctrl_rpe1.shape[0], n_predict_cells, replace=False)]
                    if X_ctrl_rpe1.shape[0] > n_predict_cells else X_ctrl_rpe1)

    pred = fit_predictor(predictor_name, X_ctrl_train, perts_train, X_pert_train,
                         noise_variant=noise_variant, seed=seed)

    common = sorted(set(X_pert_calib) & set(X_pert_test))
    print(f"[wcp]   {len(common)} common perturbations; {len(train_genes)} genes", flush=True)

    # Predicted + observed populations in both contexts.
    pred_calib = [pred.predict_samples(X_input_calib, p, n_cells=X_pert_calib[p].shape[0]) for p in common]
    obs_calib = [X_pert_calib[p] for p in common]
    pred_test = [pred.predict_samples(X_input_test, p, n_cells=X_pert_test[p].shape[0]) for p in common]
    obs_test = [X_pert_test[p] for p in common]

    # Covariate per common perturbation = its OBSERVED pseudobulk in each context.
    calib_cov = np.vstack([pseudobulk(o) for o in obs_calib])     # K562 context
    # Domain-classifier training pools: ALL perts' observed pseudobulks per dataset.
    k562_bulks = np.vstack([pseudobulk(o) for o in obs_calib])
    rpe1_bulks = np.vstack([pseudobulk(o) for o in obs_test])
    weights = domain_classifier_weights(calib_cov, k562_bulks, rpe1_bulks)
    print(f"[wcp]   importance weights: min={weights.min():.3g} max={weights.max():.3g} "
          f"ess={(weights.sum()**2)/(weights**2).sum():.1f}/{len(weights)}", flush=True)

    out = {"predictor": predictor_name, "noise_variant": noise_variant, "seed": seed,
           "n_common_perts": len(common), "n_genes": len(train_genes),
           "weight_ess": float((weights.sum()**2) / (weights**2).sum()),
           "per_score": {}}

    for sn in SCORES:
        s_calib = np.array([SCORES[sn](pc, oc) for pc, oc in zip(pred_calib, obs_calib)], float)
        s_test = np.array([SCORES[sn](pt, ot) for pt, ot in zip(pred_test, obs_test)], float)
        rec = {}
        for a in alphas:
            tau_unw = split_conformal_quantile(s_calib, a)
            cov_unw = float(np.mean(s_test <= tau_unw))
            wcp = WeightedPerturbationConformal(alpha=a).calibrate(s_calib, weights=weights)
            cov_w = wcp.coverage(s_test)
            rec[f"alpha_{a}"] = {
                "target": round(1 - a, 3),
                "unweighted_coverage": round(cov_unw, 4),
                "weighted_coverage": round(cov_w, 4),
                "unweighted_dev": round(abs((1 - a) - cov_unw), 4),
                "weighted_dev": round(abs((1 - a) - cov_w), 4),
                "tau_unweighted": round(tau_unw, 6),
                "tau_weighted": round(float(wcp.tau_), 6),
            }
        out["per_score"][sn] = rec
        a2 = rec["alpha_0.2"]
        print(f"[wcp]   {sn:22s} a=0.20 target=0.80  unweighted={a2['unweighted_coverage']:.3f}"
              f"  weighted={a2['weighted_coverage']:.3f}", flush=True)

    out["elapsed_sec"] = round(time.time() - t0, 1)
    return out


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--predictor", default="ahlmann_bilinear_ridge",
                   choices=["mean", "ahlmann_bilinear_ridge", "noisy_mean"])
    p.add_argument("--noise", default="C_per_gene_marginal")
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--n_predict_cells", type=int, default=200)
    p.add_argument("--out", default=str(ROOT / "results" / "weighted_cp_xcl.json"))
    args = p.parse_args()

    res = run(args.predictor, args.noise, seed=args.seed, n_predict_cells=args.n_predict_cells)
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(res, indent=2))
    print(f"[wcp] -> {out_path}", flush=True)

    # Summary verdict at the headline alpha.
    print("\n[wcp] SUMMARY at alpha=0.20 (target coverage 0.80):", flush=True)
    for sn, rec in res["per_score"].items():
        a = rec["alpha_0.2"]
        verdict = "RECOVERS" if a["weighted_dev"] < a["unweighted_dev"] - 0.05 else \
                  ("~same" if abs(a["weighted_dev"] - a["unweighted_dev"]) <= 0.05 else "WORSE")
        print(f"  {sn:22s} unw={a['unweighted_coverage']:.3f} (dev {a['unweighted_dev']:.3f})  "
              f"wtd={a['weighted_coverage']:.3f} (dev {a['weighted_dev']:.3f})  {verdict}", flush=True)


if __name__ == "__main__":
    raise SystemExit(main())
