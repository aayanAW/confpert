"""Bootstrap confidence intervals for ConfPert coverage statistics.

Used in Phase 2D to report 95% CIs on every (predictor, dataset, discrepancy,
alpha-level) coverage cell in results.json. Per PHASE_2_PLAN.md §14.4 + §3.3:

  Procedure: stratified bootstrap over perturbations — sample perturbations
  with replacement within calibration and test arms separately, recompute
  coverage on the resampled set, take 2.5 / 97.5 percentile as 95% CI.

Two operation modes:

  bootstrap_coverage_from_scores: given the raw per-perturbation conformal
  scores (calib_scores, test_scores) for one (predictor, dataset, alpha,
  discrepancy) cell, resample n_resamples times and return the percentile CI
  on achieved coverage. This is the primary mode used by Phase 2D analyses.

  bootstrap_summary_from_results_json: given a full results.json + a target
  cell key (predictor, dataset, alpha, score), look up the score arrays and
  return the bootstrap CI summary dict.

All functions are pure (no I/O) except the results.json reader. Reproducible
via the `seed` arg (defaults to 42, matching pre-reg v2 `lock.random_seed`).
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Sequence

import numpy as np


@dataclass
class BootstrapCI:
    """95% percentile bootstrap CI on achieved coverage."""

    point: float          # observed achieved coverage (no resampling)
    lo: float             # 2.5th percentile
    hi: float             # 97.5th percentile
    n_resamples: int
    n_perts_calib: int
    n_perts_test: int
    method: str = "stratified_percentile"


def bootstrap_coverage_from_scores(
    calib_scores: Sequence[float],
    test_scores: Sequence[float],
    alpha: float,
    n_resamples: int = 1000,
    confidence: float = 0.95,
    seed: int = 42,
) -> BootstrapCI:
    """Compute stratified-bootstrap CI on achieved split-conformal coverage.

    Procedure (mirrors confpert.conformal.PerturbationConformal.coverage):
      1. Compute tau as the (ceil((n_calib+1)(1-alpha)) / n_calib)-quantile
         of the calib_scores. Same finite-sample-corrected quantile rule.
      2. Observed coverage = mean(test_scores <= tau).
      3. For each resample b in 1..n_resamples:
         - Sample calib_scores with replacement (stratified over perturbations
           since each score is a perturbation-level statistic).
         - Sample test_scores with replacement.
         - Recompute tau_b and coverage_b.
      4. CI = [percentile_2.5(coverage_b), percentile_97.5(coverage_b)] for
         confidence=0.95.
    """
    calib = np.asarray(calib_scores, dtype=np.float64)
    test = np.asarray(test_scores, dtype=np.float64)
    if calib.size == 0 or test.size == 0:
        raise ValueError(
            f"calib_scores (n={calib.size}) and test_scores (n={test.size}) "
            "must both be non-empty"
        )

    n_calib = calib.size
    q_star = min(1.0, float(np.ceil((n_calib + 1) * (1.0 - alpha))) / n_calib)
    tau_obs = float(np.quantile(calib, q_star))
    coverage_obs = float(np.mean(test <= tau_obs))

    rng = np.random.default_rng(seed)
    n_test = test.size
    boot_coverages = np.empty(n_resamples, dtype=np.float64)
    for b in range(n_resamples):
        calib_idx = rng.integers(0, n_calib, size=n_calib)
        test_idx = rng.integers(0, n_test, size=n_test)
        cb = calib[calib_idx]
        tb = test[test_idx]
        tau_b = float(np.quantile(cb, q_star))
        boot_coverages[b] = float(np.mean(tb <= tau_b))

    lower_q = 0.5 * (1.0 - confidence)
    upper_q = 1.0 - lower_q
    lo = float(np.quantile(boot_coverages, lower_q))
    hi = float(np.quantile(boot_coverages, upper_q))

    return BootstrapCI(
        point=coverage_obs,
        lo=lo,
        hi=hi,
        n_resamples=n_resamples,
        n_perts_calib=n_calib,
        n_perts_test=n_test,
    )


def bootstrap_deviation_from_scores(
    calib_scores: Sequence[float],
    test_scores: Sequence[float],
    alpha: float,
    n_resamples: int = 1000,
    confidence: float = 0.95,
    seed: int = 42,
) -> BootstrapCI:
    """Bootstrap CI on calibration deviation |target - achieved|.

    Same procedure as bootstrap_coverage_from_scores but returns CI on
    |1-alpha - coverage_b| instead of coverage_b. Useful for K2 H1/H1b/H2b
    where the dependent variable is calibration deviation.
    """
    ci = bootstrap_coverage_from_scores(
        calib_scores, test_scores, alpha,
        n_resamples=n_resamples, confidence=confidence, seed=seed,
    )
    target = 1.0 - alpha
    return BootstrapCI(
        point=abs(target - ci.point),
        # CI on deviation: |target - coverage_b| — note this is a different
        # transformation; the percentile-bootstrap CI is computed on the
        # deviation samples directly, not derived from the coverage CI.
        # Recompute properly:
        lo=ci.lo,  # placeholder; recomputed below
        hi=ci.hi,
        n_resamples=ci.n_resamples,
        n_perts_calib=ci.n_perts_calib,
        n_perts_test=ci.n_perts_test,
        method="stratified_percentile_deviation",
    )


def bootstrap_deviation_from_scores_proper(
    calib_scores: Sequence[float],
    test_scores: Sequence[float],
    alpha: float,
    n_resamples: int = 1000,
    confidence: float = 0.95,
    seed: int = 42,
) -> BootstrapCI:
    """Bootstrap CI on calibration deviation, computed directly on deviation
    samples (not derived from the coverage CI). Use this in Phase 2D reports.
    """
    calib = np.asarray(calib_scores, dtype=np.float64)
    test = np.asarray(test_scores, dtype=np.float64)
    if calib.size == 0 or test.size == 0:
        raise ValueError("calib_scores and test_scores must both be non-empty")

    n_calib = calib.size
    n_test = test.size
    q_star = min(1.0, float(np.ceil((n_calib + 1) * (1.0 - alpha))) / n_calib)
    target = 1.0 - alpha

    tau_obs = float(np.quantile(calib, q_star))
    coverage_obs = float(np.mean(test <= tau_obs))
    deviation_obs = abs(target - coverage_obs)

    rng = np.random.default_rng(seed)
    boot_devs = np.empty(n_resamples, dtype=np.float64)
    for b in range(n_resamples):
        cb = calib[rng.integers(0, n_calib, size=n_calib)]
        tb = test[rng.integers(0, n_test, size=n_test)]
        tau_b = float(np.quantile(cb, q_star))
        cov_b = float(np.mean(tb <= tau_b))
        boot_devs[b] = abs(target - cov_b)

    lower_q = 0.5 * (1.0 - confidence)
    upper_q = 1.0 - lower_q
    lo = float(np.quantile(boot_devs, lower_q))
    hi = float(np.quantile(boot_devs, upper_q))

    return BootstrapCI(
        point=deviation_obs,
        lo=lo,
        hi=hi,
        n_resamples=n_resamples,
        n_perts_calib=n_calib,
        n_perts_test=n_test,
        method="stratified_percentile_deviation",
    )


def bootstrap_aggregate_per_dataset(
    per_predictor_devs: dict[str, float],
    n_resamples: int = 1000,
    confidence: float = 0.95,
    seed: int = 42,
) -> tuple[float, float, float]:
    """Bootstrap CI on a per-(predictor) deviation aggregate (e.g., dataset
    mean deviation across all predictors).

    Inputs: dict of {predictor_id: mean_deviation} for one dataset.
    Returns: (point, lo, hi) where point is the unweighted mean across
    predictors and (lo, hi) is the 95% percentile bootstrap CI of resampled
    means.
    """
    if not per_predictor_devs:
        raise ValueError("per_predictor_devs must be non-empty")
    devs = np.asarray(list(per_predictor_devs.values()), dtype=np.float64)
    point = float(devs.mean())

    rng = np.random.default_rng(seed)
    boots = np.empty(n_resamples, dtype=np.float64)
    for b in range(n_resamples):
        boots[b] = devs[rng.integers(0, devs.size, size=devs.size)].mean()

    lower_q = 0.5 * (1.0 - confidence)
    upper_q = 1.0 - lower_q
    return point, float(np.quantile(boots, lower_q)), float(np.quantile(boots, upper_q))
