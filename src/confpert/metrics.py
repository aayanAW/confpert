"""Six distributional discrepancies for ConfPert.

Per preregistration.md, ConfPert reports six per-(predictor, dataset, split) discrepancies:

  1. KS                     (per-gene 1D Kolmogorov-Smirnov, mean over genes)
  2. Wasserstein-1          (per-gene 1D earth-mover, mean over genes)
  3. Energy distance        (per-gene Szekely-Rizzo, mean over genes)
  4. MMD-RBF                (per-population, median-bandwidth heuristic)
  5. Bimodality coef match  (SAS bimodality coefficient classification accuracy at b > 5/9)
  6. Variance ratio         (per-gene var(X_pred)/var(X_obs), mean over genes)

All discrepancies are score functions s(X_pred, X_obs) -> R that ConfPert wraps with
finite-sample conformal coverage via split conformal calibration.

References:
  KS / Wasserstein -- standard 1D divergences (scipy).
  Energy -- Szekely & Rizzo 2004.
  MMD-RBF -- Gretton et al. 2012, with median-bandwidth heuristic per Sriperumbudur 2009.
  Bimodality coef -- SAS Ratemaker; Pfister 2013. b > 5/9 = bimodal.
  Variance ratio -- HetPert lineage; central diagnostic for mean-collapse.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import stats


# ---------------------------------------------------------------------------
# Per-gene 1D discrepancies
# ---------------------------------------------------------------------------


def ks_per_gene(X_obs: np.ndarray, X_pred: np.ndarray) -> np.ndarray:
    """Two-sample KS per gene. Returns [d_genes] in [0, 1]."""
    _validate_pair(X_obs, X_pred)
    d = X_obs.shape[1]
    out = np.empty(d, dtype=np.float64)
    for j in range(d):
        out[j] = stats.ks_2samp(X_obs[:, j], X_pred[:, j], method="asymp").statistic
    return out


def wasserstein1_per_gene(X_obs: np.ndarray, X_pred: np.ndarray) -> np.ndarray:
    """Wasserstein-1 per gene. Returns [d_genes] in expression units."""
    _validate_pair(X_obs, X_pred)
    d = X_obs.shape[1]
    out = np.empty(d, dtype=np.float64)
    for j in range(d):
        out[j] = stats.wasserstein_distance(X_obs[:, j], X_pred[:, j])
    return out


def energy_distance_per_gene(X_obs: np.ndarray, X_pred: np.ndarray) -> np.ndarray:
    """Energy distance per gene (Szekely-Rizzo). Returns [d_genes], non-negative."""
    _validate_pair(X_obs, X_pred)
    d = X_obs.shape[1]
    out = np.empty(d, dtype=np.float64)
    for j in range(d):
        out[j] = stats.energy_distance(X_obs[:, j], X_pred[:, j])
    return out


def variance_ratio_per_gene(X_obs: np.ndarray, X_pred: np.ndarray,
                            eps: float = 1e-9) -> np.ndarray:
    """Per-gene var(X_pred)/var(X_obs). 1.0 = perfect; <1 mean-collapse; >1 over-dispersed.

    Returns [d_genes]. NaN where observed variance < eps (uninformative gene).
    """
    _validate_pair(X_obs, X_pred)
    v_obs = X_obs.var(axis=0, ddof=1)
    v_pred = X_pred.var(axis=0, ddof=1)
    safe = v_obs > eps
    out = np.full(X_obs.shape[1], np.nan, dtype=np.float64)
    out[safe] = v_pred[safe] / v_obs[safe]
    return out


def bimodality_coef_per_gene(X: np.ndarray) -> np.ndarray:
    """SAS bimodality coefficient: b = (skew^2 + 1) / (kurt_excess + 3(n-1)^2/((n-2)(n-3))).

    b > 5/9 ~= 0.555 indicates bimodal-or-multimodal distribution.
    Returns [d_genes] in (0, 1].
    """
    n = X.shape[0]
    if n < 4:
        return np.zeros(X.shape[1], dtype=np.float64)
    skew = stats.skew(X, axis=0, bias=False)
    kurt = stats.kurtosis(X, axis=0, bias=False)
    denom = kurt + 3.0 * (n - 1) ** 2 / ((n - 2) * (n - 3))
    denom = np.where(np.abs(denom) < 1e-12, 1e-12, denom)
    return (skew ** 2 + 1.0) / denom


def bimodality_match(X_obs: np.ndarray, X_pred: np.ndarray,
                     threshold: float = 5.0 / 9.0) -> dict[str, float]:
    """Per-gene bimodality classification match.

    Returns {accuracy, sensitivity (recall on observed-bimodal),
             specificity (recall on observed-unimodal), n_obs_bimodal}.
    """
    _validate_pair(X_obs, X_pred)
    b_obs = bimodality_coef_per_gene(X_obs) > threshold
    b_pred = bimodality_coef_per_gene(X_pred) > threshold
    n_obs_bi = int(b_obs.sum())
    n_obs_uni = int((~b_obs).sum())
    accuracy = float((b_obs == b_pred).mean())
    sensitivity = float((b_pred[b_obs].mean()) if n_obs_bi > 0 else 0.0)
    specificity = float(((~b_pred[~b_obs]).mean()) if n_obs_uni > 0 else 0.0)
    return {"accuracy": accuracy, "sensitivity": sensitivity,
            "specificity": specificity, "n_obs_bimodal": n_obs_bi}


# ---------------------------------------------------------------------------
# Per-population MMD-RBF (Gretton 2012, median-bandwidth heuristic)
# ---------------------------------------------------------------------------


def _pairwise_sq_dists(X: np.ndarray, Y: np.ndarray) -> np.ndarray:
    """Pairwise squared Euclidean distances. Returns [n_X, n_Y]."""
    sq_X = (X ** 2).sum(axis=1, keepdims=True)
    sq_Y = (Y ** 2).sum(axis=1, keepdims=True).T
    cross = X @ Y.T
    out = sq_X + sq_Y - 2 * cross
    return np.maximum(out, 0.0)


def median_bandwidth(X: np.ndarray, Y: np.ndarray, max_pts: int = 500,
                     seed: int = 0) -> float:
    """Median-bandwidth heuristic for RBF kernel: gamma = 1 / (2 sigma^2),
    sigma = median pairwise Euclidean distance.

    Subsamples up to max_pts cells per population for numerical stability.
    """
    rng = np.random.RandomState(seed)
    n_X = X.shape[0]
    n_Y = Y.shape[0]
    if n_X > max_pts:
        X = X[rng.choice(n_X, max_pts, replace=False)]
    if n_Y > max_pts:
        Y = Y[rng.choice(n_Y, max_pts, replace=False)]
    Z = np.vstack([X, Y])
    n_Z = Z.shape[0]
    if n_Z > max_pts:
        Z = Z[rng.choice(n_Z, max_pts, replace=False)]
    sq = _pairwise_sq_dists(Z, Z)
    iu = np.triu_indices(Z.shape[0], k=1)
    median_sq = float(np.median(sq[iu]))
    return float(np.sqrt(max(median_sq, 1e-12)))


def mmd_rbf(X_obs: np.ndarray, X_pred: np.ndarray, sigma: float | None = None,
            max_pts: int = 500, seed: int = 0) -> float:
    """Squared MMD with RBF kernel (Gretton 2012 unbiased estimator).

    If sigma is None, use median-bandwidth heuristic computed on X_obs U X_pred.
    Returns scalar MMD^2 (non-negative up to estimator noise).
    """
    _validate_pair(X_obs, X_pred)
    rng = np.random.RandomState(seed)
    n_X = X_obs.shape[0]
    n_Y = X_pred.shape[0]
    if n_X > max_pts:
        X_obs = X_obs[rng.choice(n_X, max_pts, replace=False)]
        n_X = X_obs.shape[0]
    if n_Y > max_pts:
        X_pred = X_pred[rng.choice(n_Y, max_pts, replace=False)]
        n_Y = X_pred.shape[0]
    if sigma is None:
        sigma = median_bandwidth(X_obs, X_pred, max_pts=max_pts, seed=seed)
    gamma = 1.0 / (2.0 * sigma ** 2 + 1e-12)
    K_xx = np.exp(-gamma * _pairwise_sq_dists(X_obs, X_obs))
    K_yy = np.exp(-gamma * _pairwise_sq_dists(X_pred, X_pred))
    K_xy = np.exp(-gamma * _pairwise_sq_dists(X_obs, X_pred))
    # Unbiased estimator (zero diagonal)
    np.fill_diagonal(K_xx, 0.0)
    np.fill_diagonal(K_yy, 0.0)
    if n_X < 2 or n_Y < 2:
        return float("nan")
    term_xx = K_xx.sum() / (n_X * (n_X - 1))
    term_yy = K_yy.sum() / (n_Y * (n_Y - 1))
    term_xy = K_xy.mean()
    return float(term_xx + term_yy - 2.0 * term_xy)


# ---------------------------------------------------------------------------
# Aggregated report
# ---------------------------------------------------------------------------


@dataclass
class DiscrepancyReport:
    """One row per (model, perturbation) of all six discrepancies plus PCC backstop."""
    n_obs: int
    n_pred: int
    d_genes: int

    ks_mean: float
    w1_mean: float
    energy_mean: float
    mmd_rbf: float

    variance_ratio_mean: float
    n_genes_collapsed_pct: float

    bimod_accuracy: float
    bimod_sensitivity: float
    bimod_specificity: float
    n_obs_bimodal_genes: int

    pcc_pseudobulk: float

    def to_dict(self) -> dict:
        return self.__dict__


def evaluate_six(X_obs: np.ndarray, X_pred: np.ndarray,
                 mmd_max_pts: int = 500, seed: int = 0) -> DiscrepancyReport:
    """Compute all six discrepancies plus a PCC pseudo-bulk backstop for backward
    compatibility with mean-only papers.

    X_obs, X_pred: [n_cells, d_genes]. Different n_cells per population is fine.
    """
    _validate_pair(X_obs, X_pred)

    ks = ks_per_gene(X_obs, X_pred)
    w1 = wasserstein1_per_gene(X_obs, X_pred)
    en = energy_distance_per_gene(X_obs, X_pred)
    mmd = mmd_rbf(X_obs, X_pred, max_pts=mmd_max_pts, seed=seed)
    vr = variance_ratio_per_gene(X_obs, X_pred)
    bimod = bimodality_match(X_obs, X_pred)

    mu_obs = X_obs.mean(axis=0)
    mu_pred = X_pred.mean(axis=0)
    if mu_obs.std() < 1e-9 or mu_pred.std() < 1e-9:
        pcc = 0.0
    else:
        pcc = float(np.corrcoef(mu_obs, mu_pred)[0, 1])

    vr_safe = vr[~np.isnan(vr)]
    return DiscrepancyReport(
        n_obs=int(X_obs.shape[0]),
        n_pred=int(X_pred.shape[0]),
        d_genes=int(X_obs.shape[1]),

        ks_mean=float(np.nanmean(ks)),
        w1_mean=float(np.nanmean(w1)),
        energy_mean=float(np.nanmean(en)),
        mmd_rbf=float(mmd),

        variance_ratio_mean=float(vr_safe.mean()) if vr_safe.size else float("nan"),
        n_genes_collapsed_pct=float((vr_safe < 0.1).mean()) if vr_safe.size else float("nan"),

        bimod_accuracy=bimod["accuracy"],
        bimod_sensitivity=bimod["sensitivity"],
        bimod_specificity=bimod["specificity"],
        n_obs_bimodal_genes=bimod["n_obs_bimodal"],

        pcc_pseudobulk=pcc,
    )


# Six named score functions for conformal calibration.
# Each takes (X_pred, X_obs) and returns a scalar: lower = better.
# Note ordering: X_pred first per typical conformal "score(predicted, observed)" convention.

def score_ks(X_pred: np.ndarray, X_obs: np.ndarray) -> float:
    return float(np.nanmean(ks_per_gene(X_obs, X_pred)))


def score_w1(X_pred: np.ndarray, X_obs: np.ndarray) -> float:
    return float(np.nanmean(wasserstein1_per_gene(X_obs, X_pred)))


def score_energy(X_pred: np.ndarray, X_obs: np.ndarray) -> float:
    return float(np.nanmean(energy_distance_per_gene(X_obs, X_pred)))


def score_mmd_rbf(X_pred: np.ndarray, X_obs: np.ndarray) -> float:
    return mmd_rbf(X_obs, X_pred)


def score_bimodality_mismatch(X_pred: np.ndarray, X_obs: np.ndarray) -> float:
    """1 - bimodality classification accuracy. Lower = better (i.e. higher accuracy)."""
    m = bimodality_match(X_obs, X_pred)
    return 1.0 - m["accuracy"]


def score_variance_ratio_dev(X_pred: np.ndarray, X_obs: np.ndarray) -> float:
    """|VR - 1| mean over genes. Lower = better (closer to perfect variance match)."""
    vr = variance_ratio_per_gene(X_obs, X_pred)
    vr_safe = vr[~np.isnan(vr)]
    if vr_safe.size == 0:
        return float("nan")
    return float(np.mean(np.abs(vr_safe - 1.0)))


SCORES = {
    "ks": score_ks,
    "w1": score_w1,
    "energy": score_energy,
    "mmd_rbf": score_mmd_rbf,
    "bimodality_mismatch": score_bimodality_mismatch,
    "variance_ratio_dev": score_variance_ratio_dev,
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _validate_pair(X_obs: np.ndarray, X_pred: np.ndarray) -> None:
    if X_obs.ndim != 2 or X_pred.ndim != 2:
        raise ValueError(f"Expected 2D arrays, got {X_obs.shape} and {X_pred.shape}")
    if X_obs.shape[1] != X_pred.shape[1]:
        raise ValueError(
            f"Gene dim mismatch: obs has {X_obs.shape[1]}, pred has {X_pred.shape[1]}"
        )
    if X_obs.shape[0] < 2 or X_pred.shape[0] < 2:
        raise ValueError(
            f"Need >=2 cells per population: obs={X_obs.shape[0]}, pred={X_pred.shape[0]}"
        )
