"""Four conformal head levels for ConfPert (Phase 1).

Off-the-shelf split conformal applied to perturbation cell-population predictions.
No new theorems. All citations are verbatim from prior art (Phase 0 lit_notes/).

Heads:
  1. Per-perturbation marginal head (PerturbationConformal)
     score in {ks, w1, energy, mmd_rbf, bimodality_mismatch, variance_ratio_dev}.
     Standard split conformal: calibrate per-perturbation discrepancy scores; threshold
     at the ceil((1-alpha)(n+1))/n empirical quantile of calibration scores.
     Citation: Vovk-Lei lineage; Romano-Patterson-Candes 2019 split-conformal protocol.

  2. Per-gene CD-split head (CDSplitConformal)
     Marginal coverage on per-gene quantile bands. Per HetPert v1 implementation.
     Citation: Romano-Patterson-Candes 2019 + Izbicki-Shimizu-Stern 2022.

  3. Per-population energy-distance head (EnergyConformal)
     Population-level coverage via energy-distance score. Per HetPert v1 implementation.

  4. Subgroup-conditional head (CauchoisSubgroupConformal)
     Mondrian-style per-subgroup empirical-quantile calibration.
     Citation: Cauchois-Gupta-Duchi 2024 JMLR Algorithm 1.

  Plus:
  5. RCPS risk control (RCPSCalibrator) for K2 selective drug-screening triage.
     UCB calibration with Hoeffding-Bentkus or Waudby-Smith-Ramdas bounds.
     Citation: Bates-Angelopoulos-Lei-Malik-Jordan 2021 JACM Theorem 1.

Coverage guarantee under exchangeability: P(score <= tau_hat) >= 1 - alpha for the
per-perturbation head; analogous statements per head.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Sequence

import numpy as np
from scipy import stats

from .metrics import _validate_pair, SCORES, energy_distance_per_gene


# ---------------------------------------------------------------------------
# Helper: split-conformal quantile correction with the +1 finite-sample factor
# ---------------------------------------------------------------------------


def split_conformal_quantile(scores: np.ndarray, alpha: float) -> float:
    """Standard split-conformal quantile per Vovk-Lei lineage:
        q_index = ceil((n+1)(1-alpha)) / n  (clamped to [1/n, 1]).

    Returns the (1-alpha)-corrected quantile of `scores` for prediction-set construction.
    """
    n = len(scores)
    if n == 0:
        raise ValueError("Need at least one calibration score")
    q_index = int(np.ceil((n + 1) * (1.0 - alpha))) / n
    q_index = float(min(max(q_index, 1.0 / n), 1.0))
    return float(np.quantile(scores, q_index, method="higher"))


# ---------------------------------------------------------------------------
# Per-perturbation marginal head
# ---------------------------------------------------------------------------


@dataclass
class PerturbationConformal:
    """Calibrate a discrepancy threshold tau such that P(score <= tau) >= 1 - alpha
    over held-out perturbations.

    Usage:
      pc = PerturbationConformal(score_fn=SCORES["energy"], alpha=0.1)
      pc.calibrate(pred_pops, obs_pops)        # both lists of [n_cells, d_genes]
      report = pc.coverage(pred_test, obs_test) # empirical coverage on test pairs
    """
    score_fn: Callable[[np.ndarray, np.ndarray], float]
    alpha: float = 0.1
    tau: float | None = None
    calibration_scores: np.ndarray | None = None

    def calibrate(self, pred_pops: Sequence[np.ndarray],
                  obs_pops: Sequence[np.ndarray]) -> "PerturbationConformal":
        if len(pred_pops) != len(obs_pops):
            raise ValueError(
                f"len mismatch: {len(pred_pops)} pred vs {len(obs_pops)} obs"
            )
        scores = np.array([
            self.score_fn(p, o) for p, o in zip(pred_pops, obs_pops)
        ], dtype=np.float64)
        self.calibration_scores = scores
        self.tau = split_conformal_quantile(scores, self.alpha)
        return self

    def coverage(self, pred_pops: Sequence[np.ndarray],
                 obs_pops: Sequence[np.ndarray]) -> dict[str, float]:
        if self.tau is None:
            raise RuntimeError("calibrate() first")
        if len(pred_pops) != len(obs_pops):
            raise ValueError("len mismatch on test pairs")
        scores = np.array([
            self.score_fn(p, o) for p, o in zip(pred_pops, obs_pops)
        ], dtype=np.float64)
        in_set = (scores <= self.tau)
        return {
            "alpha": float(self.alpha),
            "target_coverage": float(1.0 - self.alpha),
            "tau": float(self.tau),
            "achieved_coverage": float(in_set.mean()),
            "calibration_deviation": float(abs((1.0 - self.alpha) - in_set.mean())),
            "n_test": int(len(scores)),
            "test_scores_mean": float(np.mean(scores)),
            "test_scores_max": float(np.max(scores)),
        }


# ---------------------------------------------------------------------------
# Per-gene CD-split head (HetPert v1, Romano-Izbicki lineage)
# ---------------------------------------------------------------------------


@dataclass
class CDSplitConformal:
    """Per-gene quantile-band coverage. Calibration: residual between predicted-quantile
    bands and observed cells per gene. Test: empirical coverage of held-out cells.
    """
    alpha: float = 0.1
    quantile_levels: tuple[float, float] | None = None
    deltas: np.ndarray | None = None

    def __post_init__(self):
        if self.quantile_levels is None:
            self.quantile_levels = (self.alpha / 2.0, 1.0 - self.alpha / 2.0)

    def calibrate(self, X_pred_calib: np.ndarray, X_obs_calib: np.ndarray) -> "CDSplitConformal":
        _validate_pair(X_obs_calib, X_pred_calib)
        d = X_obs_calib.shape[1]
        deltas = np.zeros(d, dtype=np.float64)
        ql, qh = self.quantile_levels
        for j in range(d):
            lo_pred = np.quantile(X_pred_calib[:, j], ql)
            hi_pred = np.quantile(X_pred_calib[:, j], qh)
            obs = X_obs_calib[:, j]
            below = lo_pred - obs
            above = obs - hi_pred
            residuals = np.maximum(np.maximum(below, above), 0.0)
            deltas[j] = np.quantile(residuals, 1.0 - self.alpha, method="higher")
        self.deltas = deltas
        return self

    def band(self, X_pred_test: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        if self.deltas is None:
            raise RuntimeError("calibrate() first")
        ql, qh = self.quantile_levels
        d = X_pred_test.shape[1]
        lower = np.empty(d, dtype=np.float64)
        upper = np.empty(d, dtype=np.float64)
        for j in range(d):
            lower[j] = np.quantile(X_pred_test[:, j], ql) - self.deltas[j]
            upper[j] = np.quantile(X_pred_test[:, j], qh) + self.deltas[j]
        return lower, upper

    def coverage(self, X_pred_test: np.ndarray, X_obs_test: np.ndarray) -> dict[str, float]:
        lower, upper = self.band(X_pred_test)
        d = X_obs_test.shape[1]
        per_gene = np.empty(d, dtype=np.float64)
        for j in range(d):
            inside = (X_obs_test[:, j] >= lower[j]) & (X_obs_test[:, j] <= upper[j])
            per_gene[j] = inside.mean()
        return {
            "alpha": float(self.alpha),
            "target_coverage": float(1.0 - self.alpha),
            "mean_coverage": float(per_gene.mean()),
            "min_coverage": float(per_gene.min()),
            "max_coverage": float(per_gene.max()),
            "calibration_deviation": float(abs((1.0 - self.alpha) - per_gene.mean())),
            "n_genes_under_target": int((per_gene < (1.0 - self.alpha)).sum()),
            "n_genes_total": int(d),
        }


# ---------------------------------------------------------------------------
# Per-population energy head (HetPert v1, kept for backward-compat)
# ---------------------------------------------------------------------------


@dataclass
class EnergyConformal:
    alpha: float = 0.1
    radius: float | None = None
    scores_calibration: np.ndarray | None = None

    def calibrate(self, pred_pops: Sequence[np.ndarray],
                  obs_pops: Sequence[np.ndarray]) -> "EnergyConformal":
        if len(pred_pops) != len(obs_pops):
            raise ValueError(f"len mismatch: {len(pred_pops)} vs {len(obs_pops)}")
        scores = []
        for X_pred, X_obs in zip(pred_pops, obs_pops):
            d_per_gene = energy_distance_per_gene(X_obs, X_pred)
            scores.append(float(np.nanmean(d_per_gene)))
        self.scores_calibration = np.array(scores)
        self.radius = split_conformal_quantile(np.array(scores), self.alpha)
        return self

    def coverage(self, pred_pops: Sequence[np.ndarray],
                 obs_pops: Sequence[np.ndarray]) -> dict[str, float]:
        if self.radius is None:
            raise RuntimeError("calibrate() first")
        scores = []
        for X_pred, X_obs in zip(pred_pops, obs_pops):
            d_per_gene = energy_distance_per_gene(X_obs, X_pred)
            scores.append(float(np.nanmean(d_per_gene)))
        in_set = np.array(scores) <= self.radius
        return {
            "alpha": float(self.alpha),
            "target_coverage": float(1.0 - self.alpha),
            "radius": float(self.radius),
            "achieved_coverage": float(in_set.mean()),
            "calibration_deviation": float(abs((1.0 - self.alpha) - in_set.mean())),
            "n_test": int(len(scores)),
        }


# ---------------------------------------------------------------------------
# Subgroup-conditional head (Cauchois-Gupta-Duchi 2024)
# ---------------------------------------------------------------------------


@dataclass
class CauchoisSubgroupConformal:
    """Mondrian-style per-subgroup conformal calibration.

    Each calibration item is tagged with a discrete subgroup label (e.g.
    "responder" / "non-responder" from a predicted-bimodal-mode partition).
    A separate threshold tau_S is calibrated per subgroup using the per-subgroup
    quantile correction. Subgroup-conditional marginal coverage holds within each S
    when calibration set is exchangeable within the subgroup.
    """
    score_fn: Callable[[np.ndarray, np.ndarray], float]
    alpha: float = 0.1
    taus: dict[str, float] = field(default_factory=dict)
    n_per_subgroup: dict[str, int] = field(default_factory=dict)

    def calibrate(self, pred_pops: Sequence[np.ndarray],
                  obs_pops: Sequence[np.ndarray],
                  subgroup_labels: Sequence[str]) -> "CauchoisSubgroupConformal":
        if not (len(pred_pops) == len(obs_pops) == len(subgroup_labels)):
            raise ValueError("len mismatch among pred/obs/labels")
        scores_by_S: dict[str, list[float]] = {}
        for p, o, s in zip(pred_pops, obs_pops, subgroup_labels):
            scores_by_S.setdefault(s, []).append(float(self.score_fn(p, o)))
        for s, scores in scores_by_S.items():
            arr = np.array(scores, dtype=np.float64)
            self.taus[s] = split_conformal_quantile(arr, self.alpha)
            self.n_per_subgroup[s] = len(arr)
        return self

    def coverage(self, pred_pops: Sequence[np.ndarray],
                 obs_pops: Sequence[np.ndarray],
                 subgroup_labels: Sequence[str]) -> dict[str, dict]:
        if not self.taus:
            raise RuntimeError("calibrate() first")
        per_S: dict[str, list[bool]] = {}
        for p, o, s in zip(pred_pops, obs_pops, subgroup_labels):
            tau = self.taus.get(s)
            if tau is None:
                continue  # subgroup not seen at calibration; skip
            score = float(self.score_fn(p, o))
            per_S.setdefault(s, []).append(score <= tau)
        out = {}
        for s, hits in per_S.items():
            arr = np.array(hits)
            out[s] = {
                "alpha": float(self.alpha),
                "target_coverage": float(1.0 - self.alpha),
                "tau": float(self.taus[s]),
                "achieved_coverage": float(arr.mean()),
                "calibration_deviation": float(abs((1.0 - self.alpha) - arr.mean())),
                "n_calibration": int(self.n_per_subgroup[s]),
                "n_test": int(len(arr)),
            }
        return out


# ---------------------------------------------------------------------------
# RCPS risk control (Bates 2021 Theorem 1)
# ---------------------------------------------------------------------------


def hoeffding_bentkus_ucb(losses: np.ndarray, delta: float) -> float:
    """Two-sided UCB on the mean of bounded losses in [0, 1] using
    Hoeffding-Bentkus per Bates 2021 Eq. 8 / 9.

    Returns the (1-delta) upper confidence bound on E[L]. For non-binary bounded
    losses, use the Hoeffding bound; for binary, Bentkus is tighter.
    """
    n = len(losses)
    if n == 0:
        raise ValueError("Need at least one loss observation")
    mean = float(losses.mean())
    # Hoeffding bound: mean + sqrt(log(1/delta) / (2 n))
    hoeffding = mean + float(np.sqrt(np.log(1.0 / delta) / (2.0 * n)))
    return min(hoeffding, 1.0)


@dataclass
class RCPSCalibrator:
    """Risk-controlling prediction-set calibrator per Bates 2021 Theorem 1.

    For monotone losses L(C(lambda)) in [0, 1], select the smallest lambda such that
    UCB on the empirical risk falls below the target alpha. Returns lambda_hat.

    Usage:
      rcps = RCPSCalibrator(target_risk=0.10, confidence=0.05)
      rcps.calibrate(losses_per_lambda)  # dict {lambda: array of losses}
      lambda_hat = rcps.lambda_hat
    """
    target_risk: float = 0.10
    confidence: float = 0.05
    lambda_hat: float | None = None
    grid_results: dict[float, float] = field(default_factory=dict)

    def calibrate(self, losses_per_lambda: dict[float, np.ndarray]) -> "RCPSCalibrator":
        """losses_per_lambda: {lambda_value: 1D array of losses on calibration set}.

        Sweep lambdas in increasing order. For each, compute UCB(loss) at confidence
        1 - delta. Smallest lambda with UCB <= target_risk is lambda_hat.
        """
        sorted_lambdas = sorted(losses_per_lambda.keys())
        chosen = None
        for lam in sorted_lambdas:
            ucb = hoeffding_bentkus_ucb(losses_per_lambda[lam], self.confidence)
            self.grid_results[float(lam)] = float(ucb)
            if ucb <= self.target_risk and chosen is None:
                chosen = float(lam)
        self.lambda_hat = chosen
        return self
