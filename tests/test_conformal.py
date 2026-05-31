"""Unit tests for confpert.conformal: four conformal head levels."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from confpert.conformal import (
    CauchoisSubgroupConformal,
    CDSplitConformal,
    EnergyConformal,
    PerturbationConformal,
    WeightedPerturbationConformal,
    RCPSCalibrator,
    hoeffding_bentkus_ucb,
    split_conformal_quantile,
)
from confpert.metrics import SCORES


@pytest.fixture
def rng():
    return np.random.RandomState(0)


def test_split_quantile_basic():
    scores = np.array([0.1, 0.2, 0.3, 0.4, 0.5])
    q = split_conformal_quantile(scores, alpha=0.2)
    # n=5, q_index = ceil(6 * 0.8)/5 = 5/5 = 1.0 -> max
    assert q == 0.5


def test_split_quantile_low_alpha():
    scores = np.arange(100).astype(np.float64)
    q = split_conformal_quantile(scores, alpha=0.05)
    # ceil(101*0.95)/100 = 96/100 -> 96th percentile
    assert q >= 95.0


def test_perturbation_conformal_marginal_coverage(rng):
    # Generate 200 calibration perts and 200 test perts from the SAME distribution.
    # Conformal coverage should be ~ 1 - alpha.
    pred_pops = []
    obs_pops = []
    for _ in range(400):
        # Both pred and obs drawn from same Gaussian, but with shift of ~variable size
        n = 50
        mu = rng.randn(8) * 0.5
        X_obs = rng.randn(n, 8) + mu
        X_pred = rng.randn(n, 8) + mu + 0.1 * rng.randn(8)  # mild misalignment
        pred_pops.append(X_pred.astype(np.float32))
        obs_pops.append(X_obs.astype(np.float32))

    pred_calib, pred_test = pred_pops[:200], pred_pops[200:]
    obs_calib, obs_test = obs_pops[:200], obs_pops[200:]

    pc = PerturbationConformal(score_fn=SCORES["energy"], alpha=0.1)
    pc.calibrate(pred_calib, obs_calib)
    cov = pc.coverage(pred_test, obs_test)
    # Expect achieved coverage in [0.85, 1.0]
    assert 0.85 <= cov["achieved_coverage"] <= 1.0, cov


def test_cd_split_marginal_coverage(rng):
    # Calibration + test populations from same distribution, large n_cells.
    X_calib_obs = rng.randn(500, 16).astype(np.float32)
    X_calib_pred = rng.randn(500, 16).astype(np.float32)
    X_test_obs = rng.randn(500, 16).astype(np.float32)
    X_test_pred = rng.randn(500, 16).astype(np.float32)

    cd = CDSplitConformal(alpha=0.1)
    cd.calibrate(X_calib_pred, X_calib_obs)
    cov = cd.coverage(X_test_pred, X_test_obs)
    assert 0.80 <= cov["mean_coverage"] <= 1.0, cov


def test_energy_conformal_basic(rng):
    pred_pops = [rng.randn(60, 8).astype(np.float32) for _ in range(80)]
    obs_pops = [rng.randn(60, 8).astype(np.float32) for _ in range(80)]
    pred_calib, pred_test = pred_pops[:40], pred_pops[40:]
    obs_calib, obs_test = obs_pops[:40], obs_pops[40:]

    ec = EnergyConformal(alpha=0.1)
    ec.calibrate(pred_calib, obs_calib)
    cov = ec.coverage(pred_test, obs_test)
    assert 0.80 <= cov["achieved_coverage"] <= 1.0


def test_subgroup_conformal_per_subgroup(rng):
    pred_pops = [rng.randn(50, 8).astype(np.float32) for _ in range(100)]
    obs_pops = [rng.randn(50, 8).astype(np.float32) for _ in range(100)]
    labels = ["responder"] * 50 + ["non_responder"] * 50
    rng.shuffle(labels)

    sc = CauchoisSubgroupConformal(score_fn=SCORES["energy"], alpha=0.1)
    sc.calibrate(pred_pops[:60], obs_pops[:60], labels[:60])
    cov = sc.coverage(pred_pops[60:], obs_pops[60:], labels[60:])
    for s, info in cov.items():
        assert 0.70 <= info["achieved_coverage"] <= 1.0, (s, info)


def test_hoeffding_bentkus_ucb_monotone():
    # UCB should decrease as n grows, given same losses.
    losses = np.array([0.05] * 100)
    ucb_100 = hoeffding_bentkus_ucb(losses, delta=0.05)
    ucb_10000 = hoeffding_bentkus_ucb(np.array([0.05] * 10000), delta=0.05)
    assert ucb_10000 < ucb_100


def test_rcps_calibrator_picks_smallest_lambda():
    # Risks decrease as lambda grows. UCB at lambda=2.0 falls below target before lambda=1.0.
    losses_per_lambda = {
        0.5: np.array([0.5] * 1000),
        1.0: np.array([0.3] * 1000),
        2.0: np.array([0.05] * 1000),
        3.0: np.array([0.01] * 1000),
    }
    rcps = RCPSCalibrator(target_risk=0.10, confidence=0.05)
    rcps.calibrate(losses_per_lambda)
    assert rcps.lambda_hat == 2.0


# ---------------------------------------------------------------------------
# Weighted split-conformal under covariate shift (Tibshirani et al. 2019)
# ---------------------------------------------------------------------------


def _shift_setup(seed=42, n=400, shift=0.6):
    """Calibration scores from the source domain; test scores shifted larger
    (the covariate-shift signature). Weights are the likelihood ratio
    dP_test/dP_calib(x) propto exp(shift * x) for a Gaussian location shift."""
    rng = np.random.default_rng(seed)
    cal = rng.normal(0.0, 1.0, size=n)
    test = rng.normal(shift, 1.0, size=n)
    w = np.exp(shift * cal)
    return cal, test, w


def test_weighted_cp_recovers_coverage_under_shift():
    cal, test, w = _shift_setup()
    head = WeightedPerturbationConformal(alpha=0.2).calibrate(cal, weights=w)
    cov = head.coverage(test)
    assert 0.75 <= cov <= 0.90  # nominal 0.80 restored (conservative upper edge)


def test_weighted_cp_beats_unweighted_under_shift():
    cal, test, w = _shift_setup()
    weighted = WeightedPerturbationConformal(alpha=0.2).calibrate(cal, weights=w).coverage(test)
    tau_unw = split_conformal_quantile(cal, alpha=0.2)
    unweighted = float(np.mean(test <= tau_unw))
    # Unweighted split-conformal undercovers under the shift; weighting recovers it.
    assert unweighted < 0.75
    assert weighted > unweighted
    assert abs(weighted - 0.80) < abs(unweighted - 0.80)


def test_weighted_cp_uniform_weights_match_unweighted():
    rng = np.random.default_rng(7)
    cal = rng.normal(size=500)
    test = rng.normal(size=4000)
    head = WeightedPerturbationConformal(alpha=0.2).calibrate(cal, weights=np.ones_like(cal))
    tau_std = split_conformal_quantile(cal, alpha=0.2)
    # With uniform weights the weighted head reduces to standard split-conformal
    # (within one order statistic of the (n+1)-corrected quantile).
    assert abs(head.tau_ - tau_std) < 0.05
    w_cov = head.coverage(test)
    std_cov = float(np.mean(test <= tau_std))
    assert abs(w_cov - std_cov) < 0.03


def test_weighted_cp_requires_calibration():
    head = WeightedPerturbationConformal(alpha=0.1)
    with pytest.raises(RuntimeError):
        head.coverage(np.array([0.1, 0.2]))


def test_weighted_cp_rejects_negative_weights():
    head = WeightedPerturbationConformal(alpha=0.1)
    with pytest.raises(ValueError):
        head.calibrate(np.array([0.1, 0.2, 0.3]), weights=np.array([1.0, -1.0, 2.0]))
