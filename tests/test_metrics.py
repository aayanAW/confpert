"""Unit tests for confpert.metrics: six distributional discrepancies."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from confpert.metrics import (
    SCORES,
    bimodality_coef_per_gene,
    bimodality_match,
    energy_distance_per_gene,
    evaluate_six,
    ks_per_gene,
    mmd_rbf,
    median_bandwidth,
    score_bimodality_mismatch,
    score_energy,
    score_ks,
    score_mmd_rbf,
    score_variance_ratio_dev,
    score_w1,
    variance_ratio_per_gene,
    wasserstein1_per_gene,
)


@pytest.fixture
def rng():
    return np.random.RandomState(0)


# Per-gene 1D ----------------------------------------------------------------


def test_ks_identical_zero(rng):
    X = rng.randn(200, 8).astype(np.float32)
    out = ks_per_gene(X, X)
    assert out.shape == (8,)
    assert np.all(out < 1e-9)


def test_w1_identical_zero(rng):
    X = rng.randn(200, 8).astype(np.float32)
    assert np.all(wasserstein1_per_gene(X, X) < 1e-9)


def test_energy_identical_zero(rng):
    X = rng.randn(200, 8).astype(np.float32)
    assert np.all(energy_distance_per_gene(X, X) < 1e-6)


def test_variance_ratio_collapsed(rng):
    X_obs = rng.randn(200, 8).astype(np.float32)
    X_pred = np.tile(X_obs.mean(axis=0, keepdims=True), (200, 1))
    vr = variance_ratio_per_gene(X_obs, X_pred)
    assert np.all(vr < 1e-6)


def test_variance_ratio_perfect(rng):
    X_obs = rng.randn(500, 8).astype(np.float32)
    X_pred = rng.randn(500, 8).astype(np.float32)
    vr = variance_ratio_per_gene(X_obs, X_pred)
    assert np.all((vr > 0.5) & (vr < 2.0))


def test_bimodality_unimodal(rng):
    X = rng.randn(2000, 4).astype(np.float32)
    b = bimodality_coef_per_gene(X)
    assert np.all(b < 0.7)


def test_bimodality_strongly_bimodal(rng):
    n = 1000
    half = n // 2
    X = np.concatenate([
        rng.randn(half, 4) - 4.0,
        rng.randn(half, 4) + 4.0,
    ]).astype(np.float32)
    b = bimodality_coef_per_gene(X)
    assert np.all(b > 5.0 / 9.0), f"bimodal data b = {b}"


def test_bimodality_match_perfect(rng):
    n = 500
    half = n // 2
    X = np.concatenate([
        rng.randn(half, 4) - 3.0,
        rng.randn(half, 4) + 3.0,
    ]).astype(np.float32)
    out = bimodality_match(X, X.copy())
    assert out["accuracy"] > 0.9


def test_bimodality_match_collapsed(rng):
    n = 500
    X_obs = np.concatenate([
        rng.randn(n // 2, 4) - 3.0,
        rng.randn(n // 2, 4) + 3.0,
    ]).astype(np.float32)
    X_pred = rng.randn(n, 4).astype(np.float32)
    out = bimodality_match(X_obs, X_pred)
    assert out["n_obs_bimodal"] == 4
    assert out["sensitivity"] < 0.5


# MMD-RBF -------------------------------------------------------------------


def test_mmd_identical_near_zero(rng):
    X = rng.randn(200, 8).astype(np.float32)
    val = mmd_rbf(X, X.copy())
    assert val < 1e-6 + 0.05  # allow tiny estimator noise


def test_mmd_distinct_positive(rng):
    X = rng.randn(200, 8).astype(np.float32)
    Y = rng.randn(200, 8).astype(np.float32) + 5.0
    assert mmd_rbf(X, Y) > 0.05


def test_median_bandwidth_positive(rng):
    X = rng.randn(100, 8).astype(np.float32)
    Y = rng.randn(100, 8).astype(np.float32)
    sigma = median_bandwidth(X, Y)
    assert sigma > 0


# Aggregator ----------------------------------------------------------------


def test_evaluate_six_smoke(rng):
    X_obs = rng.randn(300, 16).astype(np.float32)
    X_pred = rng.randn(300, 16).astype(np.float32)
    rep = evaluate_six(X_obs, X_pred)
    assert rep.n_obs == 300
    assert rep.n_pred == 300
    assert rep.d_genes == 16
    assert 0.0 < rep.variance_ratio_mean < 2.0
    assert rep.ks_mean >= 0.0
    assert rep.mmd_rbf >= -0.1  # unbiased estimator can dip below 0 slightly
    assert rep.w1_mean >= 0.0
    assert rep.energy_mean >= 0.0


def test_evaluate_mean_collapse_signature(rng):
    X_obs = rng.randn(200, 16).astype(np.float32)
    X_pred = np.tile(X_obs.mean(axis=0, keepdims=True), (200, 1))
    rep = evaluate_six(X_obs, X_pred)
    assert rep.pcc_pseudobulk > 0.99
    assert rep.variance_ratio_mean < 1e-3
    assert rep.n_genes_collapsed_pct > 0.99


# Score functions for conformal calibration ---------------------------------


def test_all_six_scores_run(rng):
    X_obs = rng.randn(100, 8).astype(np.float32)
    X_pred = rng.randn(100, 8).astype(np.float32)
    for name, fn in SCORES.items():
        s = fn(X_pred, X_obs)
        assert isinstance(s, float), f"{name} returned non-float"
        assert not np.isnan(s) or name in {"variance_ratio_dev"}, name


def test_score_signs(rng):
    X_obs = rng.randn(100, 8).astype(np.float32)
    s_match = score_ks(X_obs.copy(), X_obs)
    s_distinct = score_ks(X_obs.copy() + 5.0, X_obs)
    assert s_match < s_distinct
