"""Unit tests for confpert.predictors and confpert.noise_models."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from confpert.noise_models import NOISE_VARIANTS, apply_noise
from confpert.predictors import (
    AdditivePredictor,
    AhlmannEltzeBilinearRidge,
    MeanPredictor,
    NoisyMeanPredictor,
)


@pytest.fixture
def fake_data():
    rng = np.random.RandomState(0)
    n_ctrl = 200
    n_per_pert = 100
    d_genes = 16
    perts = ["A", "B", "C"]
    X_ctrl = rng.randn(n_ctrl, d_genes).astype(np.float32)
    X_pert = []
    p_labels = []
    for p in perts:
        shift = rng.randn(d_genes) * 0.3
        X_p = (rng.randn(n_per_pert, d_genes) + shift).astype(np.float32)
        X_pert.append(X_p)
        p_labels.extend([p] * n_per_pert)
    X_pert = np.concatenate(X_pert)
    p_labels = np.array(p_labels)
    return X_ctrl, X_pert, p_labels, perts


def test_noise_variant_a_no_noise(fake_data):
    _, X_pert, p_labels, perts = fake_data
    mu = X_pert[p_labels == "A"].mean(axis=0)
    samples = apply_noise(mu, n_cells=50, variant="A_no_noise")
    assert samples.shape == (50, mu.shape[0])
    # Variance ratio: 0
    assert samples.var(axis=0).max() < 1e-9


def test_noise_variant_b_isotropic(fake_data):
    _, X_pert, p_labels, _ = fake_data
    mu = X_pert[p_labels == "A"].mean(axis=0)
    X_train = X_pert[p_labels == "A"]
    samples = apply_noise(mu, n_cells=200, variant="B_isotropic", X_train_pert=X_train, seed=42)
    assert samples.shape == (200, mu.shape[0])
    assert samples.var() > 0.5


def test_noise_variant_c_per_gene(fake_data):
    _, X_pert, p_labels, _ = fake_data
    mu = X_pert[p_labels == "A"].mean(axis=0)
    X_train = X_pert[p_labels == "A"]
    samples = apply_noise(mu, n_cells=200, variant="C_per_gene_marginal",
                          X_train_pert=X_train, seed=42)
    assert samples.shape == (200, mu.shape[0])


def test_noise_variant_d_full_cov(fake_data):
    _, X_pert, p_labels, _ = fake_data
    mu = X_pert[p_labels == "A"].mean(axis=0)
    X_train = X_pert[p_labels == "A"]
    samples = apply_noise(mu, n_cells=200, variant="D_full_covariance",
                          X_train_pert=X_train, seed=42)
    assert samples.shape == (200, mu.shape[0])


def test_mean_predictor_round_trip(fake_data):
    X_ctrl, X_pert, p_labels, perts = fake_data
    pred = MeanPredictor(noise_variant="A_no_noise").fit(X_ctrl, p_labels, X_pert)
    out = pred.predict_samples(X_ctrl[:50], "A", n_cells=50)
    assert out.shape == (50, X_pert.shape[1])
    # Variant A means all rows identical
    assert np.allclose(out[0], out[10])


def test_mean_predictor_with_noise_c(fake_data):
    X_ctrl, X_pert, p_labels, _ = fake_data
    pred = MeanPredictor(noise_variant="C_per_gene_marginal").fit(X_ctrl, p_labels, X_pert)
    out = pred.predict_samples(X_ctrl[:50], "A", n_cells=200)
    assert out.var(axis=0).max() > 1e-3  # noise was added


def test_ahlmann_bilinear_ridge_smoke(fake_data):
    X_ctrl, X_pert, p_labels, perts = fake_data
    pred = AhlmannEltzeBilinearRidge(K=3, lam=0.1).fit(X_ctrl, p_labels, X_pert)
    mu = pred.predict_mean("A")
    assert mu.shape == (X_pert.shape[1],)
    out = pred.predict_samples(X_ctrl[:50], "A", n_cells=50)
    assert out.shape == (50, X_pert.shape[1])


def test_additive_predictor_round_trip(fake_data):
    X_ctrl, X_pert, p_labels, _ = fake_data
    pred = AdditivePredictor().fit(X_ctrl, p_labels, X_pert)
    out = pred.predict_samples(X_ctrl[:50], "A")
    assert out.shape == (50, X_pert.shape[1])


def test_noisy_mean_predictor_has_variance(fake_data):
    X_ctrl, X_pert, p_labels, _ = fake_data
    pred = NoisyMeanPredictor().fit(X_ctrl, p_labels, X_pert)
    out = pred.predict_samples(X_ctrl[:50], "A", n_cells=200)
    assert out.var(axis=0).max() > 1e-3
