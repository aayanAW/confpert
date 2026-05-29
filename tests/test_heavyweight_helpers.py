"""Tests for the train-fold-only variance + sampler helpers in
`confpert.heavyweight_helpers`.

These exist specifically to catch the data-leakage failure mode that
broke scGPT v1 (2026-05-25 incident: wrapper used test-fold variance
for sampling noise, would have inflated coverage). Every heavyweight
predictor MUST flow through these helpers.
"""
from __future__ import annotations

import numpy as np
import pytest

from confpert.heavyweight_helpers import (
    safe_train_variance,
    sample_from_train_distribution,
    make_train_test_split_by_perturbation,
)


def test_safe_train_variance_basic_correctness():
    rng = np.random.default_rng(0)
    X = rng.standard_normal(size=(200, 5)).astype(np.float64)
    # Pre-rotated so the first axis has higher variance than the others
    X[:, 0] *= 3.0
    var = safe_train_variance(X, train_indices=np.arange(0, 100))
    assert var.shape == (5,)
    assert var[0] > var[1:].max()


def test_safe_train_variance_refuses_train_test_overlap():
    X = np.zeros((50, 3))
    with pytest.raises(ValueError, match="overlapping"):
        safe_train_variance(
            X, train_indices=np.arange(0, 35), test_indices=np.arange(30, 45)
        )


def test_safe_train_variance_refuses_too_small_train_fold():
    X = np.zeros((50, 3))
    with pytest.raises(ValueError, match="min_train_n"):
        safe_train_variance(X, train_indices=np.arange(0, 5))


def test_safe_train_variance_clips_zero_variance():
    X = np.ones((100, 3))  # all rows identical -> per-gene variance 0
    var = safe_train_variance(X, train_indices=np.arange(0, 100))
    assert (var > 0).all()
    assert var.max() < 1e-3  # close to the clip floor


def test_safe_train_variance_out_of_range_index():
    X = np.zeros((50, 3))
    with pytest.raises(IndexError, match="out of range"):
        safe_train_variance(X, train_indices=np.array([0, 1, 60]))


def test_sample_from_train_distribution_returns_correct_shape():
    rng = np.random.default_rng(1)
    mu = np.array([1.0, 2.0, 3.0])
    var = np.array([0.5, 1.0, 2.0])
    out = sample_from_train_distribution(mu, var, n_cells=128, rng=rng)
    assert out.shape == (128, 3)
    # Sample mean should be roughly mu within a few SE
    se = np.sqrt(var / 128)
    assert np.all(np.abs(out.mean(axis=0) - mu) < 5 * se)


def test_sample_from_train_distribution_rejects_negative_variance():
    mu = np.array([0.0, 0.0])
    var = np.array([-1.0, 0.5])
    with pytest.raises(ValueError, match="negative"):
        sample_from_train_distribution(mu, var, n_cells=10)


def test_sample_from_train_distribution_rejects_shape_mismatch():
    mu = np.zeros(5)
    var = np.ones(3)
    with pytest.raises(ValueError, match="match"):
        sample_from_train_distribution(mu, var, n_cells=10)


def test_make_train_test_split_by_perturbation_no_overlap_invariant():
    labels = np.array(["A", "B", "C", "A", "B", "C", "ctrl", "ctrl", "ctrl"])
    train_idx, test_idx = make_train_test_split_by_perturbation(
        labels, test_perturbations={"A", "C"}
    )
    overlap = np.intersect1d(train_idx, test_idx)
    assert overlap.size == 0
    test_labels = labels[test_idx]
    assert set(test_labels.tolist()) <= {"A", "C"}


def test_make_train_test_split_by_perturbation_with_ctrl_mask_excludes_ctrl_from_test():
    """When is_ctrl_mask is given, control cells should appear only in train."""
    labels = np.array(["A", "ctrl", "B", "ctrl"])
    is_ctrl = np.array([False, True, False, True])
    # Ctrl labeled with "ctrl" string would normally be in train anyway because
    # is_ctrl is in `test_perturbations={"A"}` -> in test by string.
    # The is_ctrl_mask path is exercised when controls happen to match the
    # test-perturbation string. So construct a labels vector where "A" cells
    # include a control:
    labels = np.array(["A", "A", "B", "A", "ctrl"])
    is_ctrl = np.array([False, True, False, False, True])
    # In test_perturbations={"A"} alone, row 1 would land in test. With the
    # mask, row 1 (which is ctrl) is excluded from test.
    train_idx, test_idx = make_train_test_split_by_perturbation(
        labels, test_perturbations={"A"}, is_ctrl_mask=is_ctrl
    )
    assert 1 not in test_idx.tolist()


def test_make_train_test_split_raises_when_no_test_rows():
    labels = np.array(["A", "B", "C"])
    with pytest.raises(ValueError, match="No test rows"):
        make_train_test_split_by_perturbation(labels, test_perturbations={"Z"})


def test_safe_variance_then_sample_round_trip():
    """End-to-end leakage-prevention path: split, train-variance, sample."""
    rng = np.random.default_rng(7)
    n_total, n_genes = 200, 8
    X = rng.standard_normal(size=(n_total, n_genes)).astype(np.float64)
    labels = np.array(["A"] * 50 + ["B"] * 50 + ["C"] * 50 + ["ctrl"] * 50)
    is_ctrl = np.array([L == "ctrl" for L in labels])

    train_idx, test_idx = make_train_test_split_by_perturbation(
        labels, test_perturbations={"C"}, is_ctrl_mask=is_ctrl
    )
    train_var = safe_train_variance(X, train_indices=train_idx,
                                    test_indices=test_idx)

    mu = X[train_idx].mean(axis=0)
    samples = sample_from_train_distribution(mu, train_var, n_cells=64, rng=rng)
    assert samples.shape == (64, n_genes)
