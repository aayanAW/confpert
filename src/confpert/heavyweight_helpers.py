"""Shared utility helpers for heavyweight predictors (scGPT, scFoundation,
Geneformer, CellFM, STATE).

Every heavyweight wrapper has the same data-leakage failure mode:
estimating per-gene sampling variance from the TEST fold and using it to
generate "predicted" samples. The scGPT v1 wrapper hit this exact bug
2026-05-25 and was reverted before any launch. To prevent recurrence:

  - Use `safe_train_variance(...)` to compute per-gene variance with an
    explicit train-fold indices argument. The function refuses any input
    where the index set overlaps with the test indices.
  - Use `sample_from_train_distribution(mu_per_pert, train_variance, n_cells)`
    to generate samples from a fitted train-fold mean + variance. The
    sampler refuses to accept a variance vector whose source is None
    (forcing every caller to pass the train indices).

Anything that touches per-gene sampling noise in the Modal wrappers must
flow through these helpers — there is no other validated path.
"""
from __future__ import annotations

from typing import Iterable

import numpy as np


def safe_train_variance(
    X: np.ndarray,
    train_indices: np.ndarray | Iterable[int],
    test_indices: np.ndarray | Iterable[int] | None = None,
    min_train_n: int = 30,
    eps: float = 1e-6,
) -> np.ndarray:
    """Compute per-gene variance using TRAIN-fold rows ONLY.

    Refuses any input where `train_indices` and `test_indices` overlap or
    where `train_indices` is shorter than `min_train_n`. Returns a
    `(n_genes,)` float64 vector of per-gene variances clipped at `eps`.

    The arguments are intentionally explicit: the caller MUST pass the
    fold partition that they intend to use. We do not accept a single
    array + an implicit "use the first 80%" rule; that pattern is exactly
    how the scGPT v1 leak slipped past review.

    Parameters
    ----------
    X : np.ndarray
        Full data matrix, shape (n_total_cells, n_genes).
    train_indices : array-like of int
        Row indices in X that constitute the training fold.
    test_indices : array-like of int or None
        Row indices in X that constitute the test fold. If provided,
        an overlap with `train_indices` raises ValueError.
    min_train_n : int
        Minimum train-fold size. Smaller folds give pathological variance
        estimates and signal a likely upstream bug.
    eps : float
        Lower clip for per-gene variance, to avoid zero-variance genes
        causing degenerate sampler scale.
    """
    train_idx = np.asarray(list(train_indices), dtype=np.int64)
    if train_idx.ndim != 1:
        raise ValueError(f"train_indices must be 1D, got shape {train_idx.shape}")

    # Range-check before size-check so out-of-range indices surface as
    # IndexError rather than masked behind a "not enough rows" ValueError.
    n_total = X.shape[0]
    if train_idx.size > 0 and (train_idx.max() >= n_total or train_idx.min() < 0):
        raise IndexError(
            f"train_indices out of range [0, {n_total}); got "
            f"min={train_idx.min()}, max={train_idx.max()}"
        )

    if train_idx.size < min_train_n:
        raise ValueError(
            f"train_indices has {train_idx.size} rows < min_train_n={min_train_n}. "
            f"Variance estimate would be too noisy; this is likely an upstream bug "
            f"in the train/test split."
        )

    if test_indices is not None:
        test_idx = np.asarray(list(test_indices), dtype=np.int64)
        overlap = np.intersect1d(train_idx, test_idx, assume_unique=False)
        if overlap.size > 0:
            raise ValueError(
                f"train_indices ∩ test_indices = {overlap.size} overlapping rows. "
                f"This would leak test-fold information into the sampling variance. "
                f"First 5 overlapping rows: {overlap[:5].tolist()}"
            )

    X_train = X[train_idx]
    # Use ddof=1 (unbiased) since the train fold is typically large enough.
    var = X_train.var(axis=0, ddof=1).astype(np.float64)
    var = np.clip(var, eps, None)
    return var


def sample_from_train_distribution(
    mu_per_pert: np.ndarray,
    train_variance: np.ndarray,
    n_cells: int,
    rng: np.random.Generator | None = None,
) -> np.ndarray:
    """Generate `n_cells` samples ~ N(mu_per_pert, diag(train_variance)).

    `train_variance` MUST come from `safe_train_variance`. We don't have a
    cryptographic stamp but we can validate length + dtype + nonneg.
    """
    if mu_per_pert.ndim != 1:
        raise ValueError(f"mu_per_pert must be 1D (n_genes), got shape "
                         f"{mu_per_pert.shape}")
    if train_variance.ndim != 1 or train_variance.shape != mu_per_pert.shape:
        raise ValueError(
            f"train_variance must match mu_per_pert shape "
            f"{mu_per_pert.shape}, got {train_variance.shape}"
        )
    if (train_variance < 0).any():
        raise ValueError("train_variance has negative entries; not a valid variance")
    if n_cells <= 0:
        raise ValueError(f"n_cells must be > 0, got {n_cells}")

    if rng is None:
        rng = np.random.default_rng(42)

    std = np.sqrt(train_variance)
    eps = rng.standard_normal(size=(n_cells, mu_per_pert.size)).astype(np.float64)
    return mu_per_pert[None, :] + eps * std[None, :]


def make_train_test_split_by_perturbation(
    labels: np.ndarray,
    test_perturbations: set[str] | list[str] | tuple[str, ...],
    is_ctrl_mask: np.ndarray | None = None,
    rng: np.random.Generator | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Build train / test row-index partitions for a per-perturbation split.

    Cells whose perturbation label is in `test_perturbations` go to test;
    every other cell (including controls, unless excluded via mask) goes
    to train. Control cells should usually appear in BOTH folds (used for
    predicting `X_pert | X_ctrl`); to enforce that, omit `is_ctrl_mask`.
    To put controls only in train, pass the is_ctrl_mask.

    Returned indices are sorted ascending.
    """
    test_set = set(test_perturbations)
    labels_arr = np.asarray(labels)
    if labels_arr.ndim != 1:
        raise ValueError(f"labels must be 1D, got shape {labels_arr.shape}")

    is_test = np.array([str(L) in test_set for L in labels_arr], dtype=bool)
    test_idx = np.where(is_test)[0]
    train_idx = np.where(~is_test)[0]

    if is_ctrl_mask is not None:
        ctrl_mask = np.asarray(is_ctrl_mask, dtype=bool)
        if ctrl_mask.shape != labels_arr.shape:
            raise ValueError(
                f"is_ctrl_mask shape {ctrl_mask.shape} != labels shape "
                f"{labels_arr.shape}"
            )
        # Drop control cells from test (controls only in train)
        ctrl_in_test = ctrl_mask[test_idx]
        test_idx = test_idx[~ctrl_in_test]

    if test_idx.size == 0:
        raise ValueError(
            f"No test rows. test_perturbations={test_perturbations!r} did not "
            f"match any label. First 5 unique labels: "
            f"{list(set(labels_arr.tolist()))[:5]}"
        )
    if train_idx.size == 0:
        raise ValueError("No train rows.")

    overlap = np.intersect1d(train_idx, test_idx)
    if overlap.size > 0:  # pragma: no cover -- partition invariant
        raise AssertionError(
            f"train ∩ test = {overlap.size}; partition logic bug"
        )

    return np.sort(train_idx), np.sort(test_idx)
