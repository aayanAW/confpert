"""Tests for confpert.bootstrap stratified-percentile CIs."""
from __future__ import annotations

import numpy as np
import pytest

from confpert.bootstrap import (
    BootstrapCI,
    bootstrap_aggregate_per_dataset,
    bootstrap_coverage_from_scores,
    bootstrap_deviation_from_scores_proper,
)


def test_coverage_ci_brackets_observed_for_well_calibrated_synthetic():
    """For perfectly-calibrated synthetic scores (calib ~ test), the
    observed coverage should be inside the 95% CI."""
    rng = np.random.default_rng(0)
    calib = rng.standard_normal(100)
    test = rng.standard_normal(100)
    ci = bootstrap_coverage_from_scores(calib, test, alpha=0.10, n_resamples=500, seed=42)
    assert ci.lo <= ci.point <= ci.hi
    # Observed coverage should be roughly 1-alpha = 0.90 for well-calibrated synthetic
    assert 0.7 <= ci.point <= 1.0


def test_coverage_ci_shrinks_with_more_samples():
    """Bigger n => narrower CI."""
    rng = np.random.default_rng(1)
    small = bootstrap_coverage_from_scores(
        rng.standard_normal(20), rng.standard_normal(20),
        alpha=0.10, n_resamples=500, seed=42,
    )
    rng2 = np.random.default_rng(2)
    big = bootstrap_coverage_from_scores(
        rng2.standard_normal(500), rng2.standard_normal(500),
        alpha=0.10, n_resamples=500, seed=42,
    )
    small_width = small.hi - small.lo
    big_width = big.hi - big.lo
    assert big_width < small_width


def test_coverage_ci_deterministic_under_fixed_seed():
    """Same inputs + seed => same output."""
    rng = np.random.default_rng(3)
    c1 = rng.standard_normal(50)
    t1 = rng.standard_normal(50)
    ci_a = bootstrap_coverage_from_scores(c1, t1, alpha=0.05, n_resamples=200, seed=99)
    ci_b = bootstrap_coverage_from_scores(c1, t1, alpha=0.05, n_resamples=200, seed=99)
    assert ci_a.point == ci_b.point
    assert ci_a.lo == ci_b.lo
    assert ci_a.hi == ci_b.hi


def test_deviation_ci_brackets_observed():
    rng = np.random.default_rng(4)
    calib = rng.standard_normal(60)
    test = rng.standard_normal(60)
    ci = bootstrap_deviation_from_scores_proper(calib, test, alpha=0.20, n_resamples=400, seed=42)
    assert ci.lo <= ci.point <= ci.hi
    assert ci.point >= 0
    assert ci.method == "stratified_percentile_deviation"


def test_deviation_ci_returns_nonneg_bounds():
    rng = np.random.default_rng(5)
    calib = rng.standard_normal(40)
    test = rng.standard_normal(40) + 2.0   # shift test scores to inflate deviation
    ci = bootstrap_deviation_from_scores_proper(calib, test, alpha=0.10, n_resamples=300, seed=42)
    assert ci.lo >= 0
    assert ci.hi >= 0
    assert ci.point > 0


def test_aggregate_per_dataset_basic():
    devs = {"mean": 0.10, "ridge": 0.05, "scgen": 0.08, "cpa": 0.15, "biolord": 0.07}
    point, lo, hi = bootstrap_aggregate_per_dataset(devs, n_resamples=500, seed=42)
    assert lo <= point <= hi
    assert abs(point - np.mean(list(devs.values()))) < 1e-10


def test_aggregate_per_dataset_raises_on_empty():
    with pytest.raises(ValueError):
        bootstrap_aggregate_per_dataset({})


def test_bootstrap_ci_dataclass_fields():
    rng = np.random.default_rng(6)
    ci = bootstrap_coverage_from_scores(
        rng.standard_normal(30), rng.standard_normal(30),
        alpha=0.10, n_resamples=100, seed=42,
    )
    assert isinstance(ci, BootstrapCI)
    assert ci.n_resamples == 100
    assert ci.n_perts_calib == 30
    assert ci.n_perts_test == 30
    assert ci.method == "stratified_percentile"
