"""Tests for confpert.power Phase 2 power-analysis report."""
from __future__ import annotations

import pytest

from confpert.power import (
    PowerResult,
    full_power_report,
    k1_v2_h1_h1b_power_report,
    k2_v2_power_report,
    power_kruskal_wallis,
    power_spearman_per_dataset,
    power_two_way_anova,
)


def test_k2_v2_report_keys():
    rep = k2_v2_power_report()
    assert set(rep.keys()) == {"H2_anova", "H2b_spearman_per_dataset", "H3_kruskal_wallis"}
    for v in rep.values():
        assert isinstance(v, PowerResult)


def test_k1_v2_h1_h1b_report_keys():
    rep = k1_v2_h1_h1b_power_report()
    assert set(rep.keys()) == {"H1_capacity_per_dataset_phase2", "H1b_data_cross_dataset_phase2"}


def test_full_power_report_combines_both():
    full = full_power_report()
    assert "H1_capacity_per_dataset_phase2" in full
    assert "H2_anova" in full
    assert "H3_kruskal_wallis" in full


def test_h2_anova_is_well_powered():
    """H2 ANOVA at n=120 cells, eta^2=0.10 target, alpha=0.0167 should be
    well-powered (≥0.80)."""
    r = power_two_way_anova(n_cells=120, n_factors=2, target_eta_sq=0.10, alpha=0.0167)
    assert r.estimated_power >= 0.80, f"H2 underpowered: {r.estimated_power}"


def test_h3_kw_is_well_powered_at_delta_0_5():
    """Pre-reg v2.0 final threshold delta=0.50 should be ≥0.80 powered at n=40/group."""
    r = power_kruskal_wallis(
        n_per_group=40, k_groups=3, target_cliffs_delta=0.50, alpha=0.0167,
        n_simulations=500,
    )
    assert r.estimated_power >= 0.75, f"H3 underpowered: {r.estimated_power}"


def test_spearman_smaller_n_lower_power():
    big = power_spearman_per_dataset(
        n_predictors=20, target_rho=0.5, alpha=0.05, n_simulations=500,
    )
    small = power_spearman_per_dataset(
        n_predictors=8, target_rho=0.5, alpha=0.05, n_simulations=500,
    )
    assert big.estimated_power > small.estimated_power


def test_power_result_underpowered_flag_in_notes():
    """Spearman with very small n at strict alpha should flag UNDERPOWERED in notes."""
    r = power_spearman_per_dataset(
        n_predictors=5, target_rho=0.5, alpha=0.0167, n_simulations=300,
    )
    if r.estimated_power < 0.80:
        assert "UNDERPOWERED" in r.notes
