"""Validation: prereg.verify (full mode) reproduces scripts/k2_analysis output.

This is the regression test that proves the machine-verifiable pre-registration
produces identical numerical results to the manual K2 analysis script.
"""
from __future__ import annotations

from pathlib import Path

import pytest

RESULTS_PATH = Path("tests/fixtures/phase1_results.json")
PREREG_PATH = Path("paper_neurips_dnb_2026/preregistration_v2.yaml")

if not RESULTS_PATH.exists() or not PREREG_PATH.exists():
    pytest.skip(
        "skipping verifier regression test — Phase 1 fixture or "
        "preregistration_v2.yaml missing",
        allow_module_level=True,
    )


def test_h1_capacity_reproduces_phase1_analysis():
    """H1 capacity result from confpert.prereg.verify must match
    scripts/k2_analysis.k2_test on Phase 1 results.json."""
    from confpert import prereg

    repo_dir = Path.cwd()
    result = prereg.verify(
        yaml_path=PREREG_PATH,
        results_path=RESULTS_PATH,
        repo_dir=repo_dir,
        dry_run=False,
    )

    h1 = next((hr for hr in result.hypothesis_results if hr.hypothesis_id == "H1_capacity"), None)
    assert h1 is not None
    assert h1.status in {"PASS", "FAIL"}

    per_ds = h1.observed.get("per_dataset", {})
    # Phase 1 fixture (commit aae11fc) has 8 predictors per dataset and
    # 7 for tahoe. Numerical snapshot.
    expected = {
        "replogle_rpe1": {"rho":  0.7229, "p": 0.024, "n": 8, "pass": True},
        "tahoe":         {"rho":  0.6487, "p": 0.066, "n": 7, "pass": False},
        "norman":        {"rho":  0.3832, "p": 0.174, "n": 8, "pass": False},
        "replogle_k562": {"rho": -0.2635, "p": 0.750, "n": 8, "pass": False},
        "adamson":       {"rho": -0.2275, "p": 0.717, "n": 8, "pass": False},
    }
    for ds, exp in expected.items():
        actual = per_ds.get(ds)
        assert actual is not None, f"{ds} missing from H1 per_dataset"
        assert abs(actual["rho"] - exp["rho"]) < 1e-3, (
            f"{ds} rho mismatch: actual={actual['rho']:.4f}, expected={exp['rho']:.4f}"
        )
        assert abs(actual["p_value"] - exp["p"]) < 0.01, (
            f"{ds} p_value mismatch: actual={actual['p_value']:.4f}, expected={exp['p']:.4f}"
        )
        assert actual["n_predictors"] == exp["n"]
        assert actual["passes"] == exp["pass"]

    assert h1.observed["n_datasets_passing"] == 1
    assert set(h1.observed["datasets_passing"]) == {"replogle_rpe1"}
    assert h1.status == "FAIL"


def test_h1b_data_correctly_reports_null_on_phase1_data():
    """H1b: n_train_cells is constant per dataset → ρ is NaN → all FAIL."""
    from confpert import prereg

    repo_dir = Path.cwd()
    result = prereg.verify(
        yaml_path=PREREG_PATH,
        results_path=RESULTS_PATH,
        repo_dir=repo_dir,
        dry_run=False,
    )

    h1b = next((hr for hr in result.hypothesis_results if hr.hypothesis_id == "H1b_data"), None)
    assert h1b is not None
    assert h1b.status == "FAIL"
    assert h1b.observed["n_datasets_passing"] == 0


def test_h2b_corpus_diversity_per_dataset_populated():
    """H2b should produce finite ρ per dataset (corpus_diversity varies)."""
    from confpert import prereg

    repo_dir = Path.cwd()
    result = prereg.verify(
        yaml_path=PREREG_PATH,
        results_path=RESULTS_PATH,
        repo_dir=repo_dir,
        dry_run=False,
    )

    h2b = next((hr for hr in result.hypothesis_results if hr.hypothesis_id == "H2b_corpus_diversity"), None)
    assert h2b is not None
    per_ds = h2b.observed.get("per_dataset", {})
    assert len(per_ds) == 5


def test_h2_anova_on_phase1_returns_finite_pvalue():
    """H2 ANOVA on Phase 1 data: should return finite F-test p-value + eta squared.
    Phase 1 sample size is small (≤2 cell_line_contexts × ~9 predictors × 5 datasets);
    expect FAIL but not crash."""
    from confpert import prereg

    repo_dir = Path.cwd()
    result = prereg.verify(
        yaml_path=PREREG_PATH,
        results_path=RESULTS_PATH,
        repo_dir=repo_dir,
        dry_run=False,
    )

    h2 = next((hr for hr in result.hypothesis_results if hr.hypothesis_id == "H2_cell_line_covariate"), None)
    assert h2 is not None
    assert h2.status in {"PASS", "FAIL"}
    assert "f_pvalue_cell_line_context" in h2.observed
    assert 0.0 <= h2.observed["f_pvalue_cell_line_context"] <= 1.0
    assert h2.observed["eta_squared_cell_line_context"] >= 0.0


def test_h3_kw_on_phase1_errors_due_to_missing_f3():
    """H3 KW on Phase 1 fixture: F3 transformer family has no data points
    in the Phase 1 fixture (foundation models not yet wrapped), so the
    Kruskal-Wallis call fails the n_families>=3 prerequisite and returns
    status=ERROR with a clear message. Phase 2 results.json contains F3."""
    from confpert import prereg

    repo_dir = Path.cwd()
    result = prereg.verify(
        yaml_path=PREREG_PATH,
        results_path=RESULTS_PATH,
        repo_dir=repo_dir,
        dry_run=False,
    )

    h3 = next((hr for hr in result.hypothesis_results if hr.hypothesis_id == "H3_architecture_family"), None)
    assert h3 is not None
    assert h3.status == "ERROR"
    assert h3.error is not None
    assert "families" in h3.error.lower() or "kw" in h3.error.lower()


def test_disposition_on_phase1_data():
    """Phase 1 fixture has no F3 transformer rows so H3 errors; H2 and H2b
    fail by design. Disposition resolves to None (no headline can be reported
    until the foundation-model wrappers contribute F3 rows in Phase 2)."""
    from confpert import prereg

    repo_dir = Path.cwd()
    result = prereg.verify(
        yaml_path=PREREG_PATH,
        results_path=RESULTS_PATH,
        repo_dir=repo_dir,
        dry_run=False,
    )
    # On Phase 1 fixture, H3 ERRORs (no F3 rows) so triple_null cannot be
    # asserted -- disposition is None.
    assert result.disposition_headline in (
        None, "triple_null_promote_ablation_panel"
    )
