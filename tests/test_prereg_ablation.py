"""Tests for confpert.prereg_ablation (Pipeline A vs B counterfactual)."""
from __future__ import annotations

from pathlib import Path

import pytest

from confpert import prereg_ablation

RESULTS_PATH = Path("tests/fixtures/phase1_results.json")
PREREG_PATH = Path("paper_neurips_dnb_2026/preregistration_v2.yaml")

if not RESULTS_PATH.exists() or not PREREG_PATH.exists():
    pytest.skip(
        "skipping ablation tests — Phase 1 fixture or "
        "preregistration_v2.yaml missing",
        allow_module_level=True,
    )


def test_ablation_panel_runs_end_to_end():
    """Smoke: run panel on Phase 1 data, no exceptions."""
    rep = prereg_ablation.run_ablation_panel(
        yaml_path=PREREG_PATH, results_path=RESULTS_PATH, repo_dir=Path.cwd(),
    )
    assert isinstance(rep, prereg_ablation.AblationReport)
    assert len(rep.rows) == 5  # H1, H1b, H2, H2b, H3


def test_h1_flips_to_pass_via_p_value_relaxation_on_phase1_data():
    """The key methodological demonstration: H1 (capacity) on Phase 1 data
    flips FAIL → PASS when p threshold relaxed from 0.05 to 0.10.

    This is the empirical evidence for the pre-reg-ablation Tier-1
    contribution per pre-reg v2 §4. Norman's p=0.057 is the load-bearing
    data point — it's just above 0.05 pre-reg threshold, but a counterfactual
    analyst could have shopped to p<0.10 to flip the headline."""
    rep = prereg_ablation.run_ablation_panel(
        yaml_path=PREREG_PATH, results_path=RESULTS_PATH, repo_dir=Path.cwd(),
    )
    h1 = next((r for r in rep.rows if r.hypothesis == "H1_capacity"), None)
    assert h1 is not None
    assert h1.prereg_result == "FAIL"
    assert h1.flipped_to_pass, f"H1 should flip via shopped relaxation; got: {h1.shopped_threshold}"
    # On the Phase 1 dbcc769 fixture, the flip path is the combined
    # "p<0.10 AND n_passing>=2" relaxation: Tahoe p=0.066 + Replogle RPE1
    # p=0.024 together meet 2-of-4 once p threshold loosens.
    assert "p<0.10" in h1.shopped_threshold


def test_h2_h2b_h3_remain_fail_under_shopped_thresholds():
    """The other hypotheses (H2, H2b, H3) should NOT flip on Phase 1 data
    even with shopped thresholds — Phase 1 data effects are weak enough that
    no plausible relaxation rescues them. This strengthens the pre-reg
    discipline argument."""
    rep = prereg_ablation.run_ablation_panel(
        yaml_path=PREREG_PATH, results_path=RESULTS_PATH, repo_dir=Path.cwd(),
    )
    for hyp_id in ["H1b_data", "H2_cell_line_covariate", "H2b_corpus_diversity",
                   "H3_architecture_family"]:
        row = next((r for r in rep.rows if r.hypothesis == hyp_id), None)
        assert row is not None, f"{hyp_id} missing from ablation panel"
        assert not row.flipped_to_pass, (
            f"{hyp_id} unexpectedly flipped to PASS via {row.shopped_threshold}"
        )


def test_ablation_headline_reports_correct_flip_count():
    rep = prereg_ablation.run_ablation_panel(
        yaml_path=PREREG_PATH, results_path=RESULTS_PATH, repo_dir=Path.cwd(),
    )
    assert rep.n_flipped == 1
    assert "1 of 5" in rep.headline


def test_shopped_relaxations_table_well_formed():
    """SHOPPED_RELAXATIONS dict must cover the test types our pre-reg uses."""
    expected_test_types = {"spearman_per_dataset", "two_way_anova_type2", "kruskal_wallis"}
    assert set(prereg_ablation.SHOPPED_RELAXATIONS.keys()) >= expected_test_types
    for test_type, relaxations in prereg_ablation.SHOPPED_RELAXATIONS.items():
        assert len(relaxations) >= 2, f"{test_type} needs at least 2 relaxations"
        for r in relaxations:
            assert "name" in r, f"{test_type} relaxation missing 'name': {r}"


def test_row_format():
    """AblationRow dataclass must serialise to JSON cleanly (used in --out flag)."""
    from dataclasses import asdict
    rep = prereg_ablation.run_ablation_panel(
        yaml_path=PREREG_PATH, results_path=RESULTS_PATH, repo_dir=Path.cwd(),
    )
    for row in rep.rows:
        d = asdict(row)
        assert "hypothesis" in d
        assert "prereg_threshold" in d
        assert "shopped_threshold" in d
        assert "flipped_to_pass" in d
