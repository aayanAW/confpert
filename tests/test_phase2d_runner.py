"""Smoke tests for the Phase 2D runner (verify + ablate + power + bootstrap)."""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from confpert import phase2d_runner


REPO = Path(__file__).resolve().parents[1]
PREREG_YAML = REPO / "paper_neurips_dnb_2026" / "preregistration_v2.yaml"
RESULTS_JSON = REPO / "baselines" / "results.json"


@pytest.mark.skipif(
    not (PREREG_YAML.exists() and RESULTS_JSON.exists()),
    reason="locked pre-reg yaml or Phase 1 results.json missing",
)
def test_phase2d_runner_emits_five_top_level_keys(tmp_path):
    out = tmp_path / "phase2d_report.json"
    report = phase2d_runner.run_phase2d(
        yaml_path=PREREG_YAML,
        results_path=RESULTS_JSON,
        repo_dir=REPO,
        out_path=out,
        skip_ablation=False,
    )
    assert set(report.keys()) >= {
        "meta", "verifier_report", "ablation_report",
        "power_report", "bootstrap_ci_table",
    }
    # JSON sidecar written
    assert out.exists()
    payload = json.loads(out.read_text())
    assert payload["verifier_report"]["lock_check_passed"] is True
    # Markdown sidecar
    md = out.with_suffix(".md")
    assert md.exists()
    assert "Phase 2D Analysis Report" in md.read_text()


@pytest.mark.skipif(
    not (PREREG_YAML.exists() and RESULTS_JSON.exists()),
    reason="locked pre-reg yaml or Phase 1 results.json missing",
)
def test_phase2d_runner_skip_ablation_flag(tmp_path):
    out = tmp_path / "phase2d_skip.json"
    report = phase2d_runner.run_phase2d(
        yaml_path=PREREG_YAML,
        results_path=RESULTS_JSON,
        repo_dir=REPO,
        out_path=out,
        skip_ablation=True,
    )
    assert report["ablation_report"] == {"skipped": True}


def test_phase2d_runner_bootstrap_uses_wilson_binomial():
    """Bootstrap table cells must report 'wilson_binomial' method when valid."""
    if not (PREREG_YAML.exists() and RESULTS_JSON.exists()):
        pytest.skip("requires locked pre-reg + results.json")
    report = phase2d_runner.run_phase2d(
        yaml_path=PREREG_YAML, results_path=RESULTS_JSON, repo_dir=REPO,
        out_path=None, skip_ablation=True,
    )
    bt = report["bootstrap_ci_table"]
    methods = {cell["ci"]["method"] for cell in bt.values()
               if cell["ci"] is not None}
    assert methods == {"wilson_binomial"}
