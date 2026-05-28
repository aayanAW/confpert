"""Unit tests for confpert.prereg (Phase 2A.3 + 2A.5)."""
from __future__ import annotations

import json
import tempfile
from pathlib import Path

import pytest

from confpert import prereg


def _minimal_yaml_spec() -> dict:
    return {
        "version": "1.0",
        "project": "confpert-test",
        "phase": 2,
        "companion_md": "companion.md",
        "lock": {
            "git_commit_md": None,
            "git_commit_yaml": None,
            "sha256_md": None,
            "sha256_yaml": None,
            "osf_doi": None,
            "results_json_max_mtime": None,
            "random_seed": 42,
        },
        "k1_v2": {
            "predictors": [{"id": "x", "family": "F1_baselines"}],
            "datasets": [{"id": "d1", "cell_line_context": "K562"}],
            "splits": [{"id": "within_perturbation"}],
            "discrepancies": ["ks"],
            "alphas": [0.05],
        },
        "k2_v2": {
            "family_alpha": 0.05,
            "bonferroni_n": 3,
            "per_hypothesis_alpha": 0.0167,
            "hypotheses": {
                "H2_cell_line_covariate": {
                    "family_member": True,
                    "test": "two_way_anova_type2",
                    "success_criteria": [
                        {"test": "f_test_pvalue", "operator": "<", "threshold": 0.0167}
                    ],
                    "disposition_pass": "p",
                    "disposition_fail": "f",
                },
            },
            "outcome_table": {"H2=P,H2b=P,H3=P": "triple_pass"},
        },
    }


def _write_yaml(tmpdir: Path, spec: dict) -> Path:
    import yaml
    p = tmpdir / "preregistration_v2.yaml"
    with open(p, "w") as fh:
        yaml.safe_dump(spec, fh, sort_keys=False, default_flow_style=False)
    return p


def _write_results(tmpdir: Path, n_rows: int = 5) -> Path:
    results = {
        "rows": [
            {
                "predictor": "mean",
                "dataset": "norman",
                "alpha": 0.05,
                "scores": {"ks": {"achieved_coverage": 0.9}},
            }
            for _ in range(n_rows)
        ],
        "sha256": "test",
    }
    p = tmpdir / "results.json"
    with open(p, "w") as fh:
        json.dump(results, fh)
    return p


def test_sha256_of_file_deterministic():
    with tempfile.TemporaryDirectory() as td:
        p = Path(td) / "x.txt"
        p.write_text("hello world")
        h1 = prereg._sha256_of_file(p)
        h2 = prereg._sha256_of_file(p)
        assert h1 == h2
        assert len(h1) == 64


def test_validate_schema_passes_on_minimal_spec():
    errors = prereg._validate_schema(_minimal_yaml_spec())
    assert errors == []


def test_validate_schema_catches_missing_keys():
    spec = _minimal_yaml_spec()
    del spec["k2_v2"]
    errors = prereg._validate_schema(spec)
    assert any("k2_v2" in e for e in errors)


def test_validate_schema_catches_missing_hypothesis_criteria():
    spec = _minimal_yaml_spec()
    spec["k2_v2"]["hypotheses"]["H2_cell_line_covariate"].pop("success_criteria")
    errors = prereg._validate_schema(spec)
    assert any("success_criteria" in e for e in errors)


def test_emit_hashes_round_trip():
    """emit_hashes writes SHA-256; subsequent verify confirms self-hash."""
    with tempfile.TemporaryDirectory() as td:
        tmpdir = Path(td)
        spec = _minimal_yaml_spec()
        yaml_path = _write_yaml(tmpdir, spec)
        # Also create companion .md so emit-hashes can hash it
        (tmpdir / "companion.md").write_text("# minimal companion\n")
        prereg.emit_hashes(yaml_path)
        import yaml
        with open(yaml_path) as fh:
            stamped = yaml.safe_load(fh)
        assert stamped["lock"]["sha256_yaml"] is not None
        assert len(stamped["lock"]["sha256_yaml"]) == 64
        assert stamped["lock"]["sha256_md"] is not None


def test_verify_dry_run_passes_on_locked_yaml():
    """Lock check should pass after emit_hashes."""
    with tempfile.TemporaryDirectory() as td:
        tmpdir = Path(td)
        spec = _minimal_yaml_spec()
        yaml_path = _write_yaml(tmpdir, spec)
        (tmpdir / "companion.md").write_text("# minimal\n")
        prereg.emit_hashes(yaml_path)
        results_path = _write_results(tmpdir, n_rows=3)
        result = prereg.verify(yaml_path, results_path, repo_dir=tmpdir, dry_run=True)
        assert result.dry_run
        # No errors (the only check that requires git is mtime, which is skipped without commit)
        assert all("yaml sha256 mismatch" not in e for e in result.errors)
        # Hypothesis rows enumerated
        assert len(result.hypothesis_results) == 1
        assert result.hypothesis_results[0].status == "DRY_RUN"


def test_verify_dry_run_fails_if_yaml_modified_after_lock():
    """Sha256 mismatch should be flagged if yaml is tampered."""
    with tempfile.TemporaryDirectory() as td:
        tmpdir = Path(td)
        spec = _minimal_yaml_spec()
        yaml_path = _write_yaml(tmpdir, spec)
        (tmpdir / "companion.md").write_text("# minimal\n")
        prereg.emit_hashes(yaml_path)
        results_path = _write_results(tmpdir)
        # Tamper: edit a field that isn't sha256_yaml itself
        import yaml as _y
        with open(yaml_path) as fh:
            tampered = _y.safe_load(fh)
        tampered["project"] = "TAMPERED"
        with open(yaml_path, "w") as fh:
            _y.safe_dump(tampered, fh, sort_keys=False, default_flow_style=False)
        # Verify with dry_run=False (where hash check is enforced)
        result = prereg.verify(yaml_path, results_path, repo_dir=tmpdir, dry_run=False)
        assert any("yaml sha256 mismatch" in e for e in result.errors)


def test_disposition_table_correctness():
    """Outcome table lookup returns the right disposition string."""
    spec = _minimal_yaml_spec()
    spec["k2_v2"]["outcome_table"] = {
        "H2=P,H2b=P,H3=P": "triple_pass",
        "H2=F,H2b=F,H3=F": "triple_null",
    }
    # All three pass
    hr_pass = [
        prereg.VerificationResult(hypothesis_id="H2_cell_line_covariate", status="PASS"),
        prereg.VerificationResult(hypothesis_id="H2b_corpus_diversity", status="PASS"),
        prereg.VerificationResult(hypothesis_id="H3_architecture_family", status="PASS"),
    ]
    assert prereg._apply_disposition_table(spec["k2_v2"], hr_pass) == "triple_pass"
    # All three fail
    hr_fail = [
        prereg.VerificationResult(hypothesis_id=h, status="FAIL")
        for h in ["H2_cell_line_covariate", "H2b_corpus_diversity", "H3_architecture_family"]
    ]
    assert prereg._apply_disposition_table(spec["k2_v2"], hr_fail) == "triple_null"


def test_disposition_returns_none_on_incomplete_family():
    """If any K2-v2-family hypothesis is DRY_RUN, no disposition is computed."""
    spec = _minimal_yaml_spec()
    hr_partial = [
        prereg.VerificationResult(hypothesis_id="H2_cell_line_covariate", status="PASS"),
        prereg.VerificationResult(hypothesis_id="H2b_corpus_diversity", status="DRY_RUN"),
        prereg.VerificationResult(hypothesis_id="H3_architecture_family", status="PASS"),
    ]
    assert prereg._apply_disposition_table(spec["k2_v2"], hr_partial) is None
