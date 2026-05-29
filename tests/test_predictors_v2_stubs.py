"""Unit tests for Phase 2 predictor metadata stubs."""
from __future__ import annotations

from confpert.predictors_v2_stubs import (
    PHASE_2_PREDICTORS,
    Phase2PredictorSpec,
    heavyweight_dataset_allowlist,
    is_heavyweight_predictor,
    list_phase2_predictors,
)


def test_all_phase2_predictors_match_prereg_yaml_ids():
    """The id field of each spec must match what preregistration_v2.yaml expects."""
    expected_ids = {"scgpt", "scfoundation", "geneformer", "cellfm", "get"}
    actual_ids = set(PHASE_2_PREDICTORS.keys())
    assert actual_ids == expected_ids


def test_phase2_predictors_have_required_metadata():
    for spec in PHASE_2_PREDICTORS.values():
        assert isinstance(spec, Phase2PredictorSpec)
        assert spec.id
        assert spec.family in {"F1_baselines", "F2_latent_delta", "F3_transformer", "EXCLUDED_H3"}
        assert spec.param_count > 0
        assert spec.corpus_diversity_index >= 0
        assert spec.license
        assert spec.upstream_repo.startswith("https://")
        assert spec.model_hub_path


def test_phase2_priority_order():
    """Must-have predictors first, deferred last."""
    priority = list_phase2_predictors(include_deferred=True)
    must_have_ids = ["scgpt", "scfoundation", "geneformer", "cellfm"]
    assert [p.id for p in priority[:4]] == must_have_ids
    assert priority[-1].id == "get"
    assert priority[-1].status == "deferred_optional"


def test_heavyweight_allowlist_matches_prereg_yaml():
    """Allowlist must match preregistration_v2.yaml §1.6."""
    expected = ["norman", "replogle_k562", "replogle_rpe1", "tahoe", "frangieh", "schmidt"]
    assert heavyweight_dataset_allowlist() == expected


def test_is_heavyweight_predictor():
    assert is_heavyweight_predictor("scgpt")
    assert is_heavyweight_predictor("scfoundation")
    assert is_heavyweight_predictor("geneformer")
    assert is_heavyweight_predictor("cellfm")
    assert is_heavyweight_predictor("state")  # Phase 1 also heavyweight
    assert not is_heavyweight_predictor("mean")
    assert not is_heavyweight_predictor("ahlmann_bilinear_ridge")
    assert not is_heavyweight_predictor("noisy_mean")


def test_excluded_from_h3_consistency():
    """GET is excluded from H3 (per prereg v2.4); same for GEARS but that's not Phase 2."""
    get_spec = PHASE_2_PREDICTORS["get"]
    assert get_spec.family == "EXCLUDED_H3"
