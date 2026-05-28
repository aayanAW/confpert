"""Smoke tests for Phase 2 dataset loader stubs."""
from __future__ import annotations

import pytest

from confpert.data._phase2_stubs import (
    PHASE_2_DATASETS,
    Phase2DatasetSpec,
    list_phase2_datasets,
    load_frangieh,
    load_schmidt,
    load_datlinger,
    load_mcfaline_figueroa,
    load_walker,
    load_lara_astiaso,
)


def test_all_phase2_dataset_ids_match_prereg_yaml():
    expected = {"frangieh", "schmidt", "datlinger", "mcfaline_figueroa",
                "lara_astiaso", "walker"}
    assert set(PHASE_2_DATASETS.keys()) == expected


def test_phase2_dataset_specs_have_required_metadata():
    for spec in PHASE_2_DATASETS.values():
        assert isinstance(spec, Phase2DatasetSpec)
        assert spec.id
        assert spec.cell_line_context in {"K562", "non-K562", "primary"}
        assert spec.organism in {"human", "mouse"}
        assert spec.perturbation_type in {"genetic", "chemical"}
        assert spec.expected_n_cells > 0
        assert spec.expected_n_perturbations > 0
        assert spec.public_url.startswith("https://")
        assert spec.citation


def test_mouse_organism_dataset_present_for_cross_organism_split():
    """lara_astiaso replaces walker as the mouse-organism dataset per
    pre-reg amendment 2026-05-25 (option B).
    """
    assert PHASE_2_DATASETS["lara_astiaso"].organism == "mouse"
    active_human = {
        spec.organism for k, spec in PHASE_2_DATASETS.items()
        if k not in {"lara_astiaso", "walker"}
    }
    assert active_human == {"human"}


def test_mcfaline_is_chemical_for_held_out_perturbation_type_split():
    """McFaline-Figueroa is one chemical dataset; needed for held-out-perturbation-type split."""
    assert PHASE_2_DATASETS["mcfaline_figueroa"].perturbation_type == "chemical"


@pytest.mark.parametrize("loader,expected_exc", [
    # Real Phase 2C loaders raise FileNotFoundError when h5ad missing locally
    (load_frangieh, FileNotFoundError),
    (load_schmidt, FileNotFoundError),
    (load_datlinger, FileNotFoundError),
    (load_mcfaline_figueroa, FileNotFoundError),
    (load_lara_astiaso, FileNotFoundError),
    # Walker remains stubbed (superseded by lara_astiaso)
    (load_walker, NotImplementedError),
])
def test_loaders_raise_with_helpful_message(loader, expected_exc):
    """Real loaders raise FileNotFoundError + manual-download URL.
    Walker raises NotImplementedError with A/B/C resolution path.
    """
    with pytest.raises(expected_exc) as exc:
        loader()
    msg = str(exc.value)
    assert "https://" in msg or "URL" in msg or "Phase 2" in msg


def test_loader_status_correctly_marked():
    """5 implemented (incl. lara_astiaso), walker superseded."""
    assert PHASE_2_DATASETS["frangieh"].loader_status == "implemented"
    assert PHASE_2_DATASETS["schmidt"].loader_status == "implemented"
    assert PHASE_2_DATASETS["datlinger"].loader_status == "implemented"
    assert PHASE_2_DATASETS["mcfaline_figueroa"].loader_status == "implemented"
    assert PHASE_2_DATASETS["lara_astiaso"].loader_status == "implemented"
    assert PHASE_2_DATASETS["walker"].loader_status == "superseded"


def test_walker_error_documents_three_resolution_options():
    """Walker's NotImplementedError MUST surface the A/B/C resolution path."""
    with pytest.raises(NotImplementedError) as exc:
        load_walker()
    msg = str(exc.value)
    assert "A." in msg or "Confirm" in msg
    assert "B." in msg or "Amend" in msg or "substitution" in msg.lower()
    assert "C." in msg or "Drop" in msg


def test_resolved_urls_point_to_real_https():
    """Phase 2C URL resolution: all active public_url fields are real, not 'TBD'."""
    for ds_id in ("frangieh", "schmidt", "datlinger", "mcfaline_figueroa",
                  "lara_astiaso"):
        url = PHASE_2_DATASETS[ds_id].public_url
        assert url.startswith("https://"), f"{ds_id} public_url not https: {url}"
        assert "TBD" not in url, f"{ds_id} still has TBD URL: {url}"
        assert "acc.cgi" not in url, (
            f"{ds_id} still pointing at GEO acc.cgi placeholder, not a "
            f"resolved direct download: {url}"
        )


def test_priority_order():
    """Stub list order matches priority: Frangieh first (multi-modal RNA + protein)."""
    priority = list_phase2_datasets()
    assert priority[0].id == "frangieh"
    assert priority[1].id == "schmidt"
