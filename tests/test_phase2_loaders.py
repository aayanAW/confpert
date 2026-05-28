"""Smoke + integration tests for Phase 2 dataset loaders.

These tests verify:
  (a) loaders import without errors
  (b) loaders surface a clear FileNotFoundError when local h5ad missing
  (c) loaders accept the documented kwargs (signature stable for Modal runner)
  (d) loaders work end-to-end on a small synthetic AnnData when one is staged
      at a known cache path (skipped if not present locally)

We do NOT download real datasets in unit tests (12 GB McFaline tar would
break CI). Real-data integration tests live in `tests/integration/` and
are skipped on machines that don't have the cached h5ads.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest


# Import loaders directly; this validates module-level scoping doesn't break.
from confpert.data.frangieh import load_frangieh, FRANGIEH_URL
from confpert.data.datlinger import load_datlinger, DATLINGER_URL
from confpert.data.schmidt import load_schmidt, SCHMIDT_BUNDLE
from confpert.data.mcfaline_figueroa import load_mcfaline_figueroa, MCFALINE_URL
from confpert.data.walker import load_walker, WALKER_SUBSTITUTION_CANDIDATES
from confpert.data.norman import PerturbDataset


def test_frangieh_loader_signature_and_missing_file_message(tmp_path):
    fake = tmp_path / "frangieh_2021_rna.h5ad"
    with pytest.raises(FileNotFoundError) as exc:
        load_frangieh(h5ad_path=fake)
    assert "Zenodo" in str(exc.value) or "zenodo.org" in str(exc.value) or FRANGIEH_URL in str(exc.value)


def test_datlinger_loader_signature_and_missing_file_message(tmp_path):
    fake = tmp_path / "datlinger_2017.h5ad"
    with pytest.raises(FileNotFoundError) as exc:
        load_datlinger(h5ad_path=fake)
    assert "zenodo.org" in str(exc.value) or DATLINGER_URL in str(exc.value)


def test_schmidt_loader_signature_and_missing_bundle_message(tmp_path):
    fake = tmp_path / "schmidt_2022.h5ad"
    # auto_assemble defaults True; bundle_dir defaults to fake.parent, where
    # nothing exists. Should raise FileNotFoundError with all 4 URLs surfaced.
    with pytest.raises(FileNotFoundError) as exc:
        load_schmidt(h5ad_path=fake)
    msg = str(exc.value)
    assert "matrix" in msg.lower()
    assert "barcodes" in msg.lower() or "GEO" in msg


def test_schmidt_loader_auto_assemble_off_yields_url_list(tmp_path):
    fake = tmp_path / "schmidt_2022.h5ad"
    with pytest.raises(FileNotFoundError) as exc:
        load_schmidt(h5ad_path=fake, auto_assemble=False)
    msg = str(exc.value)
    # All 4 SCHMIDT_BUNDLE URLs should be in the error message
    for url in SCHMIDT_BUNDLE.values():
        assert url in msg


def test_mcfaline_loader_signature_and_missing_file_message(tmp_path):
    fake = tmp_path / "mcfaline_2024.h5ad"
    with pytest.raises(FileNotFoundError) as exc:
        load_mcfaline_figueroa(h5ad_path=fake)
    assert MCFALINE_URL in str(exc.value) or "sci-Plex" in str(exc.value)


def test_walker_loader_raises_with_substitution_candidates():
    with pytest.raises(NotImplementedError) as exc:
        load_walker()
    msg = str(exc.value)
    # All 3 substitution candidate URLs must be surfaced for the user.
    for candidate_url in WALKER_SUBSTITUTION_CANDIDATES.values():
        assert candidate_url in msg


def test_loaders_use_perturbdataset_schema():
    """Each implemented loader is documented to return PerturbDataset.
    Verify the dataclass import path stays stable so Modal runner works.
    """
    assert PerturbDataset is not None
    fields = {"X_ctrl", "X_pert", "gene_names", "metadata"}
    # PerturbDataset is a @dataclass with explicit fields
    pd_fields = {f.name for f in PerturbDataset.__dataclass_fields__.values()}
    assert fields.issubset(pd_fields)


def test_walker_substitution_candidates_are_zenodo_scperturb():
    """User-supervised amendment path: candidates must be public Zenodo h5ad."""
    for name, url in WALKER_SUBSTITUTION_CANDIDATES.items():
        assert url.startswith("https://zenodo.org/records/13350497/files/")
        assert url.endswith(".h5ad")


def test_frangieh_zenodo_url_format():
    assert FRANGIEH_URL.startswith("https://zenodo.org/records/")
    assert FRANGIEH_URL.endswith("FrangiehIzar2021_RNA.h5ad")


def test_datlinger_zenodo_url_format():
    assert DATLINGER_URL.startswith("https://zenodo.org/records/")
    assert DATLINGER_URL.endswith("DatlingerBock2017.h5ad")


def test_mcfaline_url_is_geo_raw_tar():
    assert MCFALINE_URL.startswith("https://ftp.ncbi.nlm.nih.gov/")
    assert "GSE225775_RAW.tar" in MCFALINE_URL


def test_schmidt_bundle_has_four_files():
    """matrix + barcodes + features + guidecalls = 4 files per cellranger output."""
    expected_files = {"matrix", "barcodes", "features", "guidecalls"}
    assert set(SCHMIDT_BUNDLE.keys()) == expected_files
    for k, url in SCHMIDT_BUNDLE.items():
        assert "GSE190604" in url
        assert url.endswith(".gz")


# ---------------------------------------------------------------------------
# Real-data integration tests — skipped if cached h5ad not present.
# Run manually with: pytest -m integration -v
# ---------------------------------------------------------------------------

@pytest.mark.integration
def test_frangieh_loads_real_h5ad_if_present():
    cache = Path("/Users/aayanalwani/FLM4S/confpert/data/frangieh_2021_rna.h5ad")
    if not cache.exists():
        pytest.skip(f"Cached Frangieh h5ad missing: {cache}")
    ds = load_frangieh(h5ad_path=cache, n_top_perturbations=5, n_hvg=64)
    assert isinstance(ds, PerturbDataset)
    assert ds.X_ctrl.shape[0] > 0
    assert len(ds.X_pert) > 0
    assert ds.X_ctrl.shape[1] == 64


@pytest.mark.integration
def test_datlinger_loads_real_h5ad_if_present():
    cache = Path("/Users/aayanalwani/FLM4S/confpert/data/datlinger_2017.h5ad")
    if not cache.exists():
        pytest.skip(f"Cached Datlinger h5ad missing: {cache}")
    ds = load_datlinger(h5ad_path=cache, n_top_perturbations=5, n_hvg=64)
    assert isinstance(ds, PerturbDataset)
    assert ds.X_ctrl.shape[0] > 0
    assert len(ds.X_pert) > 0
