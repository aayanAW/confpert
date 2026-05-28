"""Tahoe-100M loader (Vevo / Arc Virtual Cell Atlas, 2025).

Source: tahoebio/Tahoe-100M on Hugging Face (CC0-1.0). 100,648,790 cells x
50 cancer cell lines x 1,100 small-molecule perturbations x 3 doses, 24h
endpoint, single-cell-village Mosaic platform. See lit_notes/tahoe100m_2025.md.

The full corpus is 429 GB on HF; we never load it directly. Instead we rely on
``scripts/download_tahoe_subset.py`` (run once on Modal) to assemble a
PRISM-overlapping h5ad subset, typically:

    - drugs:       ~50 small molecules whose Tahoe `drug` name matches a PRISM
                   24Q2 ``Drug.Name`` (case-insensitive).
    - cell lines:  top ~10 by cell count among the matched drugs (DepMap IDs
                   intersected with PRISM screen).
    - dose:        single highest dose per drug for clean signal (configurable).
    - cells:       ~500 K total (~10 % of full corpus per compute_estimate.md).

The subset h5ad is laid out exactly like Norman / Replogle so the same
``PerturbDataset`` schema applies. Two view modes:

    - ``cell_line=None`` (default): chemical perturbation per drug, cells from
      ALL selected cell lines pooled. ``X_ctrl`` = DMSO cells across the same
      lines. This is the K3 view: train one chemical predictor that sees the
      cross-cell-line distribution, transfer to PRISM.
    - ``cell_line=<CVCL_id or DepMap_id>``: filter to one line. This is the K1
      view: drop Tahoe into the K1 sweep as a fourth single-context dataset.

Pre-reg compliance:
    The K1 dataset list locks "Tahoe-100M cross-cell-line subset"
    (preregistration.md L30). We register Tahoe under that label in the K1
    sweep. The cross-cell-line evaluation is achieved at K1 split-time by
    holding out specific cell lines for the calibration / test arms (HetPert
    pattern).
"""
from __future__ import annotations

from pathlib import Path
from typing import Optional

import numpy as np

from .norman import PerturbDataset


# Tahoe uses ``DMSO_TF`` as the in-plate vehicle control marker, plus we
# allow plain ``DMSO`` for downstream re-canonicalised h5ads. ``str.contains``
# does substring matching, so we don't anchor — ``ctrl`` matches ``ctrl_*``
# perturbations e.g. Norman ``ctrl`` or Replogle ``control_001``.
TAHOE_CTRL_REGEX = "dmso|vehicle|control|ctrl"

DEFAULT_TAHOE_PATH = "/Users/aayanalwani/FLM4S/confpert/data/tahoe_subset.h5ad"


def load_tahoe(
    h5ad_path: str | Path = DEFAULT_TAHOE_PATH,
    n_top_perturbations: int = 50,
    n_hvg: int = 512,
    min_cells_per_pert: int = 30,
    log_normalize: bool = True,
    target_sum: float = 1e4,
    cell_line: Optional[str] = None,
    drug_col: str = "drug",
    cell_line_col: str = "cell_line_id",
) -> PerturbDataset:
    """Load a pre-built Tahoe-100M subset h5ad into a ``PerturbDataset``.

    Args:
        h5ad_path: path to subset h5ad produced by
            ``scripts/download_tahoe_subset.py``. Must already be a small,
            in-memory-tractable file (~hundreds of MB to few GB).
        n_top_perturbations: keep top-N drugs by cell count.
        n_hvg: top-N highly variable genes.
        min_cells_per_pert: drop drugs with fewer than this many cells.
        log_normalize: ``sc.pp.normalize_total`` + ``log1p``.
        target_sum: target_sum for normalize_total.
        cell_line: if given, filter to this single cell line value (matches
            against ``cell_line_col``). Used for the K1 single-context view.
        drug_col: name of the perturbation column in ``adata.obs``. Tahoe
            convention is ``"drug"``.
        cell_line_col: name of the cell-line column. Tahoe convention is
            ``"cell_line_id"`` (Cellosaurus IDs e.g. ``CVCL_0023``).
    """
    try:
        import scanpy as sc
    except ImportError as exc:
        raise RuntimeError("scanpy is required. pip install scanpy") from exc

    h5ad_path = Path(h5ad_path)
    if not h5ad_path.exists():
        raise FileNotFoundError(
            f"Tahoe subset h5ad not found at {h5ad_path}. Build it on Modal:\n"
            f"  modal run --detach scripts/modal_launch.py::tahoe_subset_build\n"
            f"then `modal volume get causeflow-artifacts "
            f"/data/tahoe_subset.h5ad {h5ad_path}`."
        )

    adata = sc.read_h5ad(h5ad_path)

    if drug_col not in adata.obs.columns:
        raise ValueError(
            f"Drug column '{drug_col}' not in obs.columns: "
            f"{list(adata.obs.columns)}"
        )

    if cell_line is not None:
        if cell_line_col not in adata.obs.columns:
            raise ValueError(
                f"cell_line filter requested but '{cell_line_col}' not in "
                f"obs.columns: {list(adata.obs.columns)}"
            )
        line_mask = (adata.obs[cell_line_col].astype(str) == str(cell_line)).values
        if line_mask.sum() < 100:
            raise ValueError(
                f"cell_line='{cell_line}' yields only {int(line_mask.sum())} cells; "
                f"available lines: "
                f"{adata.obs[cell_line_col].astype(str).value_counts().head(10).to_dict()}"
            )
        adata = adata[line_mask].copy()

    labels = adata.obs[drug_col].astype(str)
    is_ctrl = labels.str.lower().str.contains(TAHOE_CTRL_REGEX, regex=True)

    n_ctrl = int(is_ctrl.sum())
    if n_ctrl < 50:
        raise ValueError(
            f"Only {n_ctrl} DMSO/control cells in '{drug_col}'. Top labels: "
            f"{labels.value_counts().head(5).to_dict()}"
        )

    pert_labels = labels[~is_ctrl]
    counts = pert_labels.value_counts()
    counts = counts[counts >= min_cells_per_pert]
    top_perts = counts.head(n_top_perturbations).index.tolist()
    if not top_perts:
        raise ValueError(
            f"No drugs meet min_cells_per_pert={min_cells_per_pert}. Top: "
            f"{pert_labels.value_counts().head(5).to_dict()}"
        )

    keep_mask = is_ctrl | labels.isin(top_perts)
    adata = adata[keep_mask.values].copy()

    sc.pp.filter_cells(adata, min_counts=200)
    sc.pp.filter_genes(adata, min_cells=5)
    if log_normalize:
        sc.pp.normalize_total(adata, target_sum=target_sum)
        sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=n_hvg)
    adata = adata[:, adata.var["highly_variable"]].copy()

    X_full = (adata.X.toarray() if hasattr(adata.X, "toarray")
              else np.asarray(adata.X)).astype(np.float32)

    labels_post = adata.obs[drug_col].astype(str)
    is_ctrl_post = labels_post.str.lower().str.contains(TAHOE_CTRL_REGEX, regex=True)
    X_ctrl = X_full[is_ctrl_post.values]

    X_pert: dict[str, np.ndarray] = {}
    for p in top_perts:
        mask = (labels_post == p).values
        if mask.sum() >= min_cells_per_pert:
            X_pert[p] = X_full[mask]

    metadata: dict = {
        "source": "Tahoe-100M (Vevo / Arc Virtual Cell Atlas 2025) -- "
                  "drug-perturbation single-cell atlas, 50 cancer lines, "
                  "1100 small molecules, subset",
        "perturbation_column": drug_col,
        "cell_line_column": cell_line_col,
        "cell_line_filter": cell_line,
        "n_top_perturbations_requested": n_top_perturbations,
        "n_perturbations_returned": len(X_pert),
        "n_ctrl": int(X_ctrl.shape[0]),
        "n_hvg": X_full.shape[1],
        "log_normalized": log_normalize,
        "target_sum": target_sum if log_normalize else None,
    }
    if cell_line_col in adata.obs.columns:
        metadata["cell_lines_present"] = sorted(
            adata.obs[cell_line_col].astype(str).unique().tolist()
        )

    return PerturbDataset(
        X_ctrl=X_ctrl,
        X_pert=X_pert,
        gene_names=adata.var_names.tolist(),
        metadata=metadata,
    )
