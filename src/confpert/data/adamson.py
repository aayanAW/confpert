"""Adamson 2016 K562 Perturb-seq loader.

Source: Adamson et al. Cell 2016 "A Multiplexed Single-Cell CRISPR Screening
Platform Enables Systematic Dissection of the Unfolded Protein Response", GSE90546.

Smaller than Norman (~6K cells, 8 perturbations) but the canonical UPR Perturb-seq
benchmark. We pull from scperturb.org (which mirrors the cleaned h5ad).

If the h5ad isn't on the Modal volume, the loader falls back to a download
attempt.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from .norman import CTRL_REGEX, PerturbDataset


ADAMSON_URLS = [
    # scPerturb mirror (preferred)
    "https://zenodo.org/records/13350497/files/AdamsonWeissman2016_GSM2406681_10X005.h5ad",
    "https://zenodo.org/records/13350497/files/AdamsonWeissman2016_GSM2406675_10X001.h5ad",
]


def load_adamson(
    h5ad_path: str | Path = "/Users/aayanalwani/FLM4S/confpert/data/adamson_2016.h5ad",
    n_top_perturbations: int = 30,
    n_hvg: int = 256,
    min_cells_per_pert: int = 30,
    log_normalize: bool = True,
    target_sum: float = 1e4,
) -> PerturbDataset:
    """Load Adamson 2016 K562 UPR Perturb-seq into a PerturbDataset.

    Will not auto-download: expects h5ad at the given path. If missing, prints
    instructions for manual download.
    """
    try:
        import scanpy as sc
    except ImportError as exc:
        raise RuntimeError("scanpy is required") from exc

    h5ad_path = Path(h5ad_path)
    if not h5ad_path.exists():
        urls = "\n  ".join(ADAMSON_URLS)
        raise FileNotFoundError(
            f"Adamson h5ad not found at {h5ad_path}. Manually download from one of:\n"
            f"  {urls}\n"
            f"or pull via Modal volume if available.")

    adata = sc.read_h5ad(h5ad_path)

    col = None
    for c in ("perturbation", "guide_id", "gene_id", "gene", "guide_identity",
              "perturbation_name", "condition", "gene_name"):
        if c in adata.obs.columns:
            col = c
            break
    if col is None:
        raise ValueError(f"No perturbation column found. obs.columns: "
                         f"{list(adata.obs.columns)}")

    labels = adata.obs[col].astype(str)
    is_ctrl = labels.str.lower().str.contains(CTRL_REGEX, regex=True)
    n_ctrl = int(is_ctrl.sum())
    if n_ctrl < 30:
        raise ValueError(f"Only {n_ctrl} control cells found in column '{col}'. "
                         f"Top labels: {labels.value_counts().head(5).to_dict()}")

    pert_labels = labels[~is_ctrl]
    counts = pert_labels.value_counts()
    counts = counts[counts >= min_cells_per_pert]
    top_perts = counts.head(n_top_perturbations).index.tolist()

    keep_mask = is_ctrl | labels.isin(top_perts)
    adata = adata[keep_mask.values].copy()

    sc.pp.filter_cells(adata, min_counts=200)
    sc.pp.filter_genes(adata, min_cells=5)
    if log_normalize:
        sc.pp.normalize_total(adata, target_sum=target_sum)
        sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=n_hvg)
    adata = adata[:, adata.var["highly_variable"]].copy()

    X_full = adata.X.toarray() if hasattr(adata.X, "toarray") else np.asarray(adata.X)
    X_full = X_full.astype(np.float32)

    labels_post = adata.obs[col].astype(str)
    is_ctrl_post = labels_post.str.lower().str.contains(CTRL_REGEX, regex=True)
    X_ctrl = X_full[is_ctrl_post.values]
    X_pert: dict[str, np.ndarray] = {}
    for p in top_perts:
        mask = (labels_post == p).values
        if mask.sum() >= min_cells_per_pert:
            X_pert[p] = X_full[mask]

    metadata = {
        "source": "Adamson et al. Cell 2016 (GSE90546) -- K562 UPR Perturb-seq",
        "perturbation_column": col,
        "n_top_perturbations_requested": n_top_perturbations,
        "n_perturbations_returned": len(X_pert),
        "n_ctrl": int(X_ctrl.shape[0]),
        "n_hvg": X_full.shape[1],
        "log_normalized": log_normalize,
    }
    return PerturbDataset(X_ctrl=X_ctrl, X_pert=X_pert,
                          gene_names=adata.var_names.tolist(), metadata=metadata)
