"""Replogle 2022 RPE1 essential Perturb-seq loader (cross-cell-line generalization).

RPE1 is a non-cancerous retinal pigment epithelial cell line; the K562 datasets
(Norman, Replogle K562 essential, Adamson) are all chronic myelogenous leukemia.
RPE1 lets us check whether the variance gap is K562-specific or generalizes to
non-cancer cell lines.

Reuses the same `PerturbDataset` schema and CTRL_REGEX as Norman/Replogle.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from .norman import CTRL_REGEX, PerturbDataset


def load_replogle_rpe1(
    h5ad_path: str | Path = "/Users/aayanalwani/FLM4S/confpert/data/replogle_rpe1.h5ad",
    n_top_perturbations: int = 50,
    n_hvg: int = 256,
    min_cells_per_pert: int = 30,
    log_normalize: bool = True,
    target_sum: float = 1e4,
    chunked: bool = False,
) -> PerturbDataset:
    """chunked=True opens h5ad in backed mode, identifies the rows we keep
    using only obs (small), then materializes only those rows. Avoids the
    macOS OOM that hit during the 2026-05-03 overnight K1 sweep on 1.2 GB RPE1.
    """
    try:
        import scanpy as sc
    except ImportError as exc:
        raise RuntimeError("scanpy is required") from exc

    h5ad_path = Path(h5ad_path)
    if not h5ad_path.exists():
        raise FileNotFoundError(f"Replogle RPE1 h5ad not found at {h5ad_path}")

    if chunked:
        import anndata as ad
        adata_backed = ad.read_h5ad(h5ad_path, backed="r")
        col_b = None
        for c in ("perturbation", "guide_id", "gene_id", "gene", "guide_identity",
                  "perturbation_name", "condition", "gene_name"):
            if c in adata_backed.obs.columns:
                col_b = c
                break
        if col_b is None:
            raise ValueError(f"No perturbation column in {h5ad_path}: "
                             f"{list(adata_backed.obs.columns)}")
        labels_b = adata_backed.obs[col_b].astype(str)
        is_ctrl_b = labels_b.str.lower().str.contains(CTRL_REGEX, regex=True)
        pert_labels_b = labels_b[~is_ctrl_b]
        counts_b = pert_labels_b.value_counts()
        counts_b = counts_b[counts_b >= min_cells_per_pert]
        top_perts_b = counts_b.head(n_top_perturbations).index.tolist()
        keep_mask = is_ctrl_b | labels_b.isin(top_perts_b)
        keep_idx = np.where(keep_mask.values)[0]
        adata = adata_backed[keep_idx].to_memory()
        adata_backed.file.close()
    else:
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
        raise ValueError(f"Only {n_ctrl} control cells found in column '{col}'.")

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

    return PerturbDataset(
        X_ctrl=X_ctrl, X_pert=X_pert,
        gene_names=adata.var_names.tolist(),
        metadata={
            "source": "Replogle Weissman 2022 *Cell* -- RPE1 essential CRISPRi Perturb-seq",
            "perturbation_column": col,
            "n_top_perturbations_returned": len(X_pert),
            "n_ctrl": int(X_ctrl.shape[0]),
            "n_hvg": X_full.shape[1],
            "log_normalized": log_normalize,
        })
