"""Norman 2019 K562 Perturb-seq loader (single + dual CRISPR perturbations).

Source: Norman et al. Cell 2019, GSE133344. We expect the cached h5ad at
~/FLM4S/hetpert/data/norman_2019.h5ad (or anywhere passed in via path).

Usage:
    ds = load_norman(n_top_perturbations=50, n_hvg=512)
    X_ctrl = ds.X_ctrl                                # [n_ctrl, d_genes]
    X_pert_a = ds.X_pert["MYC"]                       # [n_pert_cells, d_genes]
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np


# Try several common control-marker patterns (Norman uses 'control', Replogle uses
# 'control', Adamson uses 'NTC' / 'non-targeting').
CTRL_REGEX = "ctrl|control|safe|ntc|non-target|non_target|nontarget|no-target|no_target|notarget|scramble"


@dataclass
class PerturbDataset:
    X_ctrl: np.ndarray                        # [n_ctrl, d_genes]
    X_pert: dict[str, np.ndarray]             # pert -> [n_pert_cells, d_genes]
    gene_names: list[str]
    metadata: dict

    def perturbations(self) -> list[str]:
        return list(self.X_pert.keys())

    def n_cells_per_pert(self) -> dict[str, int]:
        return {p: X.shape[0] for p, X in self.X_pert.items()}


def load_norman(
    h5ad_path: str | Path = "/Users/aayanalwani/FLM4S/confpert/data/norman_2019.h5ad",
    n_top_perturbations: int = 50,
    n_hvg: int = 512,
    min_cells_per_pert: int = 30,
    log_normalize: bool = True,
    target_sum: float = 1e4,
) -> PerturbDataset:
    """Load Norman 2019 K562 Perturb-seq into a `PerturbDataset`.

    Selects top `n_top_perturbations` by cell count (excluding control), filters
    perturbations with fewer than `min_cells_per_pert` cells, log-normalizes,
    and selects top `n_hvg` highly variable genes.
    """
    try:
        import scanpy as sc
    except ImportError as exc:
        raise RuntimeError("scanpy is required. pip install scanpy") from exc

    h5ad_path = Path(h5ad_path)
    if not h5ad_path.exists():
        raise FileNotFoundError(f"Norman h5ad not found at {h5ad_path}. "
                                f"Pull from Modal volume: "
                                f"`modal volume get causeflow-artifacts "
                                f"/data/norman_2019.h5ad {h5ad_path}`")

    adata = sc.read_h5ad(h5ad_path)

    # Find perturbation column
    col = None
    for c in ("perturbation", "guide_id", "gene_id", "gene", "guide_identity",
              "perturbation_name", "condition"):
        if c in adata.obs.columns:
            col = c
            break
    if col is None:
        raise ValueError(f"No perturbation column found. obs.columns: "
                         f"{list(adata.obs.columns)}")

    labels = adata.obs[col].astype(str)
    is_ctrl = labels.str.lower().str.contains(CTRL_REGEX, regex=True)

    n_ctrl = int(is_ctrl.sum())
    if n_ctrl < 50:
        raise ValueError(f"Only {n_ctrl} control cells found in column '{col}'. "
                         f"Top labels: {labels.value_counts().head(5).to_dict()}")

    # Top-N perturbations by cell count, excluding ctrl
    pert_labels = labels[~is_ctrl]
    counts = pert_labels.value_counts()
    counts = counts[counts >= min_cells_per_pert]
    top_perts = counts.head(n_top_perturbations).index.tolist()

    # Subset
    keep_mask = is_ctrl | labels.isin(top_perts)
    adata = adata[keep_mask.values].copy()

    # Preprocess
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
        "source": "Norman et al. Cell 2019 (GSE133344) -- K562 Perturb-seq",
        "perturbation_column": col,
        "n_top_perturbations_requested": n_top_perturbations,
        "n_perturbations_returned": len(X_pert),
        "n_ctrl": int(X_ctrl.shape[0]),
        "n_hvg": X_full.shape[1],
        "log_normalized": log_normalize,
        "target_sum": target_sum if log_normalize else None,
    }

    return PerturbDataset(
        X_ctrl=X_ctrl,
        X_pert=X_pert,
        gene_names=adata.var_names.tolist(),
        metadata=metadata,
    )


class NormanDataset:
    """Backward-compat wrapper. Prefer load_norman() above."""

    def __init__(self, *args, **kwargs):
        self._ds = load_norman(*args, **kwargs)

    def __getattr__(self, item):
        return getattr(self._ds, item)
