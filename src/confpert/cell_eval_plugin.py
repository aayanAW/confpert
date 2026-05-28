"""ConfPert plug-in for ArcInstitute/cell-eval.

Per `baselines/cell_eval_audit.md` (commit `2db6b9c6`, audited 2026-05-03), 5 of
the 6 ConfPert distributional discrepancies are absent from cell-eval HEAD:
KS, Wasserstein-1, MMD-RBF, bimodality coefficient match, variance ratio. Energy
distance is present but in a "Pearson correlation across perturbations" framing
(`pearson_edistance`); ConfPert's `score_energy` (per-perturbation absolute
values) is registered separately as `confpert_energy`.

Usage:

    # In any pipeline using cell-eval, just import this module to register the
    # six confpert_* metrics with cell_eval.metrics_registry.
    import confpert.cell_eval_plugin  # noqa: F401

    from cell_eval import MetricsEvaluator
    evaluator = MetricsEvaluator(...)
    results = evaluator.compute(["confpert_ks", "confpert_w1", ...])

The `confpert_` prefix prevents naming collisions if cell-eval upstream adds a
same-named metric.

Per-perturbation scores are returned as `dict[str, float]` keyed by
perturbation name (matching cell-eval's pearson_delta convention).
"""
from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from scipy.sparse import issparse

from .metrics import (
    score_bimodality_mismatch,
    score_energy,
    score_ks,
    score_mmd_rbf,
    score_variance_ratio_dev,
    score_w1,
)

if TYPE_CHECKING:
    from cell_eval._types import PerturbationAnndataPair


def _per_pert_X(adata, mask) -> np.ndarray:
    """Materialize the cell-by-gene block for one perturbation as a dense float32."""
    X = adata.X[mask]
    if issparse(X):
        X = X.toarray()
    return np.asarray(X, dtype=np.float32)


def _per_perturbation(data: "PerturbationAnndataPair", score_fn,
                      skip_control: bool = True) -> dict[str, float]:
    """Apply a ConfPert score function per perturbation.

    score_fn must have signature score_fn(X_pred, X_obs) -> float.
    Returns {pert_name: score}; perturbations with <2 cells in either side or
    the control are skipped (or returned as NaN).
    """
    out: dict[str, float] = {}
    for p in data.perts:
        if skip_control and p == data.control_pert:
            continue
        m_real = data.pert_mask_real.get(p)
        m_pred = data.pert_mask_pred.get(p)
        if m_real is None or m_pred is None:
            out[p] = float("nan")
            continue
        X_real = _per_pert_X(data.real, m_real)
        X_pred = _per_pert_X(data.pred, m_pred)
        if X_real.shape[0] < 2 or X_pred.shape[0] < 2:
            out[p] = float("nan")
            continue
        try:
            out[p] = float(score_fn(X_pred, X_real))
        except Exception as e:
            out[p] = float("nan")
    return out


def confpert_ks(data: "PerturbationAnndataPair", embed_key: str | None = None
                ) -> dict[str, float]:
    """Per-perturbation, per-gene Kolmogorov-Smirnov, mean over genes.

    Lower is better (0 = perfect match). Computed via scipy.stats.ks_2samp.
    """
    return _per_perturbation(data, score_ks)


def confpert_w1(data: "PerturbationAnndataPair", embed_key: str | None = None
                ) -> dict[str, float]:
    """Per-perturbation, per-gene Wasserstein-1, mean over genes.

    Lower is better (0 = perfect match). Computed via scipy.stats.wasserstein_distance.
    """
    return _per_perturbation(data, score_w1)


def confpert_energy(data: "PerturbationAnndataPair", embed_key: str | None = None
                    ) -> dict[str, float]:
    """Per-perturbation, per-gene energy distance (Szekely-Rizzo), mean over genes.

    Differs from cell-eval's `pearson_edistance` (which reports a Pearson
    correlation across perturbations); ConfPert's variant returns the raw
    per-perturbation energy distance.
    """
    return _per_perturbation(data, score_energy)


def confpert_mmd_rbf(data: "PerturbationAnndataPair", embed_key: str | None = None
                      ) -> dict[str, float]:
    """Per-perturbation MMD-RBF (Gretton 2012, median-bandwidth heuristic).

    Lower is better (0 = perfect distributional match). Subsamples to 500
    cells per population for numerical stability.
    """
    return _per_perturbation(data, score_mmd_rbf)


def confpert_bimodality_match(data: "PerturbationAnndataPair",
                               embed_key: str | None = None) -> dict[str, float]:
    """1 minus per-perturbation bimodality classification accuracy (SAS coef
    threshold b > 5/9). Lower is better (0 = perfect bimodality recovery).
    """
    return _per_perturbation(data, score_bimodality_mismatch)


def confpert_variance_ratio_dev(data: "PerturbationAnndataPair",
                                 embed_key: str | None = None) -> dict[str, float]:
    """Per-perturbation, per-gene |variance ratio - 1|, mean over genes.

    Lower is better (0 = perfect variance recovery). The classic mean-collapse
    diagnostic from the HetPert lineage.
    """
    return _per_perturbation(data, score_variance_ratio_dev)


# Auto-register on import (with the standard ConfPert prefix to avoid naming
# collisions if cell-eval upstream adds a same-named metric).
def _register() -> None:
    try:
        from cell_eval._types import MetricBestValue, MetricType
        from cell_eval.metrics import metrics_registry
    except ImportError:
        # cell-eval not installed; the user can still call the metric functions
        # directly. Skip registration.
        return

    metrics_registry.register(
        name="confpert_ks",
        metric_type=MetricType.ANNDATA_PAIR,
        description=(
            "Per-perturbation per-gene Kolmogorov-Smirnov, mean over genes "
            "(scipy.stats.ks_2samp). Lower is better."
        ),
        best_value=MetricBestValue.ZERO,
        func=confpert_ks,
    )
    metrics_registry.register(
        name="confpert_w1",
        metric_type=MetricType.ANNDATA_PAIR,
        description=(
            "Per-perturbation per-gene Wasserstein-1, mean over genes "
            "(scipy.stats.wasserstein_distance). Lower is better."
        ),
        best_value=MetricBestValue.ZERO,
        func=confpert_w1,
    )
    metrics_registry.register(
        name="confpert_energy",
        metric_type=MetricType.ANNDATA_PAIR,
        description=(
            "Per-perturbation per-gene Szekely-Rizzo energy distance "
            "(scipy.stats.energy_distance). Differs from pearson_edistance, "
            "which is a Pearson correlation across perturbations. Lower is better."
        ),
        best_value=MetricBestValue.ZERO,
        func=confpert_energy,
    )
    metrics_registry.register(
        name="confpert_mmd_rbf",
        metric_type=MetricType.ANNDATA_PAIR,
        description=(
            "Per-perturbation MMD-RBF (Gretton 2012, median-bandwidth heuristic). "
            "Subsamples to 500 cells per population. Lower is better."
        ),
        best_value=MetricBestValue.ZERO,
        func=confpert_mmd_rbf,
    )
    metrics_registry.register(
        name="confpert_bimodality_match",
        metric_type=MetricType.ANNDATA_PAIR,
        description=(
            "1 - bimodality coefficient classification accuracy at b > 5/9 "
            "(SAS bimodality coefficient). Lower is better."
        ),
        best_value=MetricBestValue.ZERO,
        func=confpert_bimodality_match,
    )
    metrics_registry.register(
        name="confpert_variance_ratio_dev",
        metric_type=MetricType.ANNDATA_PAIR,
        description=(
            "Per-perturbation per-gene |var(X_pred)/var(X_obs) - 1|, mean over "
            "genes. The HetPert mean-collapse diagnostic. Lower is better."
        ),
        best_value=MetricBestValue.ZERO,
        func=confpert_variance_ratio_dev,
    )


_register()
