"""Four pre-registered noise variants for point-estimate predictors.

Per preregistration.md (commit c7046e4b), Mean predictor and Bilinear ridge produce
per-pert mean expression with no native cell-level variance. To compute
distributional discrepancies (KS, W1, energy, MMD, bimodality, variance ratio)
against them, ConfPert adds noise per a documented variant.

Variants:
  A  no_noise           point mass at the mean (variance ratio = 0 by construction)
  B  isotropic          mu + sigma_global * eps,  eps ~ N(0, I)
  C  per_gene_marginal  mu + sigma_j * eps_j      (calibrated per-gene std)
  D  full_covariance    mu + Sigma^(1/2) * eps    (calibrated covariance)

Each variant is a PARAMETER-FREE calibration step using only train-time observed
per-pert statistics. No test-time leak.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class NoiseSpec:
    name: str  # "no_noise", "isotropic", "per_gene_marginal", "full_covariance"
    description: str


NOISE_VARIANTS = {
    "A_no_noise": NoiseSpec("A_no_noise",
                            "Variant A: point mass at the per-pert mean. VR=0 floor."),
    "B_isotropic": NoiseSpec("B_isotropic",
                             "Variant B: mu + sigma_global * eps, eps ~ N(0, I)."),
    "C_per_gene_marginal": NoiseSpec("C_per_gene_marginal",
                                     "Variant C: mu + sigma_j * eps_j (per-gene std)."),
    "D_full_covariance": NoiseSpec("D_full_covariance",
                                   "Variant D: mu + Sigma^(1/2) * eps (full covariance)."),
}


def apply_noise(mu: np.ndarray, n_cells: int, variant: str,
                X_train_pert: np.ndarray | None = None,
                seed: int = 0) -> np.ndarray:
    """Apply the named noise variant to a per-pert mean vector.

    Args:
      mu: [d_genes] per-pert mean prediction.
      n_cells: number of synthetic cells to generate.
      variant: one of {"A_no_noise", "B_isotropic", "C_per_gene_marginal", "D_full_covariance"}.
      X_train_pert: [n_train_pert_cells, d_genes] observed perturbed-population from
                    train fold, used to fit noise statistics. Required for B/C/D.
      seed: RNG seed.

    Returns:
      [n_cells, d_genes] array of synthetic cells.
    """
    if variant not in NOISE_VARIANTS:
        raise KeyError(f"Unknown variant {variant}. Options: {list(NOISE_VARIANTS)}")

    d = mu.shape[0]
    rng = np.random.RandomState(seed)

    if variant == "A_no_noise":
        return np.tile(mu[None, :], (n_cells, 1)).astype(np.float32)

    if X_train_pert is None or X_train_pert.shape[0] < 2:
        raise ValueError(
            f"Variant {variant} requires X_train_pert with >=2 cells "
            f"to fit noise statistics."
        )

    if variant == "B_isotropic":
        sigma_global = float(X_train_pert.std(ddof=1))
        eps = rng.randn(n_cells, d).astype(np.float32)
        return (mu[None, :] + sigma_global * eps).astype(np.float32)

    if variant == "C_per_gene_marginal":
        sigma_j = X_train_pert.std(axis=0, ddof=1).astype(np.float32)
        eps = rng.randn(n_cells, d).astype(np.float32)
        return (mu[None, :] + sigma_j[None, :] * eps).astype(np.float32)

    if variant == "D_full_covariance":
        # Use Cholesky if PD, else regularize and try again.
        cov = np.cov(X_train_pert, rowvar=False, ddof=1).astype(np.float64)
        # Regularize
        cov = cov + 1e-4 * np.eye(d)
        try:
            L = np.linalg.cholesky(cov)
        except np.linalg.LinAlgError:
            # Eigendecomposition fallback
            w, V = np.linalg.eigh(cov)
            w = np.maximum(w, 0.0)
            L = (V * np.sqrt(w)[None, :])
        eps = rng.randn(n_cells, d).astype(np.float64)
        samples = mu[None, :].astype(np.float64) + eps @ L.T
        return samples.astype(np.float32)

    raise RuntimeError(f"unreachable variant {variant}")
