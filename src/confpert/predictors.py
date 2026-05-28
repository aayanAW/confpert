"""Wrappable predictors for ConfPert (Phase 1 lightweight set).

The eight pre-registered predictors per preregistration.md:

  Lightweight (CPU-tractable, this module):
    1. MeanPredictor                      ~0 params
    2. AhlmannEltzeBilinearRidge          ~10^4 (K=10 PCs of pseudo-bulk Y_train)

  Re-implementation of HetPert lineage (this module):
    3. AdditivePredictor                  (kept as sanity baseline; not in K1 main eight)
    4. NoisyMeanPredictor                 (kept as sanity baseline)

  Sample-producing wrappers (deferred to predictors_modal.py / per-method modules):
    5. scGen retrain
    6. CPA cpa-tools
    7. biolord
    8. sVAE+ / SAMS-VAE
    9. GEARS uncertainty mode
   10. STATE SE-600M + ST checkpoint

Each predictor implements:
  fit(X_ctrl_train, perts_train, X_pert_train) -> self
  predict_samples(X_ctrl_test, perturbation, n_cells) -> [n_cells, d_genes]

The predict_samples() return type is the universal sample-producing interface that
the conformal heads operate on.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from sklearn.linear_model import Ridge

from .noise_models import apply_noise


@dataclass
class PredictorSpec:
    name: str
    description: str
    n_params: int


# ---------------------------------------------------------------------------
# 1. Mean predictor
# ---------------------------------------------------------------------------


class MeanPredictor:
    """Per-pert mean of observed-train cells. Floor for variance ratio.

    Combined with a noise variant {A, B, C, D} at predict time.
    """
    spec = PredictorSpec(
        name="mean",
        description="Per-pert mean, with externally-applied noise variant.",
        n_params=0,
    )

    def __init__(self, noise_variant: str = "A_no_noise", seed: int = 42):
        self.noise_variant = noise_variant
        self.seed = seed
        self.means: dict[str, np.ndarray] = {}
        self.train_pops: dict[str, np.ndarray] = {}

    def fit(self, X_ctrl_train: np.ndarray, perts_train: np.ndarray,
            X_pert_train: np.ndarray) -> "MeanPredictor":
        for p in np.unique(perts_train):
            mask = perts_train == p
            X_p = X_pert_train[mask]
            self.means[p] = X_p.mean(axis=0).astype(np.float32)
            self.train_pops[p] = X_p
        return self

    def predict_samples(self, X_ctrl_test: np.ndarray, perturbation: str,
                        n_cells: int | None = None) -> np.ndarray:
        if perturbation not in self.means:
            raise KeyError(f"Unseen perturbation {perturbation}")
        if n_cells is None:
            n_cells = X_ctrl_test.shape[0]
        return apply_noise(
            self.means[perturbation],
            n_cells,
            variant=self.noise_variant,
            X_train_pert=self.train_pops.get(perturbation),
            seed=self.seed,
        )


# ---------------------------------------------------------------------------
# 2. Ahlmann-Eltze bilinear ridge (eq. 1, 3 from their 2025 paper)
# ---------------------------------------------------------------------------


class AhlmannEltzeBilinearRidge:
    """Low-rank bilinear ridge baseline per Ahlmann-Eltze, Huber, Anders 2025.

    Y_train pseudo-bulk = G W P^T + b, with:
      Y_train: [d_genes, n_perts]  pseudo-bulk per-pert means (genes are rows)
      G:       [d_genes, K]        top-K PCs of Y_train
      P:       [n_perts, K]        rows of G corresponding to perturbed genes
                                   (binary indicator times G for genes targeted)
      W:       [K, K]              learned weights via closed-form ridge with lambda
      b:       [d_genes]           per-gene control mean

    Closed form (eq. 3):
      W = (G^T G + lambda I)^{-1} G^T (Y_train - b) P (P^T P + lambda I)^{-1}

    Predict:
      Y_hat[gene, pert] = G W P_pert^T + b_gene
    where P_pert is the row of P for the new perturbation.

    This baseline is mean-only by construction (one prediction per perturbation, no
    cell-level variance). Combine with a noise variant {A, B, C, D} at predict time.
    """
    spec = PredictorSpec(
        name="ahlmann_bilinear_ridge",
        description="Bilinear ridge per Ahlmann-Eltze 2025 eq. 1,3 (K=10 PCs, lambda=0.1).",
        n_params=int(1e4),
    )

    def __init__(self, K: int = 10, lam: float = 0.1, noise_variant: str = "A_no_noise",
                 seed: int = 42):
        self.K = K
        self.lam = lam
        self.noise_variant = noise_variant
        self.seed = seed
        self.G: np.ndarray | None = None
        self.W: np.ndarray | None = None
        self.b: np.ndarray | None = None
        self.pert_to_idx: dict[str, int] = {}
        self.gene_names: list[str] | None = None
        self.means: dict[str, np.ndarray] = {}
        self.train_pops: dict[str, np.ndarray] = {}

    def fit(self, X_ctrl_train: np.ndarray, perts_train: np.ndarray,
            X_pert_train: np.ndarray, gene_names: list[str] | None = None,
            perturbed_gene_indices: dict[str, list[int]] | None = None) -> "AhlmannEltzeBilinearRidge":
        """Fit the bilinear ridge.

        gene_names: [d_genes] gene labels (HVG selection).
        perturbed_gene_indices: {pert_name: [gene_idx, ...]} indicating which gene
            indices each perturbation targets. If None, fall back to per-pert one-hot
            indicator using a fresh basis derived from per-pert observed means.
        """
        d_genes = X_pert_train.shape[1]
        unique_perts = list(np.unique(perts_train))
        self.pert_to_idx = {p: i for i, p in enumerate(unique_perts)}
        n_perts = len(unique_perts)
        self.gene_names = gene_names

        # Build Y_train pseudo-bulk: [d_genes, n_perts]
        Y_train = np.zeros((d_genes, n_perts), dtype=np.float64)
        for p, i in self.pert_to_idx.items():
            mask = perts_train == p
            X_p = X_pert_train[mask]
            Y_train[:, i] = X_p.mean(axis=0)
            self.train_pops[p] = X_p
            self.means[p] = X_p.mean(axis=0).astype(np.float32)

        # b = per-gene control mean
        b = X_ctrl_train.mean(axis=0).astype(np.float64)
        self.b = b

        # G = top-K PCs of (Y_train - b[:, None])
        Y_centered = Y_train - b[:, None]
        # SVD on the d_genes x n_perts matrix
        U, S, _ = np.linalg.svd(Y_centered, full_matrices=False)
        K_eff = min(self.K, U.shape[1])
        self.G = U[:, :K_eff]

        # P = if perturbed_gene_indices is provided, P[p_idx, k] = sum_g 1{g in pert} G[g, k]
        # else fall back to a Y-based encoding using the column basis from SVD.
        if perturbed_gene_indices is not None and gene_names is not None:
            gname_to_idx = {g: i for i, g in enumerate(gene_names)}
            P = np.zeros((n_perts, K_eff), dtype=np.float64)
            for p, p_idx in self.pert_to_idx.items():
                indices = []
                for g in perturbed_gene_indices.get(p, []):
                    if g in gname_to_idx:
                        indices.append(gname_to_idx[g])
                if not indices:
                    # Fallback: use the corresponding column in the V^T basis from SVD
                    # (the per-perturbation row in the right-singular subspace)
                    _, _, Vt = np.linalg.svd(Y_centered, full_matrices=False)
                    P[p_idx, :] = Vt[:K_eff, p_idx]
                else:
                    P[p_idx, :] = self.G[indices, :].sum(axis=0)
        else:
            # Use right-singular vectors as P
            _, _, Vt = np.linalg.svd(Y_centered, full_matrices=False)
            P = Vt[:K_eff, :].T  # [n_perts, K_eff]

        # Closed-form W per eq. 3:
        # W = (G^T G + lambda I)^{-1} G^T (Y_train - b) P (P^T P + lambda I)^{-1}
        GtG_inv = np.linalg.inv(self.G.T @ self.G + self.lam * np.eye(K_eff))
        PtP_inv = np.linalg.inv(P.T @ P + self.lam * np.eye(K_eff))
        self.W = GtG_inv @ self.G.T @ Y_centered @ P @ PtP_inv

        # Cache P for predict
        self._P = P

        return self

    def predict_mean(self, perturbation: str) -> np.ndarray:
        """Y_hat[:, pert] = G W P_pert^T + b. Returns [d_genes]."""
        if self.W is None:
            raise RuntimeError("fit() first")
        idx = self.pert_to_idx.get(perturbation)
        if idx is None:
            # Out-of-distribution perturbation: fall back to control mean.
            return self.b.astype(np.float32)
        P_pert = self._P[idx, :]  # [K_eff]
        y_hat = self.G @ self.W @ P_pert + self.b
        return y_hat.astype(np.float32)

    def predict_samples(self, X_ctrl_test: np.ndarray, perturbation: str,
                        n_cells: int | None = None) -> np.ndarray:
        if n_cells is None:
            n_cells = X_ctrl_test.shape[0]
        mu = self.predict_mean(perturbation)
        return apply_noise(
            mu,
            n_cells,
            variant=self.noise_variant,
            X_train_pert=self.train_pops.get(perturbation),
            seed=self.seed,
        )


# ---------------------------------------------------------------------------
# Sanity baselines (kept from HetPert lineage)
# ---------------------------------------------------------------------------


class AdditivePredictor:
    """X_pert ~= X_ctrl + delta(perturbation), where delta = mean(pert) - mean(ctrl)."""
    spec = PredictorSpec(
        name="additive",
        description="ctrl_cell + (mean_pert - mean_ctrl). Sanity baseline.",
        n_params=1,
    )

    def __init__(self):
        self.deltas: dict[str, np.ndarray] = {}
        self.ctrl_mean: np.ndarray | None = None

    def fit(self, X_ctrl_train: np.ndarray, perts_train: np.ndarray,
            X_pert_train: np.ndarray) -> "AdditivePredictor":
        self.ctrl_mean = X_ctrl_train.mean(axis=0).astype(np.float32)
        for p in np.unique(perts_train):
            mask = perts_train == p
            mu_p = X_pert_train[mask].mean(axis=0)
            self.deltas[p] = (mu_p - self.ctrl_mean).astype(np.float32)
        return self

    def predict_samples(self, X_ctrl_test: np.ndarray, perturbation: str,
                        n_cells: int | None = None) -> np.ndarray:
        if perturbation not in self.deltas:
            raise KeyError(perturbation)
        delta = self.deltas[perturbation]
        if n_cells is None or n_cells == X_ctrl_test.shape[0]:
            return X_ctrl_test + delta[None, :]
        # Sample with replacement if n_cells differs
        rng = np.random.RandomState(0)
        idx = rng.randint(0, X_ctrl_test.shape[0], size=n_cells)
        return X_ctrl_test[idx] + delta[None, :]


class NoisyMeanPredictor:
    """Mean + per-gene observed-std Gaussian noise. Equivalent to MeanPredictor + Variant C."""
    spec = PredictorSpec(
        name="noisy_mean",
        description="Per-pert mean + per-gene calibrated Gaussian noise.",
        n_params=1,
    )

    def __init__(self, seed: int = 42):
        self.seed = seed
        self.means: dict[str, np.ndarray] = {}
        self.stds: dict[str, np.ndarray] = {}

    def fit(self, X_ctrl_train: np.ndarray, perts_train: np.ndarray,
            X_pert_train: np.ndarray) -> "NoisyMeanPredictor":
        for p in np.unique(perts_train):
            mask = perts_train == p
            X_p = X_pert_train[mask]
            self.means[p] = X_p.mean(axis=0).astype(np.float32)
            self.stds[p] = X_p.std(axis=0, ddof=1).astype(np.float32)
        return self

    def predict_samples(self, X_ctrl_test: np.ndarray, perturbation: str,
                        n_cells: int | None = None) -> np.ndarray:
        if perturbation not in self.means:
            raise KeyError(perturbation)
        if n_cells is None:
            n_cells = X_ctrl_test.shape[0]
        rng = np.random.RandomState(self.seed)
        mu = self.means[perturbation]
        sd = self.stds[perturbation]
        eps = rng.randn(n_cells, mu.shape[0]).astype(np.float32)
        return (mu[None, :] + sd[None, :] * eps).astype(np.float32)


PREDICTORS = {
    "mean": MeanPredictor,
    "ahlmann_bilinear_ridge": AhlmannEltzeBilinearRidge,
    "additive": AdditivePredictor,
    "noisy_mean": NoisyMeanPredictor,
}
