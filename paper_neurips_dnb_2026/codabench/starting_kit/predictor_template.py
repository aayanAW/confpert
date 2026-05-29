"""Reference predictor template for ConfPert Codabench submissions.

Implement `predict_samples` and package the result as a Docker image with
`Dockerfile` (this directory). The scorer will load this module and call
`predict_samples` per (dataset, perturbation) cell.

The contract:

    def predict_samples(
        X_ctrl_test: np.ndarray,   # (n_cells, n_genes), control cells
                                    #   from the same dataset and split
        perturbation: str,         # gene name (genetic) or compound id (chemical)
        n_cells: int,              # number of cells to sample
    ) -> np.ndarray:               # (n_cells, n_genes) — predicted perturbed pop

Any model is acceptable: pretrained, train-from-scratch, deterministic baseline,
ensemble. The image must be self-contained — no internet access at scoring time.

Example implementations:

  - Trivial mean baseline (this file's default): predict mu = mean(X_ctrl_test),
    sample with low isotropic noise.
  - scGen / CPA / scGPT: load pretrained checkpoint from the image, run forward
    with the perturbation as conditioning input.
"""
from __future__ import annotations

import numpy as np


def predict_samples(
    X_ctrl_test: np.ndarray,
    perturbation: str,
    n_cells: int,
) -> np.ndarray:
    """Trivial mean baseline: replicate the control mean with small Gaussian
    noise. Override this to plug in your own predictor.

    Parameters
    ----------
    X_ctrl_test : np.ndarray
        (n_ctrl_cells, n_genes) control cells for the held-out fold.
    perturbation : str
        The held-out perturbation label (gene name or compound id).
    n_cells : int
        How many cells to sample.

    Returns
    -------
    np.ndarray
        (n_cells, n_genes) sampled predicted perturbed population.
    """
    if X_ctrl_test.ndim != 2:
        raise ValueError(f"X_ctrl_test must be 2D, got shape {X_ctrl_test.shape}")
    if n_cells <= 0:
        raise ValueError(f"n_cells must be > 0, got {n_cells}")

    mu = X_ctrl_test.mean(axis=0)
    sigma = X_ctrl_test.std(axis=0) + 1e-6
    rng = np.random.default_rng(42)
    return mu[None, :] + rng.standard_normal(size=(n_cells, mu.size)) * sigma[None, :]


if __name__ == "__main__":
    # Sanity-check round-trip
    X = np.random.RandomState(0).randn(50, 16)
    out = predict_samples(X, "MYC", n_cells=32)
    assert out.shape == (32, 16)
    print(f"OK -- predict_samples returns shape {out.shape}")
