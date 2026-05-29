"""ConfPert: Calibrated distributional uncertainty for single-cell perturbation predictors.

A model-agnostic conformal wrapper providing finite-sample marginal coverage on
six distributional discrepancies for cell-population perturbation responses:
KS, Wasserstein-1, energy distance, MMD-RBF, bimodality coefficient match, variance ratio.

Three knockout claims (per preregistration.md):
  K1 (benchmark + library): conformal coverage on 8 wrappable predictors x 4 datasets
                            x 3 splits x 6 discrepancies x 3 alpha levels.
  K2 (calibration vs capacity / data): pre-registered Spearman tests on H1 (capacity hurts)
                                       and H1b (data helps).
  K3 (drug-resistance discovery): BH-FDR-corrected discovery on PRISM via DepMap and
                                  Hallmark MSigDB orthogonal validation.

Public API:
  confpert.metrics    -- six distributional discrepancies + aggregator
  confpert.conformal  -- four conformal head levels (per-gene, per-perturbation,
                         per-population, subgroup-conditional) + RCPS risk control
  confpert.predictors -- 8 wrappable predictors + four noise variants for point estimates
  confpert.data       -- 4 dataset loaders (Norman, Replogle K562, Replogle RPE1, Adamson)
                         plus Tahoe-100M subset stub

Phase 0 commit: c7046e4b4c08106acd52c37e7ea8410e126d2ae4 (preregistration locked).
"""
__version__ = "0.2.0rc2"

# Avoid hard import errors when optional deps are missing (e.g. on Modal images
# that don't include torch). Heavy submodules are imported on demand.
__all__ = ["__version__"]
