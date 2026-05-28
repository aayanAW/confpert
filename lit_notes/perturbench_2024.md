# Wu et al. 2024 - PerturBench

**Venue:** NeurIPS Datasets and Benchmarks 2024. arXiv:2408.10609.

## Claim

PerturBench is a standardized benchmark for single-cell perturbation response prediction. It packages a modular model development and evaluation platform, five curated perturbational datasets, three task types (covariate transfer, combinatorial prediction, data scaling and imbalance), and a metric suite that combines population-level fit metrics with novel rank metrics. The authors evaluate published models (CPA, SAMS-VAE, BioLord, GEARS, scGPT-embedding) against simple baselines (Linear, Latent Additive, Decoder-Only) and report that no single architecture dominates while simpler architectures stay competitive and scale well. They explicitly flag mode collapse in popular published models and argue that rank metrics complement RMSE for diagnosing failure modes.

## Method

Benchmark organization centers on three task types:

- Covariate transfer: train on perturbation effects measured in some cell lines, predict those effects in held-out cell lines where the perturbations were never observed.
- Combo prediction: train on individual perturbation effects plus a fraction of combinations, predict effects of unseen combinations.
- Data scaling and imbalance: vary training set size and per-covariate sampling balance to test robustness.

Standard splits per task:

- Srivatsan20 and Frangieh21 covariate transfer: hold out 30 percent of perturbations such that held-out perturbations are still observed in at least one other cell type.
- Jiang24: hold out 70 percent of covariates across 9 cell states.
- McFalineFigueroa23: nested subsets (small subset of medium subset of full) for scaling curves.
- Norman19 combo: train on all singletons plus 30 percent of dual perturbations, hold out 70 percent of duals for evaluation.
- Srivatsan20 imbalance: downsample perturbations per cell type to graded imbalance levels measured by normalized entropy.

Reference baselines used as the comparison floor are Linear (one-hot encoded perturbation regression), Latent Additive (separate encoders for expression and perturbation that are summed), and Decoder Only (perturbation and covariate inputs without baseline expression). Published models are reimplemented to isolate the core modeling component, marked with asterisk.

## Datasets / experiments

Five public datasets:

- Norman19: 155 single plus 131 dual genetic perturbations, 1 cell state, used for combo prediction.
- Srivatsan20: 188 chemical perturbations, 3 cell states, used for covariate transfer and imbalance.
- Frangieh21: 248 genetic perturbations, 3 cell states.
- McFalineFigueroa23: 525 genetic perturbations, 15 cell states, used for scaling.
- Jiang24: 219 genetic perturbations, 30 cell states, large-scale covariate transfer.

Models benchmarked: CPA, SAMS-VAE, BioLord, GEARS, scGPT (embedding-based latent additive variant), plus Linear, Latent Additive (LA), and Decoder Only baselines.

Concrete leaderboard numbers:

- Srivatsan20 covariate transfer: LA (scGPT) wins with cosine LogFC 0.50 plus or minus 0.004 and RMSE 0.017 plus or minus 0.0001; cosine rank 0.13 plus or minus 0.007. BioLord variance is high at 0.18 plus or minus 0.1. Decoder-only with no perturbation info still scores cosine 0.30.
- Frangieh21: LA best on cosine (0.17 plus or minus 0.006) and RMSE (0.024 plus or minus 0.0004). CPA and SAMS-VAE collapse on rank metrics (approximately 0.42 to 0.48 where 0.5 is random). Decoder Only wins rank at 0.15 plus or minus 0.0004, which the authors flag as a mode-collapse warning.
- Norman19 combo: LA best with cosine 0.79 plus or minus 0.01 and RMSE rank 0.014 plus or minus 0.001. Linear stays competitive at cosine 0.60 plus or minus 0.02. CPA reimplementation drops to 0.52 plus or minus 0.06.
- McFalineFigueroa23 scaling: most models improve with more data, CPA fails to improve on cosine rank, Linear is surprisingly best on rank.
- Srivatsan20 imbalance: Linear acceptable when balanced and fails under imbalance, CPA robust to imbalance but shows mode collapse, scGPT embeddings buffer LA against imbalance.

## Metrics

Population-level fit:

- RMSE on predicted vs observed mean expression vectors per perturbation.
- Cosine similarity on log fold changes (cosine LogFC) between predicted and observed perturbation effects.

Rank metrics (claimed as a novel contribution):

- For a given observed perturbation, the prediction for that perturbation should be more similar to the true effect than predictions for other perturbations are.
- Rank metric counts the fraction of other perturbations whose predicted effect is closer to the true effect than the prediction for the target perturbation. Range 0 (perfect) to 1 (worst), with random approximately 0.5.
- Two rank variants: cosine LogFC rank and RMSE rank.

Hyperparameter optimization loss combines them as L_HPO = RMSE + 0.1 times rank_RMSE.

No conformal prediction, calibration, coverage, prediction interval, or distributional uncertainty metric is included. The authors explicitly note this gap in the limitations and suggest MMD or Wasserstein as future work to capture response heterogeneity.

## What we steal for ConfPert

PerturBench is the standard NeurIPS-grade leaderboard format that ConfPert plugs into. Concrete adoption:

- Three task types map cleanly onto ConfPert evaluation axes: covariate transfer (held-out cell line, distribution shift in covariates), combo prediction (held-out perturbation combinations), and data scaling and imbalance (calibration set size sensitivity).
- Adopt their split protocols verbatim: 30 percent held-out perturbations for Srivatsan20 and Frangieh21, 70 percent held-out covariates for Jiang24, nested McFalineFigueroa23 subsets, and the Norman19 30 percent train and 70 percent test combo split.
- Adopt their datasets as the K1 evaluation suite (all five) and adopt their model zoo as the base predictors that ConfPert wraps (CPA, SAMS-VAE, BioLord, GEARS, scGPT-LA, Linear, LA, Decoder Only).
- Adopt RMSE, cosine LogFC, and the rank metrics as the point-prediction baselines.
- Add coverage and prediction-set size as new ConfPert metrics columns: marginal coverage, conditional coverage stratified by perturbation and cell line, average interval width, and width-normalized coverage. Calibrated coverage is the headline new metric.
- Cite PerturBench as the standard benchmark we extend, not replace. ConfPert K1 result is a PerturBench-format leaderboard with extra calibration columns.

## What we wrap

Not a model. PerturBench is the evaluation harness. ConfPert wraps the benchmarked predictors (CPA, GEARS, BioLord, SAMS-VAE, scGPT-LA, Linear, LA, Decoder Only) with split conformal and weighted conformal layers and reports through PerturBench's metric pipeline. The repo at github.com/altoslabs/perturbench is the codebase we fork or extend.

## Failure modes / caveats

- All metrics are population-level (mean expression per perturbation), so they cannot detect single-cell heterogeneity or distributional miscalibration. The authors acknowledge this directly: future work should use MMD or Wasserstein.
- Published models are reimplemented to isolate the core modeling component, so reported numbers may differ from original-author implementations. This is flagged explicitly in their limitations.
- Matching is limited by knowledge of relevant covariates. Cell cycle is named as a confounder usually unavailable as a covariate.
- No uncertainty, calibration, coverage, or prediction-interval metric is present. This is the gap ConfPert fills.
- Mode collapse is detectable only via the rank metric divergence from RMSE (Decoder Only beating CPA on rank for Frangieh21 is the canonical example). Coverage of calibrated intervals would catch this directly.

## Code URLs

Confirmed: https://github.com/altoslabs/perturbench

Paper: https://arxiv.org/abs/2408.10609

## Verbatim quotes worth keeping

"We introduce a comprehensive framework for modeling single cell transcriptomic responses to perturbations, aimed at standardizing benchmarking in this rapidly evolving field."

"Through extensive evaluation of both published and baseline models across diverse datasets, we highlight the limitations of widely used models, such as mode collapse."

"We also demonstrate the importance of rank metrics which complement traditional model fit measures, such as RMSE, for validating model effectiveness."

"While no single model architecture clearly outperforms others, simpler architectures are generally competitive and scale well with larger datasets."

"Since our benchmarking metrics are defined on a population level, they may not fully capture heterogeneity in the perturbation response among cells. Future work in this area could include using distributional metrics such as maximum mean discrepancy (MMD) or Wasserstein distance to better capture response heterogeneity."

"For a given observed perturbation, the prediction for that perturbation should be more similar than predictions for other perturbations."
