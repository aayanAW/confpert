# Ahlmann-Eltze, Huber, Anders 2025 - DL Perturbation Response vs Linear Baselines

**Venue:** Nature Methods 2025. DOI 10.1038/s41592-025-02772-6. PMID 40759747. PMC12328236. bioRxiv 10.1101/2024.09.16.613342 (v1 Sep 2024 through v5).

## Claim
Five single-cell foundation models (scGPT, scFoundation, scBERT, Geneformer, UCE) and two task-specific deep learning models (GEARS, CPA) fail to outperform a low-rank bilinear ridge baseline on Perturb-seq for predicting unseen single-gene perturbations. None outperform the additive baseline `y_AB = y_A + y_B - y_empty` on Norman doubles. The qualitative diagnostic in Extended Data Fig. 7 documents across-perturbation variance collapse: scGPT, UCE, scBERT predict near-constant outputs across perturbations for most genes, while GEARS and scFoundation predict perturbation-conditional outputs whose variance is compressed relative to ground truth. The authors frame this as a benchmarking failure of the field: foundation-model pretraining on observational scRNA-seq does not transfer to causal perturbation response prediction, and the bar to clear is much lower than community claims suggest.

## Method
Bilinear ridge baseline (eq. 1):
`argmin_W || Y_train - (G W P^T + b) ||_F^2 + lambda * ||W||_F^2`

where `G` is a K-dimensional gene embedding obtained as the top-K principal components of the pseudo-bulked training expression matrix `Y_train` (rows = readout genes), `P` is an L-dimensional one-hot perturbation indicator matrix (rows = perturbations), `W` is the K x L coefficient matrix, and `b` is the row mean of `Y_train` (per-gene control mean). `K = 10` PCs and `lambda = 0.1` for numerical stability. Equation 3 gives the closed-form ridge solution via the normal equations. At test time, prediction for an unseen perturbation `p_test` is `y_hat = G W p_test + b`; for foundation-model embedding evaluations, `G` is replaced by the model's gene embeddings and the same ridge fit is applied on top.

Additive baseline for Norman-doubles (eq. 2):
`y_AB = y_A + y_B - y_empty`

where `y_A`, `y_B` are observed pseudo-bulk means under the corresponding single perturbations and `y_empty` is the unperturbed control mean. This baseline never sees double-perturbation training data and requires no model fitting; it is the null for genetic-interaction detection.

Interaction discovery: Efron empirical-null mixture model (locfdr v1.1-8) on the residuals between observed double-perturbation expression and additive expectation, controlled at 5 percent false discovery rate. Pseudo-bulking is computed per perturbation (mean across cells of the same gRNA assignment) before all evaluations.

## Datasets / experiments
- Norman et al. 2019: K562 CRISPRa, 100 single-gene perturbations and 124 double-gene perturbations across 225 conditions, 81,143 cells, 19,264 genes. Used for both unseen-single and unseen-double splits.
- Replogle K562 2022: 398 perturbations after intersection with foundation-model gene vocabularies. CRISPRi. Used for unseen-single evaluation.
- Replogle RPE1 2022: 629 perturbations. CRISPRi. Used for unseen-single evaluation including transfer-learning experiments.
- Adamson et al. 2016: K562, 73 perturbations. Smaller secondary benchmark.

All evaluations restricted to the top-1,000 most-highly-expressed genes (computed on the unperturbed pseudo-bulk) to avoid measurement noise dominating L2.

## Metrics
- L2 distance on pseudo-bulk means: `sum_g (y_hat_g - y_g)^2` over the top-1,000 genes.
- Pearson-delta: `cor(y_hat - y_empty, y - y_empty)`. Distributional-flavored in the sense that it operates on the change vector rather than raw expression, but it is still a per-perturbation mean-on-mean correlation; it does not exercise per-cell predictive distributions or coverage.
- True-positive rate and false-discovery proportion for Norman-doubles interaction discovery, against the additive null at FDR 5 percent.
- Across-perturbation variance collapse is documented qualitatively in Extended Data Fig. 7 by per-gene scatter of predicted vs observed perturbation effects with marginal density.

## What we steal for ConfPert
Ahlmann-Eltze 2025 is the existential motivator for ConfPert. The headline finding "DL doesn't beat linear baselines on means" is correct and non-trivial, but it is mean-only. The implication for our framing is direct: the field needs to shift from capacity contests on point estimates to calibrated distributional metrics with provable coverage, which is exactly ConfPert's pitch. We cite Extended Data Fig. 7 specifically: it documents ACROSS-perturbation variance collapse (the model predicts near-identical means across different perturbations), while ConfPert documents WITHIN-perturbation per-cell collapse (the model predicts near-deterministic distributions for any single perturbation). The two failures stack: a model that collapses across perturbations and within perturbations is a constant predictor with miscalibrated nominal coverage, and the only way to detect this is calibration auditing rather than mean MSE. We adopt their bilinear ridge as ConfPert's mean-only reference baseline so reviewers cannot accuse us of strawmanning DL by comparing only to neural baselines, and so our calibration story stands on top of an Ahlmann-Eltze-grade foundation. We also adopt their Pearson-delta as a comparator while making explicit it is still mean-based.

## What we wrap
Not a model. A critique paper. The bilinear ridge is fully specified by equations 1 and 3 plus K=10, lambda=0.1, and is implemented in their public R/Python codebase. We re-implement in our pipeline rather than depending on their package, since their code is paper-bound rather than library-bound.

## Failure modes / caveats
Their evaluation is strictly mean-on-mean at the pseudo-bulk level. Their critique is correct but partial. Three gaps ConfPert closes:
1. They evaluate point predictions, not predictive distributions. A model can match pseudo-bulk means while being arbitrarily badly calibrated at the single-cell level.
2. Their "variance collapse" diagnostic is qualitative (visual scatter in ED Fig. 7), not a quantitative calibration test with coverage targets.
3. They do not address distribution shift between training and test perturbations beyond an unseen-single split; conformal coverage under perturbation-induced shift is unaddressed.
Calibration is the missing axis. ConfPert is positioned as the constructive complement: they show DL fails on means, we show DL fails on coverage, and we provide the prediction-set machinery to make uncertainty meaningful.

## Code URLs
- GitHub: https://github.com/const-ae/linear_perturbation_prediction-Paper
- Zenodo archive: https://doi.org/10.5281/zenodo.14832393
- Results and predictions RDS files (Zenodo 16092690, includes single_perturbation_results_predictions.RDS at 5.9 GB and double_perturbation_results_predictions.RDS at 1.6 GB) for direct re-evaluation under our calibration metrics.
The bilinear ridge is fully specified by equations 1, 3 plus the two hyperparameters, so re-implementation is trivial and we do not need their codebase as a runtime dependency.

## Verbatim quotes worth keeping
1. "For most genes, the predictions of scGPT, UCE and scBERT did not vary across perturbations, and those of GEARS and scFoundation varied considerably less than the ground truth." This is THE ED Fig. 7 quote and the load-bearing diagnostic for the variance-collapse argument.
2. "None outperformed the baselines, which highlights the importance of critical benchmarking in directing and evaluating method development." Abstract.
3. "Deep learning is effective in many areas of single-cell omics. However, prediction of perturbation effects still remains an open challenge." Discussion.
4. "We hypothesize that the poor performance is partially because the pre-training data is observational." (v1 abstract; retained in revisions.) Sets up the causal-vs-observational gap that ConfPert's calibration framing further sharpens.
5. The additive baseline "has never seen the double perturbation training data and does not require any computationally expensive model training/fine-tuning." Frames the embarrassment of DL models being beaten by a closed-form expression.
6. Their bilinear ridge uses "a ridge penalty of lambda = 0.1 for numerical stability" with K=10 PCs of the pseudo-bulk training expression. Method, hyperparameters fixed not tuned per dataset.
7. "Despite the considerable computational expense of training foundation models, they have not yet achieved generalizable representations of cellular state sufficient to predict the outcome of unperformed experiments." Conclusion.

ConfPert positioning sentence to seed in our intro: Ahlmann-Eltze 2025 establishes that contemporary perturbation-response DL underperforms a 10-PC bilinear ridge on pseudo-bulk means; we extend the diagnosis from means to calibrated predictive distributions, replacing point-error contests with conformal coverage audits.
