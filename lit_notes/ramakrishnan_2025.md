# Ramakrishnan et al. 2025 - Distributional Shifts beyond Mean

**Venue:** arXiv preprint arXiv:2507.02980, submitted 1 July 2025 (q-bio.GN). Authors: Kalyan Ramakrishnan, Jonathan G. Hedley, Sisi Qu, Puneet K. Dokania, Philip H. S. Torr, Cesar A. Prada-Medina, Julien Fauqueur, Kaspar Martens. Title: "Modeling Gene Expression Distributional Shifts for Unseen Genetic Perturbations." Not peer-reviewed at time of writing. Affiliations span Oxford VGG and GSK AI/ML.

## Claim
Existing perturbation response models predict only mean expression changes and ignore the stochasticity of single-cell expression. Ramakrishnan et al. argue that gene expression is best modeled as a distribution across cells, that real Perturb-seq perturbations induce shifts in variance, skewness, and modality (not just mean), and that no prior method jointly handles distributional prediction and generalization to unseen perturbations. They propose a histogram-output model conditioned on LLM-derived gene embeddings as the first method to do both.

Exact framing of the gap (verbatim, abstract): "Existing methods only predict changes in the mean expression, overlooking stochasticity inherent in single-cell data. In contrast, we offer a more realistic view of cellular responses by modeling expression distributions."

Introduction framing (verbatim): "Recent techniques such as Perturb-seq ... enable expression profiling across thousands of single cells following genetic perturbations (gene knockouts), revealing shifts in mean expression, as well as in variance, skewness, and modality (Fig. 1)."

## Method
Per-gene marginal distribution model. For each gene g, the expression range is discretized into B fixed-width bins to form a histogram. Default B = 15. The expression range is r_g = [min(x_g) - eps_g, max(x_g) + eps_g] with eps_g = (max(x_g) - min(x_g)) / (2(B-1)).

Perturbation embedding: gene-description embeddings from GPT-3.5 concatenated with ProtT5 protein-sequence embeddings (frozen).

Architecture: an MLP Delta_theta maps the perturbation embedding e_p to per-gene-bin shifts. Predicted histogram is the softmax of the control log-histogram plus the shift: h_theta(p) = softmax(log h_ctrl + Delta_theta(e_p)). This is an additive log-shift on top of the control distribution, ensuring identity behavior under no perturbation.

Training objective: convex combination of closed-form Wasserstein-1 distance between predicted and empirical perturbed histograms and an MSE term on the mean:
L_total = lambda_wass * L_wass + lambda_mse * L_mse, with lambda_wass = 0.75, lambda_mse = 0.25.

Generalization mechanism: because perturbations are represented as continuous LLM embeddings, the model produces histograms for genes not seen during training. The shift module Delta_theta is shared across perturbations.

## Datasets / experiments
- Norman et al. 2019 Perturb-seq: 102 perturbations.
- Replogle et al. 2022 K562 essential-gene Perturb-seq: 1076 perturbations.
- Replogle et al. 2022 RPE1 essential-gene Perturb-seq: 1517 perturbations.

Evaluation protocol: 9-fold cross-validation across perturbations, with held-out perturbations functioning as unseen targets. The training cost claim is roughly 3x faster than competing distributional baselines.

## Metrics
- Negative log-likelihood (NLL) of held-out cells under the predicted histogram.
- Relative Mean Absolute Error per moment: RMAE(tau_hat, tau) = sum_g |tau_hat_g - tau_g| / sum_g |tau_g|, computed for predicted vs observed mean, variance, skewness, and excess kurtosis. Moments are computed in closed form from histogram bin probabilities (Appendix A).
- Pearson correlation between predicted and observed mean shifts (for comparison to mean-only baselines).

Distributional baselines:
- Ctrl Hist (control histogram applied to all perturbations).
- Non-Ctrl Hist (pooled perturbed-cell histogram).
- TG (truncated Gaussian fit per gene-perturbation pair).
- ZITG (zero-inflated truncated Gaussian).

Mean-based baselines for the moment-1 comparison:
- Non-Ctrl Mean.
- GEARS (Roohani et al. 2024).
- GP from Martens et al. 2024.
- MLP+Mean (authors' own MLP head trained on mean only).
- scGPT and scFoundation are referenced but excluded from quantitative evaluation, citing Ahlmann-Eltze et al. 2024.

Reported headline finding: the histogram model wins on variance, skewness, and kurtosis RMAE while remaining competitive with mean-only methods on mean RMAE and Pearson correlation, at lower training cost.

## What we steal for ConfPert
Cite as the third "field is in crisis" co-motivator alongside Ahlmann-Eltze 2025 and Csendes 2025. Ahlmann-Eltze shows that current foundation models do not beat linear baselines on mean prediction. Csendes shows that perturbation responses are heterogeneous in ways that mean predictions wash out. Ramakrishnan supplies the third leg: even when mean prediction works, the distribution still carries information about variance, skewness, and modality that every mean-only method discards by construction. The triad supports the ConfPert thesis that distributional fidelity beyond the mean is the missing axis. ConfPert's calibrated-distributional framework is the constructive response: rather than swap one point estimate for another, wrap nine existing predictors and produce calibrated distributional prediction sets with finite-sample coverage guarantees.

Specific technical pieces worth borrowing for framing and ablations:
- The per-gene marginal histogram view as a tractable distributional target. ConfPert can use it as a strawman comparator showing that uncalibrated histogram fits do not give coverage.
- The four-moment RMAE (mean, variance, skewness, excess kurtosis) as an evaluation axis. Useful as a complement to Wasserstein-1 and energy distance when arguing that ConfPert prediction sets respect higher-order structure.
- The closed-form W1 between histograms as a cheap training-time proxy and as an evaluation metric.
- The control-anchored additive log-shift parameterization is a clean inductive bias worth referencing when discussing why naive distributional outputs collapse to the mean.

## What we wrap
Not a model we wrap. Positioning citation only. Ramakrishnan's predictor is one of the natural distributional baselines, but it is per-gene marginal, treats genes independently, and is constrained to the training expression range, all of which make it a weaker target for wrapping than CPA, scGen, CellOT, GEARS, scGPT, biolord, sVAE+, STATE, or scDFM. Use it in the related-work section and in the motivation, not in the predictor zoo.

## Failure modes / caveats
From the conclusion (verbatim): "First, predictions are constrained by the expression range observed during training due to the histogram representation, which limits extrapolation and requires truncated baselines for fair comparison. Second, the method treats genes independently by modeling marginal distributions but not the regulatory dependencies between genes."

Other caveats:
- The model produces marginal per-gene histograms, not a joint distribution over the gene panel. Co-expression and regulatory dependencies are unmodeled, which is precisely the structure many downstream uses (pathway shifts, gene-gene epistasis) need.
- The histogram representation is bounded by the training expression range; out-of-range observations cannot be assigned positive mass.
- Default B = 15 bins is coarse for high-dynamic-range genes.
- Discrete representation; continuous extensions (KDE, flows) are listed as future work.
- LLM embeddings (GPT-3.5 description text, ProtT5 sequence) carry their own biases and saturation behavior; performance for genes with sparse literature or unusual sequences is not characterized.
- The paper does not provide finite-sample coverage guarantees on its histogram predictions. Calibration of the predicted distributions is not reported. This is the gap ConfPert exploits.
- No theoretical analysis of when distributional extrapolation succeeds.
- Evaluation is restricted to single-gene CRISPRi perturbations on Norman and Replogle; combinatorial perturbations and chemical perturbations are not evaluated.

## Code URLs
GitHub: https://github.com/Kalyan0821/LLMHistPert (linked from the paper). No model checkpoints are advertised in the abstract or paper page; assume training-from-scratch on the public Norman and Replogle datasets.

## Verbatim quotes worth keeping
1. "Existing methods only predict changes in the mean expression, overlooking stochasticity inherent in single-cell data. In contrast, we offer a more realistic view of cellular responses by modeling expression distributions."
2. "Gene expression is stochastic - even in genetically identical cells, levels fluctuate due to transcriptional noise, regulatory variation, and measurement uncertainty (Paulsson, 2005; Raj and van Oudenaarden, 2008). Thus, gene expression is best modeled as a distribution across cells, a crucial aspect many methods overlook."
3. "Recent techniques such as Perturb-seq ... reveal shifts in mean expression, as well as in variance, skewness, and modality."
4. "Despite their scale, these models only predict shifts in the mean and ignore distributional changes. Moreover, Ahlmann-Eltze et al. (2024) show that these sophisticated models fail to outperform simple linear baselines when evaluated on out-of-distribution (OOD) perturbations."
5. "We address the above limitations by (i) predicting gene expression distributions and (ii) generalizing to unseen perturbations. To our knowledge, no existing method jointly tackles both."
6. "First, predictions are constrained by the expression range observed during training due to the histogram representation, which limits extrapolation and requires truncated baselines for fair comparison. Second, the method treats genes independently by modeling marginal distributions but not the regulatory dependencies between genes."
