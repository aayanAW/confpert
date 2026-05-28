# Chi et al. 2025 - Unlasting

**Venue:** arXiv:2506.21107 (v2, June 2025; renamed "Doloris" in v3, April 2026). Authors: Changxi Chi, Jun Xia, Yufei Huang, Jingbo Zhou, Siyuan Li, Yunfan Liu, Chang Yu, Stan Z. Li. Affiliations: Zhejiang University; AI Lab, Westlake University; HKUST (Guangzhou).

## Claim
Unlasting frames single-cell perturbation estimation as an unpaired distribution-translation problem and proposes a Dual Diffusion Implicit Bridges (DDIB) framework that learns control and perturbed distributions separately, sharing a Gaussian latent space. The paper explicitly motivates a shift away from conditional-mean evaluation because "gene expression under the same perturbation often varies significantly across cells, frequently exhibiting a bimodal distribution that reflects intrinsic heterogeneity." It introduces distribution-aware metrics (E-distance, EMD) that respect bimodal response structure and reports state-of-the-art performance on Adamson, sci-Plex3, and Norman.

## Method
Architecture: two conditional diffusion models, source and target, sharing a Gaussian latent. Source learns the unperturbed distribution; target learns the perturbed distribution. Direct x_0 prediction parameterization: x_hat_0 = x_hat_theta(x_t, t, c, mu_c, sigma_c, P) where c is cell type, P is perturbation, and (mu_c, sigma_c) are control-group expectation and standard deviation. To preserve cell heterogeneity, the target model conditions on a noisy control sample ctrl_noisy = mu + sigma * eps rather than the bare mean, because "using only the expectations mu of unperturbed group gene expression is unreasonable, as it disregards cell heterogeneity."

GRN block: three condition-specific embeddings (G_ctrl, G_mole, G_gene) for unperturbed, molecular, and gene-knockout perturbations. Knockout encoded by zeroing the perturbed gene's embedding via a mask M_bar. Message passing uses multi-head graph attention over the GRN: F^(l+1) = (1/H) sum_h GAT^h(A, F^l). Gene-wise output head: F_GW,i = W_i (elementwise) g_tilde_P^i + b_i.

Loss: masked MSE on x_0 prediction, with mask M zeroing silent genes so the diffusion focuses on expressed genes. Equation 13: L = E[ ||M (elementwise) (x_0 - x_hat_theta(x_t, t, c, mu_c, sigma_c, P))||^2 / sum_i M_i ].

Sampling: deterministic DDIB forward then reverse. Forward maps a control cell to latent with the source ODE solver (t: 0->1); reverse maps latent to perturbed cell with the target ODE solver conditioned on (c, P) (t: 1->0). Final output is masked by a learned silent-gene mask M_hat_{c,P} and rescaled by x_max: x_hat_0 = (M_hat_{c,P} (elementwise) x^t) * x_max.

## Datasets / experiments
Adamson (CRISPR knockouts): 87 single-gene perturbations, single cell type, 5,000 genes. sci-Plex3 (chemical): 187 drugs at four dosage levels across three cell types, 2,000 genes. Norman: double gene knockouts. Held-out scenarios include double-gene perturbations and OOD drugs. Training and evaluation on one Nvidia A100.

## Metrics
Energy Distance (E-distance): D_E(X,Y) = (2/(nm)) sum_ij ||X_i - Y_j||_2 - (1/n^2) sum_ij ||X_i - X_j||_2 - (1/m^2) sum_ij ||Y_i - Y_j||_2. Earth Mover's Distance (EMD): D_EMD(X,Y) = (1/|N|) sum_{j in N} EMD(X_{:,j}, Y_{:,j}), per-gene 1D Wasserstein averaged over genes. Reported on All-genes and DE20 (top-20 differentially expressed). Concrete numbers: Adamson E-dist 0.8426 +/- 0.0446 (Unlasting) vs 1.9939 (ScLambda), 0.8705 (GRAPE); Adamson EMD 0.0397 vs 0.0906 (ScLambda), 0.0444 (GRAPE). sci-Plex3 E-dist 0.7034 vs 0.7847 (chemCPA); EMD 0.0255 vs 0.0838. Double-gene (Norman) E-dist 0.7040 vs 0.7722 (GRAPE), EMD 0.0145 vs 0.0169. OOD drug E-dist 0.7371 vs 0.8861 (chemCPA), EMD 0.0355 vs 0.0959.

## What we steal for ConfPert
ConfPert positioning: Unlasting is the closest co-motivator on the bimodality point. Their paper documents the same phenomenon ConfPert is built around (response distributions are bimodal, so conditional means are misleading) and they propose distributional metrics (E-distance, EMD) as the corrective at evaluation time. ConfPert's contribution is orthogonal and complementary: instead of replacing the predictor's loss target or the evaluation metric, ConfPert wraps any base perturbation predictor (Unlasting included) with finite-sample, distribution-free coverage guarantees over the bimodal response. Cite Unlasting as a co-motivator for "conditional mean is the wrong target under bimodal response heterogeneity," then differentiate: Unlasting answers with a better generative model and better metrics; ConfPert answers with calibrated prediction sets / bands that cover the multimodal response with user-specified probability. Borrow their Adamson/sci-Plex3/Norman split conventions and their E-distance/EMD reporting so ConfPert's tables are directly comparable. Use their Figure 5 framing of bimodal DE genes as a motivating figure for ConfPert with proper attribution.

## What we wrap
Unlasting is wrappable in principle: the inference path is a deterministic ODE solver pair, and DDIB sampling produces a per-cell perturbed expression vector that can be treated as the base predictor f(x, P) inside a split-conformal pipeline. Practical wrapping is currently blocked because v2 contains no public code link. v3 ("Doloris") points to https://github.com/ChangxiChi/Doloris, which is the same architecture; verify that the released checkpoint covers Adamson/sci-Plex3/Norman before committing to wrap. If the code release is incomplete, ConfPert can still cite Unlasting as a baseline to be wrapped once weights are available.

## Failure modes / caveats
No dedicated limitations section. Inferred risks: (1) GRN dependence -- performance is conditioned on the quality and coverage of the input gene regulatory network; perturbations of genes outside the GRN are degraded. (2) The masked MSE loss on expressed genes may underfit dropout-driven distributional structure even though the silent-gene mask model is separate. (3) Conditional-mean style supervision per gene is retained inside the diffusion target despite the bimodality argument; the bimodality argument is honored at evaluation, not at training. (4) DDIB assumes a shared latent prior across source and target; large covariate shifts between control and perturbed populations (e.g., cell-type composition shifts) may break the alignment assumption. (5) Reported metrics are means with standard deviations across seeds; no per-perturbation calibration or coverage analysis, which is exactly the gap ConfPert addresses.

## Code URLs
v2 ("Unlasting"): no public code link in the manuscript. v3 ("Doloris", same paper renamed in April 2026): https://github.com/ChangxiChi/Doloris (referenced in PDF metadata; treat as the operative repo).

## Verbatim quotes worth keeping
1. "Moreover, gene expression under the same perturbation often varies significantly across cells, frequently exhibiting a bimodal distribution that reflects intrinsic heterogeneity."
2. "In this section, we observe strong heterogeneity in single-cell data, where many differentially expressed (DE) genes exhibit bimodal distributions under the same condition."
3. "This renders expectation-based metrics unreliable, as they may obscure true expression patterns. To address this, we adopt distribution-aware evaluation metrics: Energy Distance (E-distance) and Earth Mover's Distance (EMD)."
4. "Cells observed under the same experimental conditions exhibit a bimodal distribution for many genes." (Figure 5 caption)
5. "Furthermore, using only the expectations mu of unperturbed group gene expression is unreasonable, as it disregards cell heterogeneity."
6. "Due to the noticeable heterogeneity among cells under identical conditions, including bimodal gene expression in some cases, conventional metrics may fail to fully capture the distributional characteristics."
7. "The results on publicly available datasets show that our model effectively captures the diversity of single-cell perturbations and achieves state-of-the-art performance."
