# Yu et al. 2026 - scDFM (ICLR 2026)

**Venue:** ICLR 2026 (published as a conference paper). Authors: Chenglei Yu, Chuanrui Wang, Bangyan Liao, Tailin Wu (corresponding). Affiliations: Zhejiang University and the Department of Artificial Intelligence, School of Engineering, Westlake University. arXiv:2602.07103. OpenReview: QSGanMEcUV. Code: https://github.com/AI4Science-WestlakeU/scDFM.

## Claim

scDFM is a *Distributional* Flow Matching framework (not discrete flow matching) that models the full conditional distribution of perturbed single-cell expression profiles given control states. The headline contribution is conditional flow matching plus a multi-kernel MMD regularizer that aligns the model's terminal distribution with the empirical post-perturbation distribution at the population level, paired with a Perturbation-Aware Differential Transformer (PAD-Transformer) backbone that uses gene-gene co-expression masked attention and differential attention. The headline result is a 19.6% MSE reduction over CellFlow on the Norman additive setting.

## Method

**Flow Matching foundation.** Source x0 is noisy gene expression and target x1 is post-perturbation expression. The model learns a conditional velocity field v_theta(x_t | t, c_x, c_p) where c_x is the control state and c_p is a multi-hot perturbation indicator in {0,1}^d. The state evolves under dX_t/dt = v_theta(X_t | t, c_x, c_p). Training uses a linear interpolation path pi_t(x0,x1) = (1-t)x0 + t x1 in log-normalized expression space. Conditional Flow Matching loss: L_CFM = E[||v_theta(x_t|t,c_x,c_p) - v(x_t|x0,x1,t,c_x,c_p)||_2^2].

**MMD regularizer.** A one-step endpoint estimate is formed at every training step as x_hat_1 = x_t + (1-t) v_theta(x_t|t, c_x, c_p). The terminal sample batch X_hat_1 is compared to ground-truth X_1 via a multi-kernel Gaussian RBF MMD with bandwidths sigma_l = sqrt(s_l * m), s_l in {0.5, 1.0, 2.0, 4.0}, where m is the median of off-diagonal pairwise squared distances in the batch. Final loss is L = L_CFM + lambda L_MMD with lambda = 0.5. Authors choose MMD over KL or Wasserstein because it is "directly sample-based, computationally efficient, and robust under support mismatch."

**PAD-Transformer.** Each layer has three operations: (1) perturbation injection - the perturbation embedding e_p is broadcast and concatenated with the latent then passed through an MLP adapter; (2) self-differential attention on the perturbed latent; (3) cross-differential attention from perturbed latent to control representation h_c. Differential attention follows Ye et al. 2025: A_1 = softmax(Q_1 K_1^T / sqrt(d_h)), A_2 = softmax(Q_2 K_2^T / sqrt(d_h)), then alpha_diff = A_1 - lambda A_2 with learnable lambda. Time embedding t_emb = MLP(SinCos(t)) modulates each attention via adaLN-Zero (Peebles & Xie 2023). Backbone width d=512, L=4 layers, H=8 heads, dropout 0.1.

**Tokenization.** This is a continuous-token model, not discrete tokens. Cell encoder: gene tokens come from a cross-attention encoder E_g with a sparse mask A_hat built from KNN (k=30) over absolute Pearson correlations |Cov(x_i,x_j) / (sigma(x_i) sigma(x_j))| computed on training data. Value embedder E_v maps each scalar log1p-normalized expression to R^d. h_c = E_v(c_x^(S)) + E_g(S), h_v^0 = E_v(x_t^(S)) + E_g(S). Gene-subset sampling at training (size s << G).

**Inference / sample-producing predict().** Initialize x_0 from noise q_0, iterate K=100 Euler steps (Heun optional) of x_{k+1} = x_k + Delta t v_theta(x_k | t_k, c_x, c_p), output x_hat_1 on the selected gene subset. The model produces *samples* of post-perturbation expression for a given (c_x, c_p); no calibration is performed and no coverage guarantee is attached.

**Training corpus.** Training is per-dataset (Norman or ComboSciPlex) from scratch on 5,029 / 5,000-gene HVG sets. No large-scale pretraining. NVIDIA H800 GPUs, Adam with lr 5e-5, cosine to 1e-6, batch size 96, 100,000 optimization steps.

## Datasets and experiments

**Norman.** Perturb-seq with CRISPRa in K562; ~100 single-gene activations and 124 dual activations (Norman et al. 2019). Uses the scFoundation-reprocessed version following Ahlmann-Eltze et al. 2025. Two splits: *additive* (all singles seen; some duals held out) and *holdout* (some singles plus all combinations involving them held out). Four random splits, mean reported.

**ComboSciPlex.** 63,378 single-cell transcriptomes across 32 drug treatment conditions (singles plus pairs) in A549, accessed via pertpy. Held-out drug pairs, single agents in training.

**Baselines.** Control, Additive, scGPT (trained from scratch in their setup), Geneformer (pretrained), GEARS, CPA, STATE, CellFlow.

## Metrics

L2 (mean perturbation distance), MSE, MAE on log-normalized expression; Pearson Delta (correlation matrix discrepancy); Pearson Delta-hat and Pearson Delta-hat_20 (residual correlation versus training perturbation centroid - mitigates control-baseline bias, follows Vinas Torne et al. 2025 Systema framework); DE-Spearman rho (rank correlation on statistically significant DE genes only, p<0.05); Discrimination Score (Roohani et al. 2025 - rank-based pseudobulk retrieval, PDS_p = 1 - (rank_p - 1)/(N - 1)). Notably absent: per-perturbation Wasserstein, energy distance, KS statistic, MMD as evaluation, calibration / coverage.

**Headline numbers.**
- Norman additive: scDFM MSE 0.00315 vs CellFlow 0.00392 (19.6% reduction); MAE 0.02155 vs 0.02207; DS 0.9737 (highest); Pearson Delta-hat_20 0.9260 (highest); DE-Spearman rho 0.5705.
- Norman holdout double: MSE 0.0047, DE-Spearman rho 0.5676, DS 0.9189, Pearson Delta-hat_20 0.8688.
- ComboSciPlex: MSE 0.0028 (lowest), MAE 0.0220 (lowest), DE-Spearman rho 0.8289 (highest), DS 0.8776 (CPA wins DS at 0.8980).

**Ablations (Norman holdout).** Removing MMD causes the largest L2 jump (+10.4%) and the largest DE-Spearman drop (-20.9%), confirming MMD is the load-bearing piece. Removing the gene-gene graph: +10.2% L2, -8.6% DE-Spearman. Removing cross-differential attention: +14.1% L2.

## What we steal for ConfPert

scDFM is the closest direct competitor and is *strictly a point/sample predictor*. It outputs ODE-integrated samples x_hat_1 from a learned conditional velocity field. It does not output sets, intervals, calibrated p-values, or any finite-sample guarantee. The terms "conformal", "coverage", "calibration", "prediction set", "uncertainty quantification", and "miscoverage" do not appear anywhere in the paper, abstract, method, or appendix. The closest authors come to a distributional guarantee is the *training-time* MMD regularizer, which is asymptotic and gives no finite-sample bound on test perturbations. Their evaluation metrics (L2, MSE, DE-Spearman, PDS, Pearson Delta) are all *point* metrics on means, ranks, or correlations - none are distributional discrepancies with calibrated thresholds.

This is exactly the K1 narrative pivot. ConfPert wraps scDFM (and CellFlow, CPA, GEARS, scGen, STATE, Geneformer, scGPT, CellOT) as a base sampler and produces a calibrated set C_alpha(c_p) with a finite-sample guarantee P(D(P_true, P_hat) <= q_hat) >= 1 - alpha for a chosen distributional discrepancy D (Wasserstein, MMD, energy, sliced-Wasserstein). scDFM gives best-effort samples; ConfPert turns those samples into a hypothesis-test-grade interval over distributional discrepancies. Map: scDFM = better generator; ConfPert = first finite-sample coverage layer for the discrepancy itself, agnostic to the generator.

## What we wrap

scDFM is wrappable predictor. The repo at https://github.com/AI4Science-WestlakeU/scDFM is MIT-licensed, has 28 stars, environment.yml is provided, and pretrained checkpoints for Norman and ComboSciPlex are released via Google Drive per the README. Wrapping requires only the inference path (Algorithm 3): given (c_x, c_p, gene subset I), run K=100 Euler steps of PAD-Transformer to draw a sample x_hat_1; repeat to obtain a Monte Carlo batch P_hat. Calibration set: held-out perturbations with paired ground-truth populations. Use them to estimate q_hat for D(P_hat, P_true).

## Failure modes / caveats

(1) Training is per-dataset from scratch, no foundation-scale pretraining; generalization across cell lines is untested. (2) Limitations section explicitly flags *representation and path design* (linear interpolant in log-space may be suboptimal, biological transitions likely non-linear/branching), *scalability to multi-context datasets* (untested on ARC-state and Virtual Cell Challenge data), and *structural priors* (Pearson correlation graph misses non-linear/causal dependencies). (3) DS is a rank-based pseudobulk score, not a distributional metric; high DS does not imply distributional fidelity. (4) ComboSciPlex evaluation drops scGPT/CPA in some columns and never reports CellFlow there, so the cross-domain comparison is partial. (5) MMD with median-bandwidth heuristic is sensitive to batch size 96 and may be biased on small minibatches. (6) Holdout single-perturbation Pearson Delta-hat_20 (0.8116) is much lower than additive (0.9260), confirming the model degrades on truly OOD perturbations - this is exactly where calibrated coverage would matter.

## Code URLs

- GitHub: https://github.com/AI4Science-WestlakeU/scDFM
- arXiv: https://arxiv.org/abs/2602.07103
- OpenReview: https://openreview.net/forum?id=QSGanMEcUV
- Cell-eval (their PDS implementation): https://github.com/ArcInstitute/cell-eval

## Verbatim quotes worth keeping

1. "We present scDFM, a generative framework based on conditional flow matching that models the full distribution of perturbed cells conditioned on control states."
2. "By incorporating an MMD objective, our method aligns perturbed and control populations beyond cell-level correspondences."
3. "FM alone enforces local dynamical consistency and does not guarantee that the terminal distribution of generated cells X_hat statistically aligns with the ground-truth perturbed distribution X."
4. "We choose MMD over KL divergence or Wasserstein distance because it is directly sample-based, computationally efficient, and robust under support mismatch, making it well-suited for high-dimensional single-cell settings."
5. "In the combinatorial setting, it reduces mean squared error by 19.6% relative to the strongest baseline."
6. "Removing MMD regularization causes the sharpest decline, underscoring its critical role in distribution-level fidelity."
7. "Without MMD, generated cells deviate substantially from the ground truth distribution. In contrast, our full model preserves the global geometry and population structure, demonstrating that MMD is essential for stable and biologically consistent predictions."

The paper contains zero references to conformal prediction, coverage, calibration sets, prediction intervals, or finite-sample guarantees. The only "guarantee"-adjacent language is qualitative claims about MMD enforcing population-level fidelity at training time. This is precisely the gap ConfPert fills.
