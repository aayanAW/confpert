# Yuan et al. 2026 - PerturbDiff

**Venue:** ICML 2026 submission. arXiv:2602.19685.

Authors: Xinyu Yuan, Xixian Liu, Yashi Zhang, Zuobai Zhang, Hongyu Guo, Jian Tang. Mila / U. Montreal / U. Ottawa / NRC Canada / HEC Montreal / CIFAR. Submitted Feb 23, 2026.

## Claim
PerturbDiff is a diffusion framework that operates over entire cell distributions (not single cells) by embedding empirical cell distributions as kernel mean embeddings in an RKHS, then defining a forward/reverse diffusion process directly in that Hilbert space. The MMD-based denoising objective falls out of the squared RKHS distance between kernel mean embeddings. Claim: SOTA on 14 metrics across PBMC, Tahoe100M, Replogle, with strong gains on differential-expression metrics and improved generalization to unseen perturbations driven by unobserved latent factors (microenvironment, batch effects).

## Method
**Distribution-valued random variable.** Standard methods learn a single fixed perturbed distribution P_{c,tau} conditioned on context c and perturbation tau. PerturbDiff instead treats the perturbed distribution itself as a random variable D_{c,tau} taking values in P(X) (space of distributions over cell space X subset R^|G|). They model the conditional law D_{c,tau} ~ F_theta(. | D_c, c, tau), where D_c is the analogous control distribution-valued random variable.

**RKHS embedding.** Each empirical batch B_pert = {x_1, ..., x_m} is mapped to its kernel mean embedding mu_P := E_{Z~P}[k(Z, .)] in H_k. The kernel k: X x X -> R is positive-definite. The paper does not pin down a specific kernel choice in the main text excerpts I read; standard MMD-with-RBF is implied by the cited Borgwardt et al. 2006 and by Gretton et al. 2012. (Confirm in appendix B if available.)

**Forward diffusion (Defn 4.2).** Variance schedule {beta_t} subset (0,1), self-adjoint trace-class psd covariance operator C: H_k -> H_k. Forward Markov chain in H_k:
  mu_t = sqrt(1 - beta_t) mu_{t-1} + sqrt(beta_t) Xi_t,  Xi_t ~ N_{H_k}(0, C).
This yields the closed-form marginal mu_t | mu_0 ~ N_{H_k}(sqrt(alpha_t) mu_0, (1 - alpha_t) C), with alpha_t := prod (1 - beta_s). This extends DDPM into H_k.

**Reverse / variational objective.** True reverse conditional admits a closed-form Gaussian measure whose mean is linear in mu_t and mu_0. Sampling from it is intractable since mu_0 is unknown at generation time, so they parameterize a learnable mean function mu_theta and minimize the RKHS-norm denoising loss L_t prop ||mu_0 - mu_theta(mu_t, t)||^2_{H_k}. Conditional version with classifier-free guidance: L_t prop ||mu_0 - mu_theta(mu_t, t, mu_c, c, tau)||^2_{H_k}.

**Distribution-aware MMD loss.** Key identity (Rem 4.1, kernel geometry property iii): ||mu_P - mu_Q||^2_{H_k} = MMD^2_k(P, Q). The loss therefore becomes the MMD between the predicted batch B_theta = {x^theta_1, ..., x^theta_m} from a generator network and the true perturbed batch:
  ||mu_0 - mu_theta(mu_t, t, mu_c, c, tau)||^2_{H_k} = MMD^2_k(P~_pert, P~_theta).

**Tractable Gaussian-noise approximation.** Direct sampling from H_k-valued Gaussian Xi_t is intractable. They inject additive Gaussian noise independently to each cell in the original space and recompute the empirical kernel mean embedding. Under a first-order linearization of the kernel feature map this induces an H_k-valued Gaussian random element with explicit covariance operator (Prop B.14).

**Multi-scale hybrid loss (Eqn 8).** L_total = MMD^2_k(P~_pert, P~_theta) + (1/m) sum_j ||x_j - x^theta_j||^2. Empirically MMD dominates; MSE acts as regularization on the global batch centroid.

**Per-perturbation conditioning.** Implemented via classifier-free guidance: training mixes conditional and unconditional inputs with prob p_drop, sampling uses x_theta(x_t, t, eta) = (1+w) x_theta(x_t, t, eta) - w x_theta(x_t, t, empty). Architecture: Multi-Modal Diffusion Transformer (MM-DiT, Esser et al. 2024) with two parallel token streams (control B_ctrl, perturbed B_pert) interacting via joint attention. (c, tau) and timestep t injected via AdaLN-Zero (Peebles and Xie 2023). Conditioning tokens: cell type, batch, perturbation.

**Sampling.** Since H_k is intractable, deterministic DDIM sampling is performed directly in cell space (consistent with training procedure).

**Marginal pretraining.** D_c ~ F_theta(. | c) on 61M unperturbed cells from CellxGene before finetuning on perturbation data. Improves zero-shot R^2 and low-data adaptation.

## Datasets / experiments
Three perturbation benchmarks (Tab 1): PBMC (9.7M cells, 90 perturbations, 18 cell types, 12 batches; signaling/cytokines), Tahoe100M (101.2M cells, 1137 perturbations, 50 CTs, 14 batches; drug), Replogle (0.6M cells, 2023 perturbations, 4 CTs, 56 batches; CRISPR genetic). Pretraining adds CellxGene 60.9M cells, 662 CTs, 10887 batches. Predictions over 2000 highly variable genes; merge-then-select gene unification yields 12626 informative genes (>98% overlap PBMC/Tahoe100M, >45% Replogle/CellxGene).

## Metrics
Cell-Eval (STATE) framework plus R^2 from CellFlow. 14 metrics in two families:
- Average expression accuracy: R^2, PDCorr (Pearson Delta), MAE, MSE, PDS_{L1}, PDS_{L2}, PDS_{cos} (Perturbation Discrimination Scores).
- Biologically meaningful DE patterns: DEOver (DE overlap), DEPrec (DE precision), DirAgr (direction agreement), LFCSpear (log-fold-change Spearman), AUROC, AUPRC, ES (effect size correlation).

DE genes identified by Wilcoxon rank-sum test with adjusted p < 0.05. PerturbDiff (From Scratch) wins on PBMC and Tahoe100M across nearly all metrics; second best on Replogle. PerturbDiff (Finetuned) closes the gap on Replogle and matches STATE there. Per-condition win rates over STATE on DE metrics: 96-100% for AUROC/AUPRC/DEPrec on PBMC; 99% AUROC on Tahoe100M. Ablation: removing MMD term degrades DE metrics substantially (Fig 9). Scaling: medium 114M-param model best on PBMC; 239M model overfits at large compute (Fig 7). Low-data: pretraining gives substantial gains at 1% and 5% sample ratios on PBMC (Fig 8).

Baselines: Mean, Linear (Ahlmann-Eltze 2024), CPA (Lotfollahi 2023), STATE (Adduri 2025), CellFlow (Klein 2025), Squidiff (He 2025). Unlasting excluded (no public code).

## What we steal for ConfPert
PerturbDiff is one of the saturated-subfield wrappable predictors. We do NOT compete with their RKHS functional-diffusion construction. We wrap it for finite-sample coverage. Their RKHS framing is directly relevant to ConfPert's MMD-RBF discrepancy choice: by their Rem 4.1 (kernel geometry), squared RKHS distance equals MMD^2_k, so any conformal score we build over kernel mean embeddings of predicted vs. observed batches is a calibrated MMD score. This means our discrepancy-on-distributions calibration target lines up with the natural distance their model already optimizes. Useful: their proof that adding pointwise Gaussian noise at the cell level yields an H_k-valued Gaussian with explicit covariance (Prop B.14) gives a closed-form way to characterize predictive uncertainty bands in feature space if we want to certify coverage in H_k rather than R^|G|.

## What we wrap
PerturbDiff (From Scratch) and PerturbDiff (Finetuned). Confirm public-code status: paper states "code and data will be made publicly available" with project page; references repo DeepGraphLearning/PerturbDiff. Status at submission likely pending; check before integration. Wrapping interface: model outputs B_theta = {x^theta_1, ..., x^theta_m} per (c, tau). Standard split-conformal over per-batch MMD discrepancy or per-gene residuals will work directly.

## Failure modes / caveats
- Kernel choice not pinned down in main text; sensitivity to kernel bandwidth not reported. RBF MMD is well-known to be bandwidth-sensitive; ConfPert calibration must fix and report bandwidth.
- H_k is infinite-dimensional; all training operates on empirical kernel mean embeddings via finite batches, with first-order linearization for noise injection (Prop B.14). Approximation error not quantified.
- Sampling is performed in cell space via DDIM, not in H_k; reverse process closed-form requires unknown mu_0 so they fall back to learned mu_theta.
- Scaling is non-monotonic: 239M model overfits, indicating capacity-data mismatch on Replogle scale.
- MSE term in hybrid loss adds only "simple regularization on the global batch centroid"; if MMD term collapses (e.g., sparse expression with >95% zeros in PBMC, >60% in Replogle as they note), MSE saturates on shared zeros.
- OOD generalization: claims gains on unseen perturbations but per-perturbation Fig 4 only spans held-out conditions within the same datasets (62 PBMC, 735 Tahoe100M); cross-dataset OOD not tested in main text.
- Computational cost in high-dim function space implicitly bounded by MM-DiT compute (57M-410M params, 10^16-10^18 FLOPs).

## Code URLs
Paper states code release is planned: project page linked in abstract, repository referenced as `github.com/DeepGraphLearning/PerturbDiff`. Verify availability before integration.

## Verbatim quotes worth keeping
On distributional variability (relevant to ConfPert bimodality coefficient): "responses vary systematically due to unobservable latent factors such as microenvironmental fluctuations and complex batch effects, forming a manifold of possible distributions for the same observed conditions." Also: "cell samples collected under the same (c, tau) can be influenced by unobservable latent factors ... These latents induce variability at the distribution level, corresponding to a family of plausible response distributions associated with the same (c, tau)."

On MMD as natural objective: "MMD is not externally imposed, but arises naturally from the equivalence property of squared RKHS distance between kernel mean embeddings, which coincides with MMD by construction."

On STATE overconfidence on DE (relevant to bimodality / heteroskedasticity in differential expression): "STATE predicts most genes as DE by assigning large -log_10(p_adj) values, including many true non-DE genes. This systematic overconfidence suggests that STATE fails to learn the perturbation driven DE patterns." This is an explicit miscalibration finding in a competing baseline and supports ConfPert's framing that point-prediction models in this subfield ship without coverage guarantees.

On distribution-as-random-variable framing: "instead of assuming a single perturbed distribution P_{c,tau}, we consider the perturbed distribution to be a random variable D_{c,tau} taking values in P(X). The realizations of D_{c,tau} are cell distributions induced by unobservable latent factors."

No explicit "bimodality" terminology, but the manifold-of-distributions framing is the closest formal hook for ConfPert's bimodality-coefficient-match metric: PerturbDiff explicitly models multi-modal response distributions over cell distributions, which is consistent with bimodal/multimodal marginals at the gene level.
