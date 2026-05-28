# Xu, Chen, Sun, Venkitasubramaniam, Xie 2025 - Wasserstein-Regularized Conformal Prediction under Distribution Shift

**Venue:** ICLR 2025 (published as a conference paper). arXiv:2501.13430 (v2, 6 Mar 2025).

Authors: Rui Xu (HKUST Guangzhou), Chao Chen (Harbin Institute of Technology), Yue Sun (Lehigh), Parvathinathan Venkitasubramaniam (Lehigh), Sihong Xie (HKUST Guangzhou, corresponding).

## Claim
A split-conformal regression procedure (WR-CP) that maintains coverage when test distribution shifts from calibration by upper-bounding the coverage gap with a 1-Wasserstein distance between calibration and test conformal score distributions, then minimizing that bound by decomposing it into a covariate-shift component (handled by importance weighting) and a concept-shift component (handled by a Wasserstein regularizer added to the training loss).

## Method
The paper builds the bound in stages:

1. Coverage gap is defined as |F_{P_V}(tau) - F_{Q_V}(tau)| where P_V, Q_V are calibration and test conformal score distributions and tau is the empirical (1-alpha) quantile of calibration scores. By Kolmogorov distance and Proposition 1 (Ross 2011), Coverage gap <= sqrt(2 L W(P_V, Q_V)) where L is a Lebesgue density bound on P_V and W is the 1-Wasserstein distance (Eq. 6).

2. Wasserstein decomposition (Eq. 7): introduce the pushforward Q_{V,s_P} = s_P # Q_X (apply the calibration-side score function s_P(x) = |h(x) - f_P(x)| to test covariates). Then by triangle inequality W(P_V, Q_V) <= W(P_V, Q_{V,s_P}) + W(Q_{V,s_P}, Q_V). The first term is purely covariate shift; the second is purely concept shift.

3. Theorem 1 (extension of Aolaritei et al. 2022) rewrites Wasserstein between two pushforwards under different functions as an infimum over couplings of the underlying spaces. Theorem 2 gives W(mu_f, nu_f) <= kappa W(mu, nu) for kappa-Lipschitz f. Combined: W(P_V, Q_V) <= kappa W(P_X, Q_X) + eta W(Q_{Y,f_P}, Q_Y), so Coverage gap <= sqrt(2 L (kappa W(P_X, Q_X) + eta W(Q_{Y,f_P}, Q_Y))) (Eq. 13). kappa is the Lipschitz constant of s_P; eta = max |s_P(x_1) - s_Q(x_2)| / |f_P(x_1) - f_Q(x_2)|.

4. Theorem 3 gives a finite-sample empirical bound: with n calibration and m test samples, W(mu, nu) <= W(mu_hat_n, nu_hat_m) + lambda_mu n^{-1/sigma_mu} + lambda_nu m^{-1/sigma_nu} + t_mu + t_nu with high probability (Eq. 16), where sigma > d_W is the upper Wasserstein dimension (Weed and Bach 2019).

5. Multi-source application (Section 5): if Q_{XY} is an unknown convex mixture of k training distributions D_{XY}^{(i)}, Theorem 4 gives W(Q_{V,s_P}, Q_V) <= sum_i sum_j w_i w_j W(D_{V,s_P}^{(i)}, D_V^{(j)}). With uniform weights, the surrogate becomes (1/k) sum_i W(D_{V,s_P}^{(i)}, D_V^{(i)}).

6. WR-CP training objective (Eq. 19): min_theta sum_i E_{D_{XY}^{(i)}}[l(h_theta(x), y)] + beta sum_i W(D_{V,s_P}^{(i)}, D_V^{(i)}). Algorithm 1 implements training plus an inference phase that runs Tibshirani et al. 2019-style importance-weighted conformal prediction at test time, with weights dQ_X(x)/dP_X(x). Densities are obtained by Gaussian KDE with bandwidth grid search; conformal score distributions are estimated point-wise to keep gradients differentiable. When beta = 0, WR-CP reduces to IW-CP.

## Datasets / experiments
Six regression datasets, each split into k source distributions:

- Airfoil self-noise (Brooks and Marcolini 2014), k = 3.
- Seattle-loop (Cui et al. 2019), traffic speed, k = 10.
- PeMSD4, PeMSD8 (Guo et al. 2019), traffic speed, k = 10.
- Japan-Prefectures, U.S.-States (Deng et al. 2020), epidemic spread forecasting, k = 10.

Architecture: MLP (input, 64, 64, 1). 10 sampling trials per dataset; 10k generated test sets per trial. Calibration sampled from union of k subsets.

Headline numerical results: average coverage gap reduced to 3.2% across alpha in [0.1, 0.9]; prediction sets 37% smaller than worst-case CP (WC-CP). Spearman correlation between W-distance and coverage gap is highest among (W, TV, KL, ExpDiff) on 5 of 6 datasets (Table 1: W gets 0.59, 0.84, 0.90, 0.84, 0.77, 0.57). Pareto fronts (Figure 5) show beta sweeps trade off coverage gap vs set size; beta = 0 collapses to IW-CP.

Baselines: vanilla CP, IW-CP (Tibshirani et al. 2019), CQR (Romano et al. 2019), WC-CP (Gendler et al. 2021; Cauchois et al. 2024; Zou and Liu 2024).

## Metrics
Coverage gap |empirical coverage - (1 - alpha)| swept over alpha in [0.1, 0.9]; prediction set size; normalized 1-Wasserstein distance between calibration and test score distributions; Spearman rank correlation between distance metrics and coverage gap.

## What we steal for ConfPert
ConfPert evaluates on Norman 2019, Replogle K562, Replogle RPE1, Adamson, Tahoe-100M, and VCC, with a held-out cell-line split that is exactly a joint distribution shift: covariate shift across cell lines / batches plus concept shift in the perturbation-response map f. WR-CP is the right formal foundation here because:

- Eq. 13 cleanly attributes coverage gap to (a) covariate Wasserstein between cell-line embeddings and (b) concept Wasserstein on the response function. ConfPert can report both terms as a diagnostic for any held-out-cell-line evaluation rather than reporting only marginal coverage.
- The IW + W-regularizer training recipe transfers: weight calibration single cells by KDE likelihood ratio between source and target cell-line embeddings (e.g., scVI latents), then add the W-regularizer on conformity scores per source domain to control the concept-shift term.
- The multi-source convex-hull guarantee (Section 5, Cauchois et al. 2024 result) gives ConfPert a coverage statement when the test cell line is a mixture of training cell lines, which is realistic for batch-level shift even within a fixed cell line.
- Theorem 3 supplies a finite-sample empirical bound with explicit dependence on n, m, and the upper Wasserstein dimension, which we need because perturbation calibration sets are often small (a few hundred guide-cell pairs per perturbation).

Cite as the formal basis for ConfPert split type 3 (held-out cell line / batch) coverage claims, and as the regularizer template for the perturbation-prediction backbone.

## What we wrap
A procedure (training-time regularizer plus inference-time importance-weighted CP), not a model. Backbone hypothesis h_theta is arbitrary; ConfPert can keep its perturbation-response model and just add the W-regularization term on per-source score distributions during fine-tuning, then run IW-CP at test time.

## Failure modes / caveats
- Section H.1: KDE-based likelihood ratio is the bottleneck. WR-CP under-performs on airfoil (highest feature dimension d = 5, smallest |S^P_{XY}| = 500). For ConfPert, scVI / single-cell embeddings can be hundreds of dimensions, which makes raw-space KDE infeasible. Need a low-dim representation (e.g., 10-30 latent dims) before density estimation.
- Section H.2: bandwidth grid search is expensive in high dimensions; FastKDE (O'Brien et al. 2016) or normalizing-flow density ratios are recommended.
- Section H.3: results assume calibration distribution is a uniform mixture of training distributions; non-uniform calibration is untested. ConfPert calibration drawn from a single source cell line is exactly this untested regime.
- Worst-case CP still wins on guaranteed coverage; WR-CP gives a tunable accuracy-efficiency trade-off but loses formal worst-case coverage guarantee when beta is large enough to bias h_theta.
- Empirical bound Eq. 16 has rates n^{-1/sigma}, m^{-1/sigma} where sigma > d_W (upper Wasserstein dimension). In high-dim covariate spaces this rate degrades sharply, so the bound is meaningful only after dimension reduction.
- Importance weighting can enlarge Wasserstein distance on some datasets (Seattle-loop, Figure 3), so IW-CP alone is not monotonically beneficial.

## Code URLs
GitHub: https://github.com/rxu0112/WR-CP (released by authors).

## Verbatim quotes worth keeping
"Conformal prediction yields a prediction set with guaranteed 1 - alpha coverage of the true target under the i.i.d. assumption, which may not hold and lead to a gap between 1 - alpha and the actual coverage."

"Prior studies bound the gap using total variation distance, which cannot identify the gap changes under distribution shift at a given alpha."

"We first propose a Wasserstein distance-based upper bound of the coverage gap and analyze the bound using probability measure pushforwards between the shifted joint data and conformal score distributions, enabling a separation of the effect of covariate and concept shifts over the coverage gap."

"Equation (13) highlights how covariate and concept shifts impact the coverage gap. While the values of W(P_X, Q_X) and W(Q_{Y,f_P}, Q_Y) are inherent properties of given data and cannot be altered, the parameters kappa and eta are linked to the model h, allowing minimizing kappa and eta via optimizing h."

"Experiments on six datasets prove that WR-CP can reduce coverage gaps to 3.2% across different confidence levels and outputs prediction sets 37% smaller than the worst-case approach on average."
