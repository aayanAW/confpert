# Tibshirani, Foygel Barber, Candes, Ramdas 2019 - Weighted Conformal Prediction under Covariate Shift

**Venue:** NeurIPS 2019. arXiv:1904.06019 (v3, 6 Jul 2020).

## Claim

Weighted-quantile conformal prediction maintains finite-sample marginal coverage at level 1 - alpha when the test covariate distribution P_X-tilde differs from the training P_X, provided the likelihood ratio w(x) = dP_X-tilde / dP_X is known (or accurately estimable from unlabeled test covariates), with conditional Y|X unchanged.

## Method

**Setup (Eq. 6).** Training (X_i, Y_i) iid from P = P_X x P_{Y|X} for i = 1,...,n; test point (X_{n+1}, Y_{n+1}) ~ P-tilde = P_X-tilde x P_{Y|X}, independently. Same conditional, shifted covariate.

**Weighted quantile (Eq. 7).** Replace the empirical distribution of nonconformity scores with a weighted version sum_i p_i^w(x) delta_{V_i^{(x,y)}} + p_{n+1}^w(x) delta_infinity, with weights

    p_i^w(x) = w(X_i) / [sum_{j=1}^n w(X_j) + w(x)],   i = 1,...,n
    p_{n+1}^w(x) = w(x) / [sum_{j=1}^n w(X_j) + w(x)].

**Coverage theorem (Corollary 1).** With P_X-tilde absolutely continuous w.r.t. P_X and w = dP_X-tilde/dP_X, define

    C-hat_n(x) = { y in R : V_{n+1}^{(x,y)} <= Quantile(1 - alpha; sum_i p_i^w(x) delta_{V_i^{(x,y)}} + p_{n+1}^w(x) delta_infty) }.

Then P{Y_{n+1} in C-hat_n(X_{n+1})} >= 1 - alpha. Remark 3 notes the same holds with w only known up to a normalizing constant (the constant cancels in Eq. 7).

**Weighted exchangeability (Definition 1).** V_1,...,V_n are weighted exchangeable with weights w_1,...,w_n if their joint density factors as f(v_1,...,v_n) = prod_i w_i(v_i) * g(v_1,...,v_n) for some symmetric g. Lemma 2: independent draws Z_i ~ P_i with each P_i absolutely continuous w.r.t. P_1 are weighted exchangeable with w_1 = 1, w_i = dP_i/dP_1, which encompasses the covariate-shift case.

**General weighted quantile lemma (Lemma 3) and Theorem 2.** For weighted exchangeable Z_1,...,Z_{n+1} with weights w_1,...,w_{n+1},

    p_i^w(z_1,...,z_{n+1}) = sum_{sigma:sigma(n+1)=i} prod_j w_j(z_{sigma(j)}) / sum_sigma prod_j w_j(z_{sigma(j)}),

and the corresponding weighted-quantile conformal band has marginal coverage >= 1 - alpha (Theorem 2). Corollary 1 follows because w_i = w(X_i) (or w_i = 1 for training and w_{n+1} = w on the augmented sample) reduces the combinatorial sum to Eq. 7.

**Split version (Eq. 10).** Pre-fit mu_0 on a held-out set; for x in R^d,

    C-hat_n(x) = mu_0(x) +/- Quantile(1 - alpha; sum_i p_i^w(x) delta_{|Y_i - mu_0(X_i)|} + p_{n+1}^w(x) delta_infty).

Coverage holds conditional on the pre-fit set.

**Estimated weights (Eq. 12).** Fit a probabilistic classifier p-hat(x) on (X_i, C_i) with C_i = 0 for training points and C_i = 1 for unlabeled test points; use w-hat(x) = p-hat(x) / (1 - p-hat(x)). Density-ratio estimation by classification, conditional-odds form. The paper uses logistic regression and random forests (with p-hat clipped to [0.01, 0.99] to prevent infinite weights, footnote 3).

## Datasets / experiments

**Airfoil (UCI).** N = 1503 NASA airfoil self-noise observations, response Y = scaled sound pressure level, d = 5 covariates (log frequency, angle of attack, chord length, free-stream velocity, suction side log displacement thickness). 5000 random partitions into D_pre (25%, fit mu_0 by linear regression), D_train (25%, calibration), D_test (50%, exchangeable test), and D_shift (25% of D_test resampled with replacement with probabilities proportional to w(x) = exp(x^T beta), beta = (-1, 0, 0, 0, 1); exponential tilting on covariates 1 and 5).

**Headline results (Figs. 2-3, alpha = 0.1, target 90%).**

- No shift, ordinary split conformal: avg coverage 90.2%.
- Shift D_shift, ordinary (unweighted) split conformal: avg coverage 82.2% (substantial undercoverage).
- Shift, oracle weights w: avg coverage 90.8% (recovered).
- Shift, logistic-regression w-hat: avg coverage ~91.0%.
- Shift, random-forest w-hat: avg coverage ~91.0%.
- Median interval lengths: oracle-weighted intervals are longer than unweighted ones at matched effective sample size, attributable to the smaller effective sample size n-hat = ||w(X_{1:n})||_1^2 / ||w(X_{1:n})||_2^2 (Gretton et al. 2009; Reddi et al. 2015).
- Sanity check (Fig. 4, no shift): weighted methods are nearly indistinguishable from unweighted, so the procedure does not over-inflate when shift is absent.

## Metrics

Marginal coverage P{Y_{n+1} in C-hat_n(X_{n+1})} averaged over 5000 random splits, plus histogram of empirical coverages, and median interval length.

## What we steal for ConfPert

ConfPert needs validity for cross-cell-line generalization (e.g., calibrate on Norman K562, deploy on RPE1) where cell-state covariate distribution shifts but the response model Y|X (perturbation effect given expression) is approximately preserved. This paper is the textbook procedure: fit a domain classifier on K562 vs RPE1 covariates to obtain w-hat by Eq. 12, then plug w-hat into the weighted split-conformal interval (Eq. 10). Cite as the foundation for ConfPert's cross-distribution coverage claim. Pair with arXiv:2501.13430 (Wasserstein-regularized CP) which complements this when the likelihood ratio is ill-defined (overlap failure) by shifting from importance reweighting to Wasserstein robustness.

The graphical-models extension (Section 4) is directly relevant: when X (cell expression) is high-dimensional but Z (cell-line label, low-dimensional) drives the shift, estimate the ratio on Z alone via P-tilde_Z / P_Z (a marginal probability ratio), avoiding high-dimensional density estimation.

## What we wrap

A procedure, not a model. ConfPert wraps any base regressor (residual networks, FLM, baseline OLS) and produces calibrated intervals under cross-distribution evaluation. Complement to RCPS (risk-controlling sets) and OT-CP (optimal transport CP).

## Failure modes / caveats

- **Estimated w-hat:** Theoretical 1 - alpha coverage holds only with the true w. With w-hat, coverage is approximate; the paper demonstrates empirical near-validity but proves no finite-sample bound for the estimated case.
- **Effective sample size:** Heavy weighting collapses the calibration set to n-hat = ||w(X_{1:n})||_1^2 / ||w(X_{1:n})||_2^2, inflating intervals. For ConfPert this means severe K562 -> RPE1 mismatch will widen intervals proportional to the tail of w.
- **No nontrivial upper bound under weighting (Remark 5):** the largest jump in the weighted CDF can equal max_i p_i^w, which can be near 1 under heavy concentration, so unlike the unweighted case there is no 1 - alpha + 1/(n+1) upper-tail bound without further conditions.
- **Overlap requirement:** P_X-tilde must be absolutely continuous w.r.t. P_X. If the test distribution puts mass where the training distribution does not (new RPE1-only cell states), the method breaks. Hand off to OT-CP / Wasserstein-CP.
- **Infinite weights from classifier saturation:** if p-hat(x) -> 1, w-hat -> infinity. Clip p-hat to [0.01, 0.99] or use a calibrated classifier.
- **Assumes Y|X invariant:** if the perturbation-response mechanism itself differs across cell lines (regulatory rewiring), the i.i.d. Y|X assumption fails and weighted CP alone is insufficient.

## Code URLs

R code for the airfoil experiments: http://www.github.com/ryantibs/conformal/ (referenced in Section 2.3).

## Verbatim quotes worth keeping

> "We extend conformal prediction methodology beyond the case of exchangeable data. In particular, we show that a weighted version of conformal prediction can be used to compute distribution-free prediction intervals for problems in which the test and training covariate distributions differ, but the likelihood ratio between these two distributions is known---or, in practice, can be estimated accurately with access to a large set of unlabeled data (test covariate points)." (Abstract.)

> "The key realization is the following: if we know the ratio of test to training covariate likelihoods, dP-tilde_X/dP_X, then we can still perform a modified version conformal inference, using a quantile of a suitably weighted empirical distribution of nonconformity scores." (Section 2.)

> "Random variables V_1,...,V_n are said to be weighted exchangeable, with weight functions w_1,...,w_n, if the density f of their joint distribution can be factorized as f(v_1,...,v_n) = prod_i w_i(v_i) * g(v_1,...,v_n), where g is any function that does not depend on the ordering of its inputs." (Definition 1.)

> "However, the histogram is more dispersed than it is when there is no covariate shift... This is because, by using a quantile of the weighted empirical distribution of nonconformity scores, we are relying on a reduced effective sample size." (Section 2.3, on the dispersion of weighted-coverage histograms.)

> "When the likelihood ratio dP-tilde_X/dP_X is not known, it can be estimated given access to unlabeled data (test covariate points), which we showed empirically, on a low-dimensional example, can still yield correct coverage." (Discussion.)

ConfPert context: cross-dataset and cross-cell-line evaluation in ConfPert requires coverage under distribution shift. This paper is the canonical NeurIPS reference. Pair with arXiv:2501.13430 (Wasserstein-regularized CP) for full coverage of the shift regime.
