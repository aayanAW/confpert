# Romano, Patterson, Candes 2019 - Conformalized Quantile Regression

**Venue:** NeurIPS 2019. arXiv:1905.03222.

## Claim
CQR wraps any quantile regression algorithm in a split-conformal calibration step on a held-out set, producing prediction intervals that are adaptive to heteroscedasticity and that satisfy the marginal, finite-sample, distribution-free coverage guarantee P{Y_{n+1} in C(X_{n+1})} >= 1 - alpha under exchangeability.

## Method
Split the n training pairs (X_i, Y_i) into a proper training set indexed by I_1 and a calibration set indexed by I_2. On I_1, fit two conditional quantile estimators q_hat_alpha_lo and q_hat_alpha_hi at levels alpha_lo = alpha/2 and alpha_hi = 1 - alpha/2 using a quantile regression algorithm A (the paper uses quantile random forests [Meinshausen 2006] and quantile neural networks [Taylor 2000] minimizing the pinball / check loss rho_alpha(y, y_hat) = alpha(y - y_hat) if y - y_hat > 0 else (1 - alpha)(y_hat - y)).

CQR conformity score (eq. 9, signed, accounts for both under- and over-coverage):
E_i = max{q_hat_alpha_lo(X_i) - Y_i, Y_i - q_hat_alpha_hi(X_i)}, i in I_2.

Final interval (eq. 10):
C(X_{n+1}) = [q_hat_alpha_lo(X_{n+1}) - Q_{1-alpha}(E, I_2), q_hat_alpha_hi(X_{n+1}) + Q_{1-alpha}(E, I_2)],

where the inflated empirical quantile (eq. 11) is
Q_{1-alpha}(E, I_2) := the (1 - alpha)(1 + 1/|I_2|)-th empirical quantile of {E_i : i in I_2}.
This is the finite-sample-corrected index k = ceil((1-alpha)(|I_2|+1)) in the order statistics of E_i.

Algorithm 1 split steps:
1. Randomly split {1, ..., n} into disjoint I_1, I_2.
2. {q_hat_alpha_lo, q_hat_alpha_hi} <- A({(X_i, Y_i) : i in I_1}).
3. Compute E_i for i in I_2 via eq. 9.
4. Compute Q_{1-alpha}(E, I_2).
5. Output C(x) = [q_hat_alpha_lo(x) - Q_{1-alpha}(E, I_2), q_hat_alpha_hi(x) + Q_{1-alpha}(E, I_2)].

Theorem 1 (marginal coverage): If (X_i, Y_i) for i = 1,...,n+1 are exchangeable, then P{Y_{n+1} in C(X_{n+1})} >= 1 - alpha. If the E_i are almost surely distinct, the upper bound P{Y_{n+1} in C(X_{n+1})} <= 1 - alpha + 1/(|I_2| + 1) also holds, so the procedure is nearly perfectly calibrated. Proof reduces, conditionally on I_1, to the inflated-empirical-quantile lemma (Lemma 2, Appendix A) applied to the exchangeable conformity scores E_i, i in I_2 union {n+1}.

Theorem 2 (asymmetric two-tailed variant, controls each tail separately): With E_i^lo := q_hat_alpha_lo(X_i) - Y_i and E_i^hi := Y_i - q_hat_alpha_hi(X_i), define
C(X_{n+1}) := [q_hat_alpha_lo(X_{n+1}) - Q_{1-alpha_lo}(E^lo, I_2), q_hat_alpha_hi(X_{n+1}) + Q_{1-alpha_hi}(E^hi, I_2)].
Then under exchangeability, P{Y_{n+1} >= q_hat_alpha_lo(X_{n+1}) - Q_{1-alpha_lo}(E^lo, I_2)} >= 1 - alpha_lo and the symmetric upper-tail statement, so for alpha = alpha_lo + alpha_hi, P{Y_{n+1} in C(X_{n+1})} >= 1 - alpha. Cost: slightly wider intervals.

Base learners: CQR-RF uses quantile regression forests (the lower / upper nominal quantiles are tunable hyper-parameters, cross-validated to minimize average interval length). CQR-NN uses a 3-layer fully connected ReLU network (64-64), Adam (lr 5e-4, batch 64, weight decay 1e-6), dropout 0.1, early stopping. To save compute, a single 2-output network jointly predicts (q_lo, q_hi) trained with the pinball loss summed over both quantiles. Tuning the nominal quantiles below the target alpha_lo and alpha_hi does NOT invalidate Theorem 1 because the calibration step is what enforces coverage.

## Datasets / experiments
11 benchmarks: MEPS_19, MEPS_20, MEPS_21 (medical expenditure panel survey), blog_data (BlogFeedback), bio (protein tertiary structure), bike (bike sharing), community (communities and crimes), STAR (Tennessee student-teacher achievement ratio), concrete (compressive strength), facebook_1, facebook_2 (Facebook comment volume). Features standardized; response rescaled by mean absolute value. 20 random 80/20 train/test splits per dataset; |I_1| = |I_2|; alpha = 0.1; 2200 total experiments.

Table 1 averages across all 11 datasets and 20 splits (length / coverage):
- Ridge: 3.06 / 90.03
- Ridge Local: 2.94 / 90.13
- Random Forests: 2.24 / 89.99
- Random Forests Local: 1.82 / 89.95
- Neural Net: 2.16 / 89.92
- Neural Net Local: 1.81 / 89.95
- CQR Random Forests: 1.41 / 90.33
- CQR Neural Net: 1.40 / 90.05
- *Quantile Random Forests (no conformalization): 2.23 / 92.62 (overcovers)
- *Quantile Neural Net (no conformalization): 1.49 / 88.51 (undercovers)

CQR wins on 10/11 datasets vs both standard and locally adaptive split conformal. CQR-RF is overly conservative on the two Facebook datasets due to ties among conformity scores (Theorem 1 upper bound does not apply with ties). Asymmetric Theorem 2 increases CQR-NN length 1.40 -> 1.58 and CQR-RF length 1.41 -> 1.57 (CQR-RF coverage rises 90.33 -> 90.99) for the per-tail guarantee.

## Metrics
Marginal coverage at alpha = 0.1; average prediction interval length; per-dataset breakdowns in Figures 3-6 with heteroscedasticity-stratified visual inspection of the locally varying interval widths as the implicit conditional-coverage diagnostic.

## What we steal for ConfPert
The split-conformal protocol with the finite-sample-corrected quantile index Q_{1-alpha} = ceil((1-alpha)(|I_2|+1))-th order statistic. The CQR conformity score (eq. 9) E_i = max{q_lo(X_i) - Y_i, Y_i - q_hi(X_i)} is the per-gene quantile-band correction inside ConfPert's per-gene-marginal CD-split head (already implemented in HetPert). We cite Theorem 1 for marginal coverage and Theorem 2 for the per-tail variant when ConfPert needs asymmetric per-gene bands (e.g., when up- and down-regulation directions of a gene's perturbation response have asymmetric noise). CQR is UNIVARIATE Y in R; the multi-output / multivariate / distributional extension is acknowledged as future work in the paper itself ("We expect the ideas behind conformal quantile regression to be applicable in the related setting of conformal predictive distributions" - Conclusion). This is the explicit gap that justifies why ConfPert needs DCP (Chernozhukov 2021) for predictive distributions and OT-CP / multivariate-Gaussian / energy-distance / MMD scoring for joint multivariate coverage over the gene vector. CQR alone gives only per-gene-marginal coverage, not joint coverage on the full ~5000-dim expression vector.

## What we wrap
CQR does not produce samples - it produces a single interval [L(x), U(x)] per query. Not a wrappable predictor in the sample-based sense. But we lift the quantile-conformity score idea into ConfPert as the per-gene scoring rule inside our CD-split head, and we reuse the inflated-empirical-quantile finite-sample correction across every conformal head we ship.

## Failure modes / caveats
- Univariate Y. No multivariate or full-distribution guarantee.
- Quantile crossing (q_hat_lo > q_hat_hi at some x) for neural-net base learners; the paper reports post-hoc rearrangement [Chernozhukov, Fernandez-Val, Galichon 2010] reduces CQR-NN length 1.40 -> 1.35 with coverage essentially unchanged. Quantile RF cannot cross by construction.
- Sample splitting: equal-sized I_1 and I_2 in the experiments; data efficiency is halved relative to a no-split estimator. Cross-conformal or jackknife+ variants exist but are not the default.
- Marginal coverage only. Conditional coverage P{Y in C(X) | X = x} >= 1 - alpha is NOT guaranteed; only assessed indirectly through heteroscedasticity-stratified plots.
- Ties in conformity scores break the upper bound 1 - alpha + 1/(|I_2|+1) (Facebook datasets show resulting overcoverage).
- Validity assumes exchangeability; covariate shift requires the weighted variant of [Tibshirani et al. 2019, ref 50].

## Code URLs
https://github.com/yromano/cqr (footnote 1, p. 2: "Source code implementing CQR is available online at https://github.com/yromano/cqr.")

## Verbatim quotes worth keeping
1. (Conclusion, p. 11, the seed of DCP) "We expect the ideas behind conformal quantile regression to be applicable in the related setting of conformal predictive distributions [49]. In this extension of conformal prediction, the aim is to estimate a predictive probability distribution, not just an interval. We see intriguing connections between our work and a very recent, independently written paper on conformal distributions [17]."
2. (Theorem 1, p. 5) "If (X_i, Y_i), i = 1, . . . , n + 1 are exchangeable, then the prediction interval C(X_{n+1}) constructed by the split CQR algorithm satisfies P{Y_{n+1} in C(X_{n+1})} >= 1 - alpha. Moreover, if the conformity scores E_i are almost surely distinct, then the prediction interval is nearly perfectly calibrated: P{Y_{n+1} in C(X_{n+1})} <= 1 - alpha + 1/(|I_2| + 1)."
3. (Section 4, p. 4, eq. 9 motivation) "The conformity score thus accounts for both undercoverage and overcoverage."
4. (Section 6.2, p. 9, on validity) "Turning to the issue of valid coverage, all methods based on conformal prediction successfully construct prediction bands at the nominal coverage rate of 90%, as the theory suggests they should."
5. (Introduction, p. 1, the design goal we inherit) "When the data is heteroscedastic, getting valid but short prediction intervals requires adjusting the lengths of the intervals according to the local variability at each query point in predictor space."
