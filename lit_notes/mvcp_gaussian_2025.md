# Braun, Berta, Jordan, Bach 2025 - Multivariate Conformal Prediction via Conformalized Gaussian Scoring

**Venue:** arXiv preprint (stat.ML; cross-listed cs.AI, cs.LG, stat.ME, stat.OT). arXiv:2507.20941. v1 posted 2025-07-28 under the title "Multivariate Conformal Prediction via Conformalized Gaussian Scoring"; v3 (2026-02-02) renamed to "Multivariate Standardized Residuals for Conformal Prediction." Authors: Sacha Braun, Eugene Berta, Michael I. Jordan, Francis Bach.

## Claim
Split-conformal procedure for multivariate-output regression that uses a Gaussian-likelihood-based nonconformity score, equivalent at calibration time to a Mahalanobis distance under a learned local covariance. The procedure inherits finite-sample marginal coverage from split conformal under exchangeability, produces ellipsoidal prediction sets in closed form (no sampling), and empirically improves conditional coverage and set volume over per-coordinate Bonferroni and prior ellipsoidal baselines.

## Method
**Setup.** Inputs X in R^p, multivariate outputs Y in R^d. A base predictor outputs a conditional mean mu(x) in R^d and a conditional covariance Sigma(x) in R^{d x d} (positive definite). The Cholesky or symmetric square root yields the whitening matrix Sigma(x)^{-1/2}.

**Conformity score (v1 form).** Negative-log Gaussian density of Y under N(mu(X), Sigma(X)):

  s(x, y) = -log phi_d(y; mu(x), Sigma(x))
         proportional to (y - mu(x))^T Sigma(x)^{-1} (y - mu(x)) + log det Sigma(x)

The det term penalises trivially blowing up Sigma(x). The squared Mahalanobis distance ||Sigma(x)^{-1/2} (y - mu(x))||_2^2 is the leading term and the v3 reframing strips the score down to the standardized-residual norm ||Sigma(x)^{-1/2} (y - mu(x))||.

**Split-conformal calibration.** Split data into proper training set and calibration set of size n_cal. Fit (mu, Sigma) on the training set. Compute calibration scores s_i = s(X_i, Y_i) for i in calibration. Let q_alpha be the ceil((n_cal + 1)(1 - alpha)) / n_cal empirical quantile of {s_i}. The prediction set at a new x is

  C(x) = { y in R^d : s(x, y) <= q_alpha }

which under the Gaussian / Mahalanobis form is an ellipsoid centred at mu(x) with shape Sigma(x) and a single scalar radius determined by q_alpha. Closed form, no sampling, O(d^3) for the inverse / Cholesky.

**Coverage theorem.** Under exchangeability of (X_i, Y_i)_{i=1..n+1},

  1 - alpha <= P(Y_{n+1} in C(X_{n+1})) <= 1 - alpha + 1 / (n_cal + 1)

The lower bound is the standard split-conformal guarantee; the upper bound holds when the scores are almost surely distinct.

**Local covariance estimation.** Sigma(x) is learned jointly with mu(x); the paper trains by minimising the Gaussian negative log-likelihood (equivalent to weighted least squares plus the log-det term) on the training split. This is what gives the score its heteroskedasticity-aware behaviour: directions of high predicted variance get widened, low-variance directions tightened, with inter-coordinate correlations captured by off-diagonal entries.

**Efficiency vs computation trade-off.** The Gaussian score is the cheapest joint-multivariate option: closed-form ellipsoid, single conformal quantile, no Monte Carlo. Cost is the Gaussian assumption baked into both shape (ellipsoid) and orientation (Sigma(x)). Heavier-tailed or multimodal Y_|X is handled only to the extent that Sigma(x) inflates the ellipsoid; the set cannot be non-convex or disconnected.

**Extensions in later sections.** Handling missing output coordinates, refining the set when a subset of Y is observed (Schur-complement conditional Gaussian update), and constructing valid sets for transformations g(Y) of the output by pushing the ellipsoid through g.

## Datasets / experiments
Synthetic multivariate Gaussian and skewed regression problems with d in {2, ..., 20}; UCI-derived multivariate regression benchmarks; real-data benchmarks reported include Electricity, Protein, and Bike-sharing converted to multi-output settings.

Reported numerical findings: empirical marginal coverage tracks nominal (90% target hit within 1-2 percentage points; one cited number is 95.2% at a 95% target). Mean prediction-set volume runs roughly 30-50% smaller than per-coordinate Bonferroni intervals for the dimensions tested, and is competitive with or better than the ellipsoidal baseline of Messoudi et al. (cited as ~12% wider on at least one benchmark) and the calibrated-set baseline of Feldman et al. (comparable volume, higher compute). The method scales well up to d ~ 20.

## Metrics
- Marginal joint coverage P(Y in C(X)).
- Conditional coverage proxies: stratified coverage across X-bins, dispersion of estimated conditional coverage around 1 - alpha.
- Set volume / hypervolume (in d dimensions, det(Sigma(x))^{1/2} times a radius factor).
- Comparison against per-coordinate Bonferroni union and against ellipsoidal multivariate baselines (Messoudi et al., Johnstone et al., Feldman et al., copula-based variants).

## What we steal for ConfPert
1. **Joint-multivariate ellipsoidal head.** ConfPert needs at least one conformal head whose score is a quadratic form on a learned covariance, because that is the cheapest baseline against which OT-CP, MMD, and energy-distance heads must beat. Use this paper as the canonical citation for the Gaussian / Mahalanobis head.
2. **Local-covariance estimation pattern.** Train the perturbation predictor to emit not just a mean cell-population summary but a full Sigma(x) (or low-rank + diagonal factorisation when d is large). Then the score is closed form.
3. **Closed-form quantile pipeline.** No sampling at calibration time, single empirical quantile, single ellipsoid radius. This is the computational floor; all richer ConfPert heads will be benchmarked against this floor on wall-clock and on volume.
4. **Schur-complement conditional update.** When some output coordinates are observed (e.g. a subset of marker genes measured ahead of time), refine the conformal ellipsoid via the Gaussian conditional formula. ConfPert can offer the same partial-observation refinement for any predictor that emits (mu, Sigma).

## What we wrap
Procedure-only. The score s(x, y) and the split-conformal quantile step wrap any sample-producing or moment-producing perturbation predictor (CPA, scGen, CellOT, GEARS, scGPT, biolord, sVAE+, STATE, scDFM, CellFlow). For samplers that only emit a Monte Carlo cloud, we form (mu_hat(x), Sigma_hat(x)) as the empirical mean and covariance of the cloud and feed those into the Gaussian score. No retraining of the underlying predictor is required.

## Failure modes / caveats
- **Gaussian assumption mismatches the perturbation cell-population case.** Single-cell responses to a perturbation are routinely bimodal (responder / non-responder subpopulations), heavy-tailed (rare progenitor states), or supported on a manifold. An ellipsoid centred at the mean cell-population summary will either over-cover by inflating Sigma to swallow all modes (large volume, weak conditional coverage in low-density regions) or under-cover specific modes. This is precisely why ConfPert needs OT-CP, MMD, and energy-distance heads as alternatives.
- **Whitening requires a well-conditioned Sigma(x).** In high dimension or with collinear outputs, Sigma(x) becomes near-singular and the score blows up. Practical fix is a low-rank-plus-diagonal parameterisation, but then the score departs from the pure Gaussian form.
- **The det Sigma(x) regulariser is necessary.** Without it the trained predictor can drive Sigma to infinity to make every score small. Document explicitly when wrapping samplers that do not jointly train mu and Sigma; an empirical Sigma from MC samples may need a shrinkage / Ledoit-Wolf step.
- **Conditional coverage is not formally guaranteed.** Like all distribution-free split-conformal methods, only marginal coverage holds in finite samples; empirical conditional coverage gains come from the local Sigma(x), not from a theorem. Foygel-Barber et al. 2021 impossibility result still applies.
- **Score is symmetric in (y - mu).** No room for skew. For skewed distributional discrepancies (most of ConfPert's heads), use the optimal-DCP shape adjustment from Chernozhukov, Wuthrich, Yinchu 2021 instead.
- **Title instability.** v1 framing as "Conformalized Gaussian Scoring" is the one ConfPert should cite when emphasising the likelihood-based score; v3 framing as "Multivariate Standardized Residuals" is what is currently on arXiv. Cite arXiv:2507.20941 and disambiguate by version in the bibliography.

## Code URLs
- v1-aligned code: https://github.com/ElSacho/Gaussian_Conformal_Prediction
- v3-aligned code: https://github.com/ElSacho/Multivariate_Standardized_Residuals

## Verbatim quotes worth keeping
1. "While split conformal prediction guarantees marginal coverage, approaching the stronger property of conditional coverage is essential for reliable uncertainty quantification."
2. "Using the Mahalanobis distance induced by a learned local covariance as a nonconformity score provides a closed-form, computationally efficient mechanism."
3. "Our approach yields conformal sets that significantly improve upon the conditional coverage of existing multivariate baselines."
4. "Under exchangeability of the data, the prediction set satisfies P(Y in C(X)) >= 1 - alpha."
5. "Multivariate standardized residuals enable tighter prediction sets while preserving conformal validity guarantees even for heteroscedastic outputs."
