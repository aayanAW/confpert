# Izbicki, Shimizu, Stern 2022 - CD-split and HPD-split

**Venue:** JMLR 2022 (v23, paper 20-797). arXiv:2007.12778v3 (Oct 6, 2021). Authors: Rafael Izbicki, Gilson Shimizu, Rafael B. Stern. Title: "CD-split and HPD-split: efficient conformal regions in high dimensions."

## Claim
Two split-conformal procedures that build prediction regions for regression which (i) keep finite-sample marginal coverage 1 - alpha, (ii) achieve local validity (CD-split, CD-split+) over a data-driven feature partition, (iii) under consistency of the conditional density estimator, asymptotically converge to the oracle highest predictive density (HPD) set and therefore satisfy asymptotic conditional coverage. The conformity score is built directly from an estimated conditional density f(y|x), so the resulting region can be a union of intervals when Y|X is multimodal. The paper introduces a third variant CD-split+ that replaces the original CD-split partition with a Voronoi partition on a "profile distance" and is the recommended practical method. HPD-split removes the partition tuning entirely by using an approximately X-independent score.

## Method

**Setup.** Split conformal. Train set D' fits density estimator f-hat = B(D'). Calibration set D yields split residuals U_i = g-hat(X_i, Y_i). Region: { y : g-hat(X_{n+1}, y) >= U_{floor(alpha)} }.

**CD-split residual (Def 11):** U_i := f-hat(Y_i|X_i). Region (Def 13, eq 1): C(x) = { y : f-hat(y|x) >= U_{floor(alpha)}(x) }, a level set of f-hat. The cutoff U_{floor(alpha)}(x) is computed only over calibration points in the same partition cell as x (Def 7), giving local validity.

**CD-split partition (Def 12):** x_a, x_b share a cell iff q-hat_alpha(x_a), q-hat_alpha(x_b) of f-hat(Y|X) lie in the same element of a partition I of R^+.

**CD-split+ partition (Def 16, 20).** Profile distance d^2(x_a,x_b) := integral ( H-hat(z|x_a) - H-hat(z|x_b) )^2 dz, where H(z|x) is the conditional cdf of Z = f(Y|X). Centroids by k-means++ on discretized cdf vectors; A is the Voronoi partition. Two points are close iff their split-residual conditional cdfs match, so partitions stay informative in high d.

**HPD-split residual (Def 14):** U_i := H-hat(f-hat(Y_i|X_i)|X_i) = integral_{y : f-hat(y|X_i) <= f-hat(Y_i|X_i)} f-hat(y|X_i) dy. Under a perfect estimator U_i ~ Unif(0,1) independent of X. Region (Def 15, eq 2): C(x) = { y : f-hat(y|x) >= H-hat^{-1}(U^{hpd}_{floor(alpha)}|x) }. No partition needed.

**Algorithms.** Alg 1 CD-split: split; fit f-hat; q-hat_alpha on D; build A(x_{n+1}); U = alpha-quantile of { f-hat(y_i|x_i) : (x_i,y_i) in A(x_{n+1}) }; return { y : f-hat(y|x_{n+1}) >= U }. Alg 2 HPD-split: split; fit f-hat; numerically integrate H-hat; U = alpha-quantile of { H-hat(f-hat(y_i|x_i)|x_i) }; return { y : H-hat(f-hat(y|x_{n+1})|x_{n+1}) >= U }. Alg 3 CD-split+: as Alg 1 but A is Voronoi via k-means++ on H-hat(.|x_i) with J = ceil(n/100).

**Theorems.**
- Theorem 6: convergence to the HPD set implies asymptotic conditional coverage.
- Theorem 8: under exchangeability of D, the generic split-with-partition method yields marginal validity and local validity w.r.t. partition A.
- Theorem 21: CD-split and CD-split+ satisfy local validity (and marginal validity) under exchangeability.
- Theorem 22: HPD-split satisfies marginal validity.
- Theorem 27: under Assumptions 23 - 26 (i.i.d.; consistency of f-hat with rates eta_n, rho_n = o(1); H(u|x) C^1 with bounded derivative and lower bound near alpha; Y bounded), HPD-split converges to the HPD set and satisfies asymptotic conditional validity.
- Theorem 28: CD-split converges to HPD and is asymptotically conditional if |I| = o(1) for every cell.
- Theorem 29: CD-split+ converges to HPD and is asymptotically conditional if |A| = o(1) for every cell of the Voronoi partition.

## Datasets / experiments
Synthetic regression at d = 20: Homoscedastic N(0.3 x_1, 1); Bimodal mixture from Lei and Wasserman 2014 with irrelevant features; Heteroscedastic N(0.3 x_1, 1+0.3|x_1|); Asymmetric Gamma. 5,000 replications, 1 - alpha = 0.90. Density estimator: FlexCode + random forest. Dimensionality sweep d in [0, 5000] at n = 1000. Classification: logistic-model simulation, then MNIST with conv-net density estimator. Real benchmark: photometric redshift on Happy A (74,950 SDSS DR12 galaxies, 64,950 train / 5,000 calib / 5,000 test, Gaussian mixture density net). Table 1 numbers at 1 - alpha = 0.80: bright coverage HPD 0.800, CD-split+ 0.795, Dist 0.802, Quantile 0.808, Reg 0.807; faint coverage HPD 0.788, CD-split+ 0.792, Dist 0.746, Quantile 0.754, Reg 0.658; faint avg size HPD 0.065, CD-split+ 0.066, Dist 0.064, Quantile 0.074, Reg 0.045. CD-split+ and HPD-split are the only methods reaching nominal coverage on multimodal faint galaxies. CD-split+ has the smallest conditional-coverage deviation in nearly every scenario; HPD-split is competitive but slightly worse. CD-split+ region size is flat in d up to d = 5000.

## Metrics
Marginal coverage P(Y in C(X)); conditional coverage absolute deviation E[ |P(Y in C(X)|X) - (1-alpha)| ]; average region size |C(X)|; size-stratified coverage violation (Angelopoulos 2020) for classification; CDE loss for the density estimator. Tuning sweep over partition size shows trade-off: larger partition reduces region size but conditional-coverage deviation has a U-shape (Figure 5).

## What we steal for ConfPert
ConfPert's "CD-split per-gene" head is Definition 13 applied per gene g: KDE-estimate marginal predictive density f-hat_g(y|x) from samples of the perturbation predictor, set U_i^g = f-hat_g(Y_{i,g}|X_i) on calibration cells, return per-gene level sets { y : f-hat_g(y|x) >= U_{floor(alpha)}^g }. Theorem 21 supplies the per-gene marginal-validity citation; Theorem 28 / Theorem 6 supply the asymptotic-conditional-coverage citation whenever the per-gene density estimator is consistent. HPD-split (Definition 15) is the right ConfPert variant for genes with bimodal predictive distributions (cell-fate genes): score is approximately X-independent, no partition needed, no tuning. CD-split+ (Definition 20, Voronoi over profile distance) is the reference if we condition on perturbation embeddings, since profile distance dodges the curse of dimensionality on perturbation space.

## What we wrap
Not a model. A wrapper procedure on top of any sample-producing predictor: KDE-estimate the per-gene predictive density from samples, then apply CD-split / CD-split+ / HPD-split per gene. The joint / multivariate-Y extension is not provided here (use DCP, OT-CP, MMD-conformal for the joint distributional discrepancy heads; Izbicki 2022 governs only the per-gene marginal head).

## Failure modes / caveats
- Conditional coverage is asymptotic only (Definition 4: sup |P(Y in C|X=x) - (1-alpha)| = o(1) on a high-probability set), not finite-sample. Marginal coverage is finite-sample.
- Requires Assumption 24 consistency of f-hat in sup-norm in expectation; if the per-gene KDE fails to converge (heavy tail, bandwidth misset), conditional coverage drops though marginal still holds.
- Assumption 25: H(u|x) needs strictly positive derivative near the alpha-quantile; a flat plateau in the predictive cdf around the level breaks the proof of Theorem 27.
- Joint application across all G genes is exponentially loose: applying CD-split jointly on Y in R^G suffers curse of dimensionality on the y grid; the per-gene marginal head sidesteps this. For joint coverage use DCP / OT-CP / MMD instead.
- KDE bandwidth selection. The paper does not prescribe a bandwidth method; for f-hat the authors use FlexCode (Izbicki and Lee 2017) plus random forest in simulations and a Gaussian mixture density network on Happy A. They argue empirically (Section 5.1, Figure 6) that "conditional density estimates with a smaller [CDE] loss lead to smaller prediction bands with a better conditional coverage" and recommend selecting f-hat by minimizing CDE loss. ConfPert should treat per-gene KDE bandwidth as a CDE-loss-minimization problem on a held-out fold.
- Region-of-arbitrary-shape interpretation: an interval is easier to communicate; using mixture-density estimators restricts the region to a union of a fixed number of intervals (Section 7).

## Code URLs
R package implementing CD-split, CD-split+, HPD-split: https://github.com/rizbicki/predictionBands . Conditional density backend FlexCode: https://github.com/rizbicki/FlexCode . Python conditional density tooling: Dalmasso et al. 2020 (cdetools).

## Verbatim quotes worth keeping
- "CD-split however contains many tuning parameters, and their role is not clear. In this paper, we provide new insights on CD-split by exploring its theoretical properties. In particular, we show that CD-split converges asymptotically to the oracle highest predictive density set and satisfies local and asymptotic conditional validity."
- "We also propose HPD-split, a variation of CD-split that requires less tuning, and show that it shares the same theoretical guarantees as CD-split."
- "Theorem 6. If a conformal method converges to the highest predictive density set, then it satisfies asymptotic conditional coverage."
- "Theorem 21. If the instances in D are exchangeable, then CD-split and CD-split+ satisfy local validity (Definition 3) with respect to the partition in Definition 12. In particular, CD-split and CD-split+ also satisfy marginal validity (Definition 1)."
- "Theorem 27. Under Assumptions 23 to 26, HPD-split converges to the hpd set and satisfies asymptotic conditional validity."
- "The CD-split residual is U_i := f-hat(Y_i | X_i)." (Definition 11)
- "The HPD-split residual is given by U_i := H-hat( f-hat(Y_i|X_i) | X_i ) = integral_{y : f-hat(y|X_i) <= f-hat(Y_i|X_i)} f-hat(y|X_i) dy." (Definition 14)
- "R code for implementing CD-split, CD-split+ and HPD-split is available at: https://github.com/rizbicki/predictionBands ."
