# Chernozhukov, Wuthrich, Yinchu 2021 - Distributional Conformal Prediction

**Venue:** PNAS 2021 118(48). DOI 10.1073/pnas.2107794118.

## Claim
The paper proves finite-sample unconditional and asymptotic conditional coverage for prediction intervals built by conformalizing the conditional rank (probability integral transform) of any conditional CDF estimator, and constructs an "optimal" shape-adjusted variant that asymptotically attains the shortest possible (1 - alpha) prediction interval.

## Method
**Setup.** Data {(Y_t, X_t)}, t = 1,...,T, Y_t continuous scalar, X_t a p-vector of predictors. Goal: predict Y_{T+1} given X_{T+1}. Let F(y, x) = P(Y_t <= y | X_t = x). The probability integral transform implies U_t := F(Y_t, X_t) ~ Uniform(0,1) and U_t is independent of X_t.

**Baseline conformity score.** psi(y, x) = |F(y, x) - 1/2| applied to a plug-in estimator F-hat. Equivalently, the rank U-hat_t = F-hat(Y_t, X_t) is permuted; ranks live on (0,1) and have constant scale, unlike residuals.

**Algorithm 1 (Full DCP).** For each test value y in a grid Y_trial: (1) augment data Z^(y) = {(Y_t, X_t)}_{t=1..T} union (y, X_{T+1}); (2) refit F-hat^(y) on the augmented sample; (3) form ranks U-hat^(y)_t = F-hat^(y)(Y_t, X_t) for t <= T and U-hat^(y)_{T+1} = F-hat^(y)(y, X_{T+1}); (4) compute V-hat^(y)_t = psi(U-hat^(y)_t); (5) p-value p-hat(y) = (T+1)^{-1} sum_{t=1..T+1} 1{V-hat^(y)_t >= V-hat^(y)_{T+1}}; (6) return {y in Y_trial : p-hat(y) > alpha}.

**Algorithm 2 (Split DCP).** Split {1,...,T} into T1 (fit) and T2 (calibration). (1) fit F-hat on T1; (2) compute V-hat_t = psi(F-hat(Y_t, X_t)) for t in T2; (3) Q-hat_{T2} := the (1 - alpha)(1 + 1/|T2|) empirical quantile of {V-hat_t}; (4) return {y : psi(F-hat(y, X_{T+1})) <= Q-hat_{T2}}. Because F-hat(., X_{T+1}) is monotone, the output is an interval.

**Optimal DCP (Algorithm S1).** Replace psi with psi_*(y, x) = |F(y, x) - b(x, alpha) - (1 - alpha)/2|, where b(x, alpha) in arg min_{z in [0, alpha]} Q(z + 1 - alpha, x) - Q(z, x). Lemma 2 shows this score yields the conformal interval [Q(b(x, alpha), x), Q(b(x, alpha) + 1 - alpha, x)], which equals the population shortest (1 - alpha) interval. For symmetric unimodal F, b(x, alpha) = alpha/2 and optimal DCP coincides with baseline DCP.

**Theorems.**
- Theorem 1 (finite-sample unconditional validity): With exchangeable data and a permutation-invariant CDF estimator, P(Y_{T+1} in C-hat^full_{1-alpha}(X_{T+1})) >= 1 - alpha exactly.
- Theorem 2 (asymptotic unconditional validity, possibly misspecified F-hat -> F^* with F^* not necessarily F): coverage = 1 - alpha + o(1).
- Theorem 3 (asymptotic conditional validity under correct specification F^* = F): P(Y_{T+1} in C-hat^full_{1-alpha}(X_{T+1}) | X_{T+1}) = 1 - alpha + o_P(1).
- Theorem 4 (optimal DCP): asymptotic conditional validity plus mu(C-hat^conf_{1-alpha}(X_{T+1})) <= mu(C^opt_{1-alpha}(X_{T+1})) + o_P(1) (length convergence to optimum).
- Theorem 5 (stronger): symmetric difference mu(C-hat^conf_{1-alpha} Delta C^opt_{1-alpha}) = o_P(1) (set convergence), under additional consistency of b-hat.

**Assumption 1.** (i) consistency of the rank estimator under a transform phi: (T+1)^{-1} sum phi(|V-hat_t - V_t|) = o_P(1); (ii) Glivenko-Cantelli for the empirical rank CDF; (iii) Lipschitz/continuity of G. No rate is required.

## Datasets / experiments
**Stock returns (Section 5.1).** CRSP value-weighted daily portfolio, 1926-07-01 to 2021-06-30. Y_t = daily return; X_t = realized volatility = sqrt(sum of squared returns over last 22 days). Time-series setup with 5 rolling exercises (50% fit, 10% test, slid forward by 10%). Equal split |T1| = |T2|. 90% nominal. Twenty bins of X_t.
- Finding: mean-based CP (CP-OLS) drops to ~50% conditional coverage in high-volatility bins. DCP-QR and DCP-QR* hold close to 90% across all bins. DCP-DR undercovers in high-volatility regimes.

**CPS wages (Section 5.2).** 2012 CPS data (R hdm package). N = 29217. Y_i = hourly wage (skewed). dim(X_i) = 100 after two-way interactions (gender, marital, education, region, experience, experience^2). 80/20 train/test split, equal split inside, 20 random repetitions. 90% nominal.
- Unconditional coverage Table 1(a): all eight methods exactly 0.90.
- Dispersion of conditional coverage (x100) Table 1(b): DCP-QR 1.80, DCP-QR* 1.71, DCP-DR 3.08, CQR 2.21, CQR-m 2.36, CQR-r 2.30, CP-OLS 11.13, CP-loc 4.11.
- Average interval length Table 2: DCP-QR 34.22, DCP-QR* 29.61, DCP-DR 33.69, CQR 34.52, CQR-m 34.84, CQR-r 34.63, CP-OLS 33.84, CP-loc 32.66. Optimal DCP is shortest of the well-calibrated methods; CP-OLS and CP-loc are shorter only because they undercover.

**Motivating analytic example.** Y_t = X_t + X_t * eps_t with X_t ~ U(0,1), eps_t ~ N(0,1). Mean-based CP gives a fixed-length interval not adaptive to heteroskedasticity; DCP gives length 2x * Q_{|eps|}(1 - alpha), correctly scaling with x.

## Metrics
Marginal coverage rate, conditional coverage rate (binned by X or via logistic regression of indicator on X), conditional and average interval length. Dispersion of predicted coverage probabilities around 1 - alpha as a single-number conditional-coverage summary.

## What we steal for ConfPert
1. **Split conformal scaffold (Algorithm 2)** as the calibration template: split fit/calibration, compute scores on calibration set, take the (1 - alpha)(1 + 1/|T2|) empirical quantile, threshold scores at test time. ConfPert wraps this around six distributional discrepancies instead of psi(F(y, x)).
2. **Theorem 1** (finite-sample marginal coverage under exchangeability) is the canonical citation for K1 (general framework). Cite for the marginal coverage guarantee under any score, including our distributional scores, as long as the scoring function is permutation-invariant in the calibration set.
3. **Probability-integral-transform construction** as the design pattern: build a score whose distribution under exchangeability is pivotal (here uniform on (0,1)). For ConfPert, the analog is sample-based discrepancies whose null distribution under matched perturbation is data-driven but exchangeable across calibration units.
4. **Optimal-shape adjustment (Lemma 2 / Algorithm S1)** as a roadmap for length-optimal conformal sets when the score is asymmetric. We will cite this when arguing that ConfPert can be tuned to shortest-length distributional balls under skewed discrepancy distributions.
5. **Assumption 1** as the regularity template: o_P(1) consistency in a weak norm with no rate, plus Glivenko-Cantelli on the score's empirical CDF. ConfPert's asymptotic guarantees should imitate this minimal-rate structure.

## What we wrap
DCP's split-conformal calibration loop (Algorithm 2) is wrapped as the outer loop of ConfPert. The inner score psi(F(y, x)) = |F(y, x) - 1/2| is replaced by the six distributional discrepancies (KS, W1, energy, MMD-RBF, bimodality coefficient match, variance ratio) computed between predicted and held-out cell-population samples. Quantile-regression and distribution-regression sample-producers from Section B (or, in our case, CPA / scGen / CellOT / GEARS / scGPT / biolord / sVAE+ / STATE / scDFM / CellFlow) substitute for F-hat. The optimal-DCP variant is wrapped only when a discrepancy's null distribution is asymmetric.

## Failure modes / caveats
- Conditional validity (Theorem 3) needs F^* = F; under misspecification only marginal coverage survives. Foygel Barber et al. 2021 prove that distribution-free conditional validity is impossible for non-trivial sets, so DCP's "conditional" guarantee is conditional on a learnable distribution class.
- Y_t is assumed continuous and scalar. The paper does not develop multivariate Y. Extension to multivariate cell-population summaries is non-trivial: F(y, x) is replaced by some scalar score and the PIT no longer applies directly.
- Full DCP requires refitting F-hat on every grid point y; computationally heavy unless split DCP is used.
- Time-series extension needs beta-mixing for the rank consistency. For our perturbation setting we should justify exchangeability across cells/perturbations explicitly.
- Theorem 4's length-optimality has only an upper-bound (mu <= mu_opt + o_P(1)); Theorem 5's set-convergence needs strong consistency of b-hat plus density-bounded-below conditions.
- Y_trial grid choice: paper recommends [-max|Y_t|, max|Y_t|], justified by P(|Y_{T+1}| >= max|Y_t|) <= 1/(T+1).

## Code URLs
- Replication code (R): https://github.com/kwuthrich/Replication_DCP
- arXiv preprint (v3, 2021-08-21): https://arxiv.org/abs/1909.07889

## Verbatim quotes worth keeping
1. "We propose a robust method for constructing conditionally valid prediction intervals based on models for conditional distributions such as quantile and distribution regression."
2. "Unlike regression residuals, ranks are independent of the predictors, allowing us to construct conditionally valid prediction intervals under heteroskedasticity."
3. "object conditional validity in the sense of (2) cannot be achieved in a distribution-free way for non-trivial predictions. By Vovk (2012); Lei and Wasserman (2014); Foygel Barber et al. (2021), any prediction set satisfying (2) for every probability distribution of (X_t, Y_t) has infinite Lebesgue measure with non-trivial probability."
4. "Constructing prediction intervals based on deviations from quantile estimates is similar to working with deviations from mean estimates, as the deviations are measured in absolute levels. By contrast, exploiting the probability integral transform, our approach is generic and relies on permuting ranks, which naturally have the same scaling on (0, 1)."
5. "Theorem 1 provides a model-free unconditional performance guarantee in finite samples, allowing for arbitrary misspecification of the model of the conditional CDF. On the other hand, it has a major theoretical drawback. Even with iid data, it provides no guarantee at all on conditional validity."
6. "DCP-QR* produces the shortest prediction intervals among of all methods. This demonstrates the practical advantage of the shape adjustment when the conditional distribution is skewed."
