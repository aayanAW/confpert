# Angelopoulos and Bates 2021 - Gentle Intro to Conformal Prediction

**Venue:** Foundations and Trends-style monograph. arXiv:2107.07511 (v6, December 8, 2022). Updated through 2024.

## Claim
This monograph establishes that conformal prediction is a model-agnostic, distribution-free procedure that converts any heuristic notion of uncertainty from a pre-trained predictor into prediction sets with finite-sample, non-asymptotic marginal coverage of at least 1 minus alpha, requiring only exchangeability of calibration and test data and a user-chosen score function.

## Method
The monograph centers on **split (inductive) conformal prediction** (Section 1) as the workhorse. Given a pre-trained model f-hat, calibration data of size n, a score function s(x,y) where larger means worse agreement, and a target miscoverage alpha, the procedure is: (1) compute scores s_i = s(X_i, Y_i) on calibration data; (2) take q-hat as the ceiling((n+1)(1-alpha))/n empirical quantile of those scores; (3) form C(X_test) = {y : s(X_test, y) <= q-hat}.

**Theorem 1 (Conformal coverage guarantee, Vovk-Gammerman-Saunders):** if (X_i, Y_i)_{i=1..n} and (X_test, Y_test) are i.i.d. (or exchangeable), then P(Y_test in C(X_test)) >= 1 - alpha. The two-sided version (Eq. 1) gives 1 - alpha <= P(Y_test in C(X_test)) <= 1 - alpha + 1/(n+1). Conditional on the calibration set, the coverage random variable is Beta(n+1-l, l) with l = floor((n+1)*alpha), which yields a square-root-n concentration.

**Conformity score families covered (Section 2):**
- **Naive softmax score** (Section 1, Figure 2): s = 1 - f-hat(X)_Y. Smallest average set size but undercovers hard inputs and overcovers easy ones.
- **Adaptive Prediction Sets, APS** (Section 2.1, following Romano-Sesia-Candes and Angelopoulos-Bates-Malik-Jordan): cumulative softmax mass up to the true class. Eq. 3.
- **RAPS** (Regularized APS) is referenced as the practical improvement layered on APS, citing reference [4].
- **Conformalized Quantile Regression, CQR** (Section 2.2, Romano-Patterson-Candes): score s(x,y) = max(t-hat_{alpha/2}(x) - y, y - t-hat_{1-alpha/2}(x)); prediction set is [t-hat_{alpha/2}(x) - q-hat, t-hat_{1-alpha/2}(x) + q-hat] (Eq. 4). Inherits asymptotic conditional coverage from quantile regression.
- **Normalized residual / scalar uncertainty score** (Section 2.3): s(x,y) = |y - f-hat(x)| / u(x); set is f-hat(x) +/- u(x) * q-hat (Eq. 5). u(x) can be a learned standard deviation, ensemble variance, MC-dropout variance, input-perturbation variance, or adversarial sensitivity.
- **Conformalized Bayes / density-based score** (Section 2.4, Hoff [11]): s(x,y) = -f-hat(y|x); the set is the superlevel set of the posterior predictive (Eq. 6). Bayes-optimal among 1-alpha sets under the modeling assumptions.

**Extensions (Section 4):**
- **Group-balanced conformal** (4.1, Proposition 1): stratify scores by discrete group g, compute q-hat^{(g)} per group; satisfies P(Y_test in C(X_test) | X_{test,1} = g) >= 1 - alpha.
- **Class-conditional conformal** (4.2, Proposition 2): stratify by Y; provisional inclusion uses q-hat^{(y)}.
- **Conformal Risk Control, CRC** (4.3, Theorem 2, Angelopoulos-Bates-Fisch-Lei-Schuster [17]): for any bounded, monotone-in-lambda loss l(C_lambda(x), y) <= B, choose lambda-hat = inf{lambda : R-hat(lambda) <= alpha - (B - alpha)/n}; then E[l(C_{lambda-hat}(X_test), Y_test)] <= alpha. Recovers miscoverage when l is the indicator. Eq. 12. The Learn then Test framework (Appendix A, [18]) handles non-monotone risks via FWER-controlled p-value testing.
- **Outlier detection** (4.4, Proposition 3): unsupervised score s(X), threshold at d(n+1)(1-alpha)e-th score; controls Type-I error P(C(X_test)=outlier) <= alpha.
- **Weighted conformal under covariate shift** (4.5, Theorem 3, Tibshirani-Foygel Barber-Candes-Ramdas [25]): weights w(x) = dP_test(x)/dP(x); reweighted quantile of calibration scores.
- **Conformal under distribution drift** (4.6, Theorem 4, Foygel Barber-Candes-Ramdas-Tibshirani [26]): weighted quantile with non-uniform weights w_i; coverage at least 1 - alpha - 2 * sum(w-tilde_i * epsilon_i) where epsilon_i is TV distance to the test distribution. Effective sample size n_eff = (sum w_i)^2 / sum w_i^2 governs variance.

**Full and cross-conformal (Section 6):** full conformal (Theorem 5) loops over candidate y, refits f-hat^y on the augmented set of size n+1, and includes y if s_{n+1}^y <= q-hat^y; computationally expensive but data-efficient. **Cross-conformal** (Vovk [41]) and **CV+ / Jackknife+** (Barber-Candes-Ramdas-Tibshirani [42]) interpolate between split and full, using K-fold or leave-one-out fits while reusing all data.

## Datasets / experiments
The monograph is pedagogical with worked examples rather than benchmark tables.
- **ImageNet** classification (Figures 1, 2, 17): fox-squirrel set examples; selective classification at alpha = 0.1.
- **MS COCO** multilabel classification (Section 5.1, Figure 13) with FNR control at alpha = 0.1 via CRC.
- **Gut polyp tumor segmentation** (Section 5.2, Figure 14) with pixel-FNR control at alpha = 0.1.
- **Yandex Weather Prediction / Shifts Project** (Section 5.3, Figure 15): time-series temperature with K-window weighted conformal, ensemble of 10 CatBoost models as f-hat and ensemble variance as u-hat.
- **Jigsaw Multilingual Toxic Comment Classification (WILDS)** (Section 5.4, Figure 16): outlier detection with Detoxify BERT; at alpha = 0.1 the system flags 70% of true toxic comments.
- **Calibration set size guidance** (Section 3.2, Table 1): at alpha = 0.1 and delta = 0.1, n = 22 gives epsilon = 0.1 slack, n = 1000 typically yields coverage in [0.88, 0.92], n = 244,390 for epsilon = 0.001.

## Metrics
Marginal coverage (Eq. 1). Conditional coverage P(Y_test in C(X_test) | X_test) >= 1 - alpha (Eq. 7); known to be impossible distribution-free in the general continuous case (Vovk [14], Lei-Wasserman [3]). **Feature-Stratified Coverage (FSC)** and **Size-Stratified Coverage (SSC)** as the two practical proxies (Section 3.1). Set size histograms and spread for adaptivity. Empirical coverage averaged over R random calibration/validation splits with cached scores (Section 3.3, Figure 12). Beta-distribution-based correctness checks via Appendix C.

## What we steal for ConfPert
1. The **split-conformal procedure (Section 1, Eq. 2; Figure 2 code template)** is the literal backbone of ConfPert: each of our six distributional discrepancies (KS, W1, energy distance, MMD-RBF, bimodality match, variance ratio) instantiates a different score function s(x, y_pred-distribution, y_true-distribution) over predicted vs. observed cell-population responses, then we calibrate via the ceiling((n+1)(1-alpha))/n quantile.
2. **Score-function design vocabulary (Section 2):** ConfPert's six discrepancies map cleanly onto the score families catalogued here. KS and W1 are CDF-based residuals; energy distance and MMD-RBF are kernel/density-based scores in the spirit of Section 2.4's conformalized-Bayes treatment; variance ratio is a normalized scalar-uncertainty score in the style of Section 2.3.2 with u(x) = sigma-hat(x); bimodality match is a structured score akin to APS (Section 2.1).
3. **Conformal Risk Control (Section 4.3, Theorem 2, Eq. 11-12)** is the lift for ConfPert knockout K2 (drug-screening triage). We replace miscoverage with a triage-relevant loss (false-negative rate on hits, expected efficacy gap, etc.) and tune lambda via Eq. 12. Since drug-screening triage losses are monotone in the inclusion threshold, Theorem 2 applies directly without needing Learn then Test.
4. **Weighted conformal (Section 4.5, Theorem 3) and distribution-drift (Section 4.6, Theorem 4)** are the basis for K1's robustness story when calibration cell lines differ from deployment cell lines (covariate shift in cell-state covariates).
5. **FSC and SSC metrics (Section 3.1)** are how we evaluate ConfPert's per-perturbation, per-cell-type, and per-magnitude conditional validity for K3 (PRISM downstream).
6. **Calibration-set sizing (Section 3.2, Table 1, Beta(n+1-l, l) result)** sets the experimental design floor: n >= 1000 calibration perturbations per ConfPert deployment, with explicit slack-vs-confidence accounting.

## What we wrap
N/A in this paper. It is a pedagogical monograph, not a model. We wrap downstream cell-population perturbation predictors (e.g., scGen-, CPA-, Geneformer-style outputs) using the procedures synthesized here.

## Failure modes / caveats
- **Heteroscedasticity** breaks naive residual scores; this motivates CQR (Section 2.2) and normalized scalar-uncertainty scores (Section 2.3). Naive |y - f-hat(x)| produces fixed-width intervals that overcover easy inputs and undercover hard ones.
- **Conditional vs marginal coverage gap** (Section 3.1, Figure 10): marginal coverage is satisfied even when all errors concentrate in one subgroup. Conditional coverage is provably impossible distribution-free for arbitrary continuous distributions (Vovk [14], Lei-Wasserman [3]). Vanishing-width conditional intervals require effective support of X_test smaller than sample-size squared (Sesia-Candes [88]).
- **Distribution shift** breaks Theorem 1's exchangeability premise: covariate shift requires weighted CP (Theorem 3) and known likelihood ratio; unknown drift only gives the slack bound 2 * sum(w-tilde_i * epsilon_i) in Theorem 4.
- **Effective sample size collapse** under aggressive weighting: n_eff = (sum w_i)^2 / sum w_i^2 can be tiny if too many weights are small, inflating coverage variance as 1/sqrt(n_eff).
- **Score informativeness is not free**: Theorem 1 holds even if the score is random noise, but the resulting sets are useless. Practical utility is determined entirely by the score function.
- **Adaptivity is not implied** by the coverage guarantee. SSC/FSC must be checked separately.
- **Discreteness and ties** require minor randomization to achieve the upper bound in Eq. 1; usually ignored in practice.

## Code URLs
Companion code and Jupyter notebooks: https://github.com/aangelopoulos/conformal-prediction. Author blog and tutorial video: http://angelopoulos.ai/blog/posts/gentle-intro/. Companion library cited in the monograph: MAPIE (scikit-learn-compatible). Awesome Conformal Prediction repository ([130] in the bibliography) curates current resources.

## Verbatim quotes worth keeping
1. "Conformal prediction is a user-friendly paradigm for creating statistically rigorous uncertainty sets/intervals for the predictions of such models. Critically, the sets are valid in a distribution-free sense: they possess explicit, non-asymptotic guarantees even without distributional assumptions or model assumptions."
2. "Although the guarantee always holds, the usefulness of the prediction sets is primarily determined by the score function."
3. "Conformal prediction can be seen as a method for taking any heuristic notion of uncertainty from any model and converting it to a rigorous one."
4. "It is extremely important to keep in mind that the conformal prediction procedure with the smallest average set size is not necessarily the best. A good conformal prediction procedure will give small sets on easy inputs and large sets on hard inputs in a way that faithfully reflects the model's uncertainty."
5. "In the most general case, conditional coverage is impossible to achieve."
