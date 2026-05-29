# Bates, Angelopoulos, Lei, Malik, Jordan 2021 - Risk-Controlling Prediction Sets (RCPS)

**Venue:** JACM 2021. arXiv:2101.02703 (v3, 4 Aug 2021).

## Claim
RCPS extends conformal prediction from binary coverage (probability that Y is in C(X)) to general bounded, monotone (set-nested) loss functions, giving distribution-free finite-sample high-probability control of the population risk R(T) = E[L(Y, T(X))] at a user-chosen level alpha with confidence 1 - delta, by calibrating a one-dimensional threshold lambda on a holdout set via an upper confidence bound.

## Method

**Definition 1 (verbatim, p.2):** "Let T be a random function taking values in the space of functions X to Y' (e.g., a functional estimator trained on data). We say that T is a (alpha, delta)-risk-controlling prediction set if, with probability at least 1 - delta, we have R(T) <= alpha."

**Setting (Sec 2.1):** i.i.d. data split into training and calibration of size n. Nested family {T_lambda}_{lambda in Lambda} satisfying (Eq. 1): "lambda_1 < lambda_2 implies T_{lambda_1}(x) subset T_{lambda_2}(x)." Loss satisfies (Eq. 2): "S subset S' implies L(y, S) >= L(y, S')." Risk R(T) = E[L(Y, T(X))]; assume some lambda_max with R(lambda_max) = 0.

**UCB calibration (Eq. 3-4):** Assume pointwise (1 - delta) UCB R-hat^+(lambda) with P(R(lambda) <= R-hat^+(lambda)) >= 1 - delta. Choose lambda-hat = inf{lambda : R-hat^+(lambda') < alpha for all lambda' >= lambda}.

**Theorem 1 (verbatim, p.5):** "Suppose (3) holds pointwise for each lambda, and that R(lambda) is continuous. Then for lambda-hat chosen as in (4), P(R(T_{lambda-hat}) <= alpha) >= 1 - delta. That is, T_{lambda-hat} is a (alpha, delta)-RCPS."

Proof (Appendix Thm A.1): define lambda^* = inf{lambda : R(lambda) <= alpha}. If R(lambda-hat) > alpha then monotonicity forces lambda-hat < lambda^*, which by definition forces R-hat^+(lambda^*) < alpha; but R(lambda^*) = alpha, so the pointwise UCB fails with probability at most delta. Monotonicity of R(lambda) avoids any uniform/union bound over Lambda.

**Concentration inequalities (Sec 3):**
- Simple Hoeffding (Eq. 5): R-hat + sqrt(log(1/delta)/(2n)). Bounded loss <= 1.
- Hoeffding-Bentkus (Thm 3): min of tighter Hoeffding (Prop 3, h_1(t; R) = t log(t/R) + (1-t) log((1-t)/(1-R))) and Bentkus (Prop 4: P(R-hat <= t) <= e * P(Binom(n, R) <= ceil(nt))); invert via Prop 2.
- Waudby-Smith-Ramdas (WSR, Prop 5, Thm 4): hedged-capital martingale K_i = prod (1 - nu_j(L_j - R)) with predictable nu_i = min(1, sqrt(2 log(1/delta) / (n * sigma-hat^2_{i-1}))). Variance-adaptive; recommended for non-binary bounded losses.
- Pinelis-Utev (Thm 5): unbounded losses with bounded coefficient of variation.
- CLT (Eq. 10, Thm 6): R-hat + Phi^{-1}(1 - delta) * sigma-hat / sqrt(n); asymptotic only.

**Calibration size (Sec 3.4):** "around 1,000 calibration points suffice" for alpha = 0.1; "a few thousand" for alpha = 0.01; "about 10,000" for alpha = 0.001. Dependence on delta is log(1/delta).

**Algorithm 1 (Greedy Sets, p.12):** Constructs the nested family T_lambda(x) by sweeping a threshold zeta from a large value down to -lambda in steps d-zeta, at each step adding to T the y' in the complement with conditional risk density rho-hat_x(y', T) > zeta. Output: T_lambda(x). Theorem 7 (and its generalization Theorem 8 to losses of the form L(y; S) = integral over S^c of ell(y, z) d-mu(z)) shows the greedy sets are optimal: among any T' with R(T') <= R(T_lambda), E[|T_lambda(X)|] <= E[|T'(X)|], i.e., greedy gives smallest expected set size.

## Datasets / experiments
1. **ImageNet classification with class-varying loss (Sec 5.1):** L_y ~iid Unif(0, 1), pretrained ResNet-152, 30,000 calibration / 20,000 test, alpha = 0.1, delta = 0.1, 100 splits. Fig. 7: realized risk concentrated near 0.1 (controlled), set sizes typically 2-6.
2. **MS COCO multi-label classification (Sec 5.2):** loss = false negative rate 1 - |y intersect S| / |y|, TResNet base, 4,000 calibration / 1,000 test, alpha = 0.1, delta = 0.1, 1,000 splits. Fig. 9: RCPS controls risk at ~0.1; conformal baseline overcovers (risk ~0.025) with much larger sets.
3. **ImageNet hierarchical classification (Sec 5.3):** WordNet hierarchy, depth D = 14, ResNet-18, alpha = 0.05, delta = 0.1. Fig. 11 shows risk near 0.05 with low tree-height predictions.
4. **Polyp segmentation (Sec 5.4):** Kvasir + Hyper-Kvasir + CVC-ColonDB + CVC-ClinicDB + ETIS-Larib (1,781 segmented polyps), 1,000 calibration / rest test, PraNet base, alpha = 0.1, delta = 0.1; per-object FNR loss; Fig. 13 shows controlled risk and set sizes comparable to average polyp size.
5. **Protein structure prediction on CASP-13 (Sec 5.5):** AlphaFold v1 distograms, l_1 projective distance loss, 35 calibration / 36 test, alpha = 2 Angstrom, delta = 0.1, CLT bound only (small sample). Fig. 15: realized risk near 2 A.

(No PASCAL VOC; image segmentation experiment is the polyp pipeline.)

## Metrics
Risk = E[L(Y, T(X))]. Coverage is the special case L(Y, S) = 1{Y not in S}. Concrete losses studied: class-varying mistake cost L_y * 1{y not in S}, multi-label false negative rate, hierarchical distance d_H scaled by tree depth D, image segmentation per-object missed-pixel fraction sum_{y' in h(y)} |y' \ S| / |y'| / |h(y)|, l_1 projective distance for protein distograms, ranking 0-1 sign-mismatch loss for U-statistic-based extensions (Prop 7, Theorem 9). Reported quantities are realized risk distribution over splits and prediction set size distribution.

## What we steal for ConfPert (K2)

K2 maps directly onto Theorem 1:
- T_lambda(x) is the emit-or-abstain predictor: emit a perturbation prediction if calibrated energy distance to truth < lambda, else abstain (encoded as the largest set).
- Loss L(Y, T_lambda(X)) = clipped excess energy distance above tau when emitting, zero when abstaining; bounded to [0, B] satisfies Sec 3.1.
- Nesting (Eq. 1) and monotone risk hold by construction: larger lambda emits more often, weakly increasing risk.
- Calibration: holdout of triage candidates, energy distance against held-out true single-cell perturbation distributions per drug, lambda grid, lambda-hat per Eq. 4.
- Concentration: WSR for non-binary bounded loss; HB if loss is near-binary; CLT only as small-sample fallback.

The load-bearing technical move we lift is the monotonicity-driven proof of Thm A.1: no uniform bound over Lambda is needed because R(lambda) is monotone, so the search collapses to the boundary lambda^*. Our energy-distance risk is monotone in the emission threshold by construction.

## What we wrap

The RCPS framework wraps any sample-producing perturbation predictor f-hat (e.g., a single-cell diffusion or flow-matching model). K2 statement: at user-chosen target risk delta, only emit a perturbation prediction if the calibrated population-level energy distance is below tau; the (alpha, delta)-RCPS guarantee is given by Theorem 1 applied to the bounded energy-distance excess loss, with WSR as the concentration tool and Eq. 4 as the threshold-selection rule. Per Remark 2, the predictor f-hat may be trained on shifted data; only calibration and deployment must be exchangeable.

## Failure modes / caveats
- Bounded-loss: energy-distance excess must be clipped to [0, B]; choice of B trades tightness vs information loss. Pinelis-Utev or CLT relax this but require c_v bound or asymptotics.
- Monotonicity of R(lambda) is essential; must be checked or enforced by construction.
- Continuity of R(lambda) assumed; Remark 3: removable with minor modifications.
- Conservatism: Hoeffding loose; WSR tighter but still finite-sample conservative; CLT tightest but asymptotic only.
- Grid search over lambda; per-lambda recompute is O(n) per grid point, fine for K2.
- Distribution shift between calibration and deployment breaks the guarantee.
- Non-i.i.d. data (correlated cells from the same donor) violates Thm 1's hypothesis.

## Code URLs
- Project page: https://angelopoulos.ai/blog/posts/rcps/
- arXiv: https://arxiv.org/abs/2101.02703
- GitHub repo (linked from project site, "public GitHub repository" in Sec 3 and Sec 5): https://github.com/aangelopoulos/rcps

## Verbatim quotes worth keeping

1. (Definition 1, p.2) "We say that T is a (alpha, delta)-risk-controlling prediction set if, with probability at least 1 - delta, we have R(T) <= alpha."

2. (Theorem 1, p.5) "Suppose (3) holds pointwise for each lambda, and that R(lambda) is continuous. Then for lambda-hat chosen as in (4), P(R(T_{lambda-hat}) <= alpha) >= 1 - delta."

3. (p.5, on why no uniform bound is needed) "Note that we are able to turn a pointwise convergence result into a result on the validity of a data-driven choice of lambda. This is due to the monotonicity of the risk function; without the monotonicity, we would need a uniform convergence result on the empirical risk in order to get a similar guarantee."

4. (Remark 2, p.5) "UCB calibration gives an RCPS even if the data used to fit the initial predictive model comes from a different distribution. The only requirement is that the calibration data and the test data come from the same distribution."

5. (Sec 7 Discussion, p.22) "Concentration is a more general tool and can apply to a wider range of problems... in contrast to the standard train/validation/test split paradigm which only estimates global uncertainty (in the form of overall prediction accuracy), RCPS allow the user to automatically return valid instance-wise uncertainty estimates for many prediction tasks."
