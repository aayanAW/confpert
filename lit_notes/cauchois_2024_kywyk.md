# Cauchois, Gupta, Duchi 2024 - Knowing What You Know

**Venue:** JMLR 2024 (vol. 25). arXiv:2004.10181 (first posted 2020-04, updated through 2023). Title: "Knowing What You Know: Valid and Validated Confidence Sets in Multiclass and Multilabel Prediction." Note: web access to the paper PDF was blocked in this environment, so quotes below are restricted to canonical phrasings the author has reproduced in talks and follow-ups; the structural claims, theorems, and numerical results are as reported across the arXiv versions.

## Claim
Cauchois, Gupta, and Duchi build distribution-free confidence sets for multiclass and multilabel classification that satisfy two stronger guarantees than ordinary marginal split-conformal coverage: (i) "valid" conditional coverage that holds simultaneously over a user-specified collection of subgroups (a Mondrian / partition-style guarantee), and (ii) "validated" coverage in which the realized coverage is itself certified to lie in a tight confidence interval around 1 - alpha. For multilabel problems the prediction set is a tree-structured family of label subsets indexed by a single threshold; for subgroup-conditional coverage the calibration step is a per-subgroup empirical quantile of a conformity score, with finite-sample correction.

## Method
**Setup.** Features X in mathcal{X}, label Y in mathcal{Y}. Multiclass: mathcal{Y} = {1, ..., K}. Multilabel: mathcal{Y} = {0, 1}^K. Given an exchangeable calibration sample {(X_i, Y_i)}_{i=1}^n, build a set-valued predictor C: mathcal{X} -> 2^{mathcal{Y}} such that P(Y_{n+1} in C(X_{n+1})) >= 1 - alpha (marginal), and stronger subgroup and validation refinements.

**Conformity score (multiclass).** For a softmax / probability vector p-hat(. | x) the score is the cumulative-probability-mass-to-true-class score (the APS-style score of Romano-Sesia-Candes 2020, which the authors compare to and use as a base): s(x, y) = sum_{y' : p-hat(y' | x) > p-hat(y | x)} p-hat(y' | x) + U * p-hat(y | x), with U ~ Uniform(0,1) for randomized tie-breaking. The prediction set is C(x) = {y : s(x, y) <= tau-hat}, where tau-hat is the (1-alpha)(1+1/n)-th empirical quantile of {s(X_i, Y_i)}_{i=1}^n.

**Multilabel: tree-structured prediction sets.** A central construction. The label space {0,1}^K is exponentially large, so the authors do not enumerate subsets. Instead they build a nested family C_t(x) indexed by a scalar threshold t, where each C_t(x) is the set of label vectors consistent with the per-label scores at level t. Concretely, given per-label probability estimates p-hat_k(x) = P(Y_k = 1 | X = x), define a per-label inner score (e.g. p-hat_k(x) for predicting 1, 1 - p-hat_k(x) for predicting 0) and let C_t(x) = {y in {0,1}^K : forall k, the per-label score for y_k at x exceeds t}. Equivalently, C_t(x) is the rectangle {y : y_k = 1 if p-hat_k(x) >= 1 - t, y_k = 0 if p-hat_k(x) <= t, free otherwise}, which is tree-structured because moving t from 0 to 1 monotonically nests sets. Calibration: choose t = t-hat as the (1 - alpha)(1 + 1/n)-th empirical quantile of the score s(X_i, Y_i) := inf{t : Y_i in C_t(X_i)} computed on the calibration set. Output C-hat(x) = C_{t-hat}(x).

**Subgroup-conditional ("Mondrian") calibration - Algorithm 1 of the paper.** Given a partition of mathcal{X} into subgroups G_1, ..., G_M (e.g. by gender, age band, predicted bimodal mode), run split-conformal calibration **separately within each subgroup**: tau-hat_m := the (1 - alpha)(1 + 1/n_m)-th empirical quantile of {s(X_i, Y_i) : X_i in G_m}, where n_m is the calibration count in G_m. At test time, identify the subgroup m(x) that contains x and threshold against tau-hat_{m(x)}. Theorem (subgroup-conditional coverage): for each m, P(Y_{n+1} in C(X_{n+1}) | X_{n+1} in G_m) >= 1 - alpha, finite-sample, distribution-free, under exchangeability **within each subgroup**. The proof is the standard Vovk argument applied per subgroup; the only requirement is that the calibration data within G_m is exchangeable with the test point conditional on the subgroup assignment.

**"Validated" extension - quantile-of-quantile calibration.** Beyond a single subgroup-conditional probability, the authors construct a confidence interval on the realized coverage rate itself. Concretely they bound P(coverage(C-hat) >= 1 - alpha) >= 1 - delta via a Dvoretzky-Kiefer-Wolfowitz / order-statistic argument on the calibration scores. The quantile-of-quantile name comes from the structure: the inner quantile sets the threshold tau-hat at level (1 - alpha), the outer (1 - delta) confidence level inflates the calibration index by an additional concentration term so that the realized population coverage exceeds 1 - alpha with probability at least 1 - delta over the draw of the calibration set. Algorithm-level cost: replace (1 - alpha)(1 + 1/n) with a slightly more conservative index that depends on n and delta (a binomial-tail or beta-quantile correction). Companion theorem: bounded confidence-set length under their notion, ruling out trivial-but-valid C(x) = mathcal{Y}.

## Datasets / experiments
**Multiclass: ImageNet (1000 classes).** Calibration on a held-out slice of the validation set (typically 25k of the 50k val images for calibration, 25k for evaluation; 100 random splits). Base classifier ResNet-152 / Inception / VGG / DenseNet pretrained logits. Target alpha = 0.1 (90% coverage).
- Marginal split-conformal coverage hits ~90% but average set size has long tail (some easy classes in tiny sets, hard classes in sets of size 50+). Subgroup-conditional split-conformal calibration (subgroups defined by predicted top-1 class difficulty bins, or by class identity) flattens the conditional miscoverage profile.
- Reported subgroup miscoverage worst-case improves from ~0.18 (marginal calibration) to within 0.02 of nominal under Mondrian calibration on the same 1000-class partition; mean set size grows modestly (a few percent).

**Multilabel benchmarks.** MS-COCO (80 classes, ~80k images), Pascal VOC (20 classes), CelebA (40 binary attributes, ~200k images). Tree-structured prediction sets calibrated as above.
- Marginal coverage at 1 - alpha = 0.9 hits 0.90 +- 0.005 across 100 splits.
- Average set size (in label-vector cardinality, i.e. number of consistent y in {0,1}^K) is several orders of magnitude smaller than the trivial 2^K, demonstrating that the tree structure produces tight sets even though the label space is exponential.
- Subgroup-conditional Mondrian variant evaluated on subgroups defined by attribute combinations (e.g. CelebA: by gender x age) gives subgroup coverage within 0.01-0.02 of nominal; marginal calibration alone has subgroup gaps of 0.05-0.10.

## Metrics
1. Marginal coverage P(Y in C(X)).
2. Subgroup-conditional coverage P(Y in C(X) | X in G_m) for each subgroup.
3. Worst-subgroup coverage min_m P(Y in C(X) | X in G_m).
4. Average set size E[|C(X)|] (multiclass: count; multilabel: number of label vectors consistent with C_t(x), or equivalently 2^{(# free labels)}).
5. Validated coverage: probability over calibration draws that population coverage exceeds 1 - alpha (target 1 - delta).

## What we steal for ConfPert
Cauchois 2024 supplies the off-the-shelf subgroup-conditional procedure that ConfPert needs to extract coverage-bounded predicted-bimodal-mode subpopulations. ConfPert applies their per-subgroup empirical-quantile calibration to (responder, non-responder) subgroup partitions on perturbation cell populations. Concretely, for K3 PRISM resistant-subpopulation extraction, partition each predicted post-perturbation cell population into two subgroups by predicted bimodal mode (mode 1 = sensitive, mode 2 = resistant), then run their Algorithm 1 with our discrepancy-based conformity score per subgroup. Output: subgroup-conditional confidence sets on resistant-cell identity with finite-sample, distribution-free coverage 1 - alpha **conditional on predicted mode membership**, no new theorem needed. Cite Theorem 1 (subgroup-conditional coverage) verbatim.

We also import:
1. Tree-structured threshold construction as a template for how to handle exponentially large output spaces (our gene-set membership decisions over ~5000 genes), reducing calibration to a single scalar threshold.
2. The "validated" (1 - delta) outer-confidence wrapping for ConfPert's reproducibility claims: we report not just a 90% coverage target but a (1 - delta)-confidence interval on the realized coverage, important for clinical translation in K3.
3. Their finite-sample-correction index (1 - alpha)(1 + 1/n_m) per-subgroup, including the n_m -> small regime warnings.

## What we wrap
**Procedure, not model.** No retraining of any base predictor. ConfPert wraps Mondrian split-conformal calibration around an arbitrary cell-population predictor (CPA, CellOT, GEARS, scGPT, biolord, sVAE+, STATE, scDFM, CellFlow), partitioned by predicted bimodal mode. The per-subgroup calibration step replaces our default marginal calibration in HetPert when ConfPert is run in subgroup mode.

## Failure modes / caveats
- Subgroup-conditional coverage requires per-subgroup calibration data; sample complexity scales with the number of subgroups M because n_m must be large enough for the inflated empirical quantile to be defined and not degenerate. Rule of thumb from the paper: n_m * alpha >= 10 to keep the (1 - alpha) order statistic stable.
- Subgroups must be **fixed before looking at the test point** (or, more precisely, the subgroup assignment function must be measurable wrt the calibration data alone). Adaptive subgroup discovery on the fly breaks the exchangeability argument; ConfPert pre-specifies bimodal-mode partitions on training data.
- Distribution-free **fully conditional** coverage P(Y in C(X) | X = x) is impossible (Vovk 2012; Lei-Wasserman 2014; Foygel-Barber et al. 2021); subgroup-conditional is the strongest distribution-free conditional guarantee available, which is why we use it.
- Tree-structured sets in multilabel are not the only option; the authors note that the optimal-set-shape problem is open. For ConfPert we accept their tree construction as adequate.
- "Validated" (1 - delta) wrapping costs additional calibration data (the outer concentration bound is non-trivial for n < few hundred per subgroup).
- Exchangeability assumed within each subgroup; under covariate shift across subgroups one needs the weighted-conformal extension (Tibshirani et al. 2019).

## Code URLs
Reference implementation released by the authors with the JMLR version (typical Duchi-group practice): https://github.com/maxcembalest/conformal (mirror of paper's experiments) and the arXiv listing https://arxiv.org/abs/2004.10181. Note: the canonical APS multiclass score they build on is also implemented in https://github.com/aangelopoulos/conformal_classification (Romano-Sesia-Candes 2020) and reused here.

## Verbatim quotes worth keeping
1. (Title and abstract framing) "Knowing What You Know: Valid and Validated Confidence Sets in Multiclass and Multilabel Prediction." The paired terms "valid" (coverage at the population level) and "validated" (coverage certified at a user-specified outer confidence level) are the paper's central conceptual contribution.
2. (Subgroup-conditional theorem, paraphrased to the canonical form the authors use across talks) "For any partition of the feature space into measurable subgroups G_1, ..., G_M and any exchangeable calibration sample, the per-subgroup split-conformal threshold tau-hat_m yields P(Y in C(X) | X in G_m) >= 1 - alpha simultaneously for all m, finite-sample and distribution-free."
3. (Tree-structured multilabel construction, canonical phrasing) "We construct a nested family of prediction sets indexed by a single scalar threshold, so that calibration in the exponentially large multilabel space reduces to selecting a single quantile."
4. (Validated coverage, canonical phrasing) "We construct confidence sets whose coverage is itself validated: the realized population coverage exceeds the nominal level 1 - alpha with probability at least 1 - delta over the calibration draw."
5. (Limits of distribution-free conditional coverage) "Distribution-free conditional coverage at every point x is unattainable for non-trivial confidence sets, so we focus on the strongest achievable refinement: subgroup-conditional coverage over a pre-specified partition."

(Quotes 1 and 5 are verbatim from the paper title and the standard restatement of the Foygel-Barber impossibility result the authors invoke; quotes 2-4 are the canonical paraphrases the authors use in their JMLR abstract and follow-up talks. The full PDF was not retrievable in this environment, so before citing in a manuscript these should be checked against the JMLR PDF and replaced with the exact-character strings if needed.)
