# ConfPert Preregistration

**Project:** ConfPert: model-agnostic conformal calibration for distributional single-cell perturbation predictors.
**Author:** Aayan Alwani.
**Document created:** 2026-05-03.
**Lock target:** committed via git before any Phase 1 framework code, baseline runs, or hyperparameter selection.
**Purpose:** lock K1, K2, K3 success criteria, K1 baseline-freeze conditions, and noise-variant reporting requirements before any tuning. Reviewers cannot challenge pre-registration that is git-stamped before the first run.

This document supersedes any later informal hypothesis revisions. Post-hoc reframing of any criterion below is treated as exploratory and reported as such, not as a primary finding.

---

## K1 lock conditions

**Baselines.** All eight wrappable predictors are trained or loaded once with default hyperparameters from each method's published or repo-default configuration. No predictor-specific hyperparameter tuning before lock-in. Results frozen in `baselines/results.json` with a SHA-256 hash recorded in this document the day Phase 1 begins.

The eight predictors, stratified by parameter count:

1. Mean predictor (~0 params)
2. Bilinear ridge per Ahlmann-Eltze 2025 eq. 1, 3 (~10^4)
3. scGen (~10^6)
4. CPA / cpa-tools (~5 x 10^6)
5. biolord (~2 x 10^7)
6. sVAE+ / SAMS-VAE (~10^7)
7. GEARS in uncertainty mode (~5 x 10^7)
8. STATE with SE-600M public checkpoint (~6 x 10^8)

scDFM, CellFlow, Squidiff, PerturbDiff, CellOT, Unlasting are NOT in the K1 sweep. They are cited as related work only. Reasons: per-dataset retrain compute cost, K2 axis-clarity preservation, and saturated-subfield rationale documented in the synthesis.

**Datasets.** Norman 2019, Replogle 2022 K562 essential, Replogle 2022 RPE1 essential, Tahoe-100M cross-cell-line subset. Adamson 2016 in appendix only, not in main K1 sweep. Dataset preprocessing follows the loaders specified in the existing HetPert codebase plus a Tahoe-100M loader to be implemented at Phase 1 start, lock the loader version with a commit hash before first run.

**Splits.** Three per dataset:
- Within-perturbation 70/30 cell-level split.
- Held-out perturbation 5-fold cross-validation.
- Cross-cell-line: K562 calibration, RPE1 test (where applicable to Replogle and Tahoe-100M subset).

**Discrepancies.** Six per (predictor, dataset, split) cell:
- Kolmogorov-Smirnov (per-gene 1D, mean over genes)
- Wasserstein-1 (per-gene 1D, mean over genes)
- Energy distance (per-gene, mean over genes)
- MMD-RBF (per-population, median-bandwidth heuristic)
- Bimodality coefficient match (SAS, threshold $b > 5/9$)
- Variance ratio (per-gene, mean over genes)

**Coverage targets.** $\alpha \in \{0.05, 0.10, 0.20\}$, three nominal levels.

**Noise variants for point-estimate predictors (Mean, Bilinear ridge).** All four reported, no cherry-picking:
- Variant A: no noise (point mass).
- Variant B: isotropic Gaussian, pooled per-pert std.
- Variant C: per-gene marginal Gaussian.
- Variant D: full empirical-covariance Gaussian.

**Reporting.** Every (predictor, dataset, split, discrepancy, $\alpha$, noise-variant for predictors 1-2) cell is reported. No selective reporting. No hypothesis-driven cell selection.

---

## K2 hypotheses, both pre-registered

### H1 (capacity hypothesis)

**Claim:** calibration error increases with parameter count at fixed nominal coverage.

**Operationalization:** Spearman rank correlation between $\log(\text{params})$ and the absolute calibration deviation $|1 - \alpha - \text{achieved coverage}|$ across the eight predictors at each $\alpha \in \{0.05, 0.10, 0.20\}$.

**Success criterion:** Spearman $\rho > 0.5$ with $p < 0.05$ (permutation test on predictor labels with $\geq 10{,}000$ permutations) on at least three of the four datasets, computed at the average over the three $\alpha$ levels.

### H1b (data hypothesis)

**Claim:** calibration error decreases with training-set size at fixed nominal coverage.

**Operationalization:** Spearman rank correlation between $\log(N_{\text{training cells per perturbation, averaged}})$ and the absolute calibration deviation, expected sign negative.

**Success criterion:** Spearman $\rho < -0.5$ with $p < 0.05$ (permutation test, $\geq 10{,}000$ permutations) on at least three of the four datasets at average over $\alpha$ levels.

### Outcome dispositions

- Both H1 and H1b land $\to$ richer story: capacity hurts AND data helps. Combined design directive in Section 7.
- Only H1 lands $\to$ Ahlmann-Eltze-style headline: bigger is worse for calibration.
- Only H1b lands $\to$ data-scale-driven calibration improvement; pairs with VCC's Tahoe-100M push.
- Neither lands $\to$ null result. Publish honestly. Lean on K1 + K3 plus the design-directive paragraph for "training-objective-family analysis is the next research direction."

### Failure plan

If neither $\rho$ exceeds threshold, no post-hoc reframing to a third axis as a pre-registered finding. Third-axis patterns (e.g., training-objective-family clustering) report as exploratory observation in the discussion section, explicitly labeled "post-hoc, not pre-registered."

---

## K3 success criterion (BH-FDR-corrected)

**Pipeline:**

1. Train ConfPert wrappers on Tahoe-100M observational arm. Predict cell-population responses for the 49 selective-and-predictable non-oncology compounds in Corsello 2020 PRISM (those with Pearson $r > 0.4$ in their primary screen).
2. For each compound, extract calibrated bimodal subpopulations using Cauchois 2024 subgroup-conditional sets on the predicted-bimodal-mode partition.
3. Cross-reference resistant-subpopulation gene signatures against:
   - DepMap CRISPR essentiality (compound x gene grid).
   - PRISM secondary screen viability AUC (compound x cell-line grid).
   - Hallmark MSigDB pathways (compound x pathway grid).
4. Repeat for uncalibrated baselines (raw GEARS predict-mean, raw STATE without conformal calibration, raw CPA without conformal calibration).

**Multiple-comparisons correction:** Benjamini-Hochberg FDR at $q = 0.05$ over the joint compound x Hallmark-pathway grid for the enrichment test, plus a separate BH at $q = 0.05$ over the compound x DepMap-gene grid for the dependency-alignment test. A signature counts only if both corrected tests pass for that compound.

**Success criterion:** at least three drug-cell-line pairs where ConfPert's calibrated bimodal subpopulation produces:
- Hallmark enrichment $q < 0.05$ (BH over the full compound x pathway grid)
- DepMap selective-dependency alignment $q < 0.05$ (BH over the compound x gene grid)
- AND uncalibrated baselines fail at least two of those three at the same corrected thresholds.

**Failure plan:** if K3 fails the BH-corrected criterion, report honestly. Do not switch to looser uncorrected thresholds post-hoc. K1 plus K2 carry the paper.

---

## Hash record

Phase 1 lock-in must record the following hashes in this document or a follow-on `lock_record.md`:

- Git commit hash of this preregistration.md at the moment Phase 1 begins.
- SHA-256 of `baselines/results.json` after first lock-in.
- Git commit hash of the dataset loader code at lock-in time.
- Git commit hash of `ArcInstitute/cell-eval` at the time of the audit, plus the present/missing discrepancy list in `baselines/cell_eval_audit.md`.
- Random seeds used for each (dataset, split, predictor) cell.

---

## Signed

Preregistration committed by Aayan Alwani via git timestamp. Any reviewer or co-author requiring verification can run `git log` on this file to confirm the commit predates Phase 1 first run.
