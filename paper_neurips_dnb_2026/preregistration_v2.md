# ConfPert K2 v2 Pre-Registration

**Project:** ConfPert Phase 2: Distribution-Aware, Pre-Registered Benchmarking of Single-Cell Perturbation Predictors.
**Author:** Aayan Alwani (Independent Researcher).
**Document created:** 2026-05-25.
**Supersedes:** none; this v2 pre-registration is *additive* to `../preregistration.md` (Phase 1, git commit `c7046e4b`, 2026-05-03 01:30 EDT). Phase 1 hypotheses (K1, K2 H1/H1b, K3 v1) are immutable. Phase 2 hypotheses (K2 v2 H2/H2b/H3 + intervention + K3 v2 sharpenings) are locked here before Phase 2 first run.
**Lock target:** committed via git + machine-readable YAML schema + OSF mirror before any Phase 2 new run (any wrapper invocation against any dataset not in the Phase 1 cell of `baselines/results.json`).
**Purpose:** lock K1 v2 expansion conditions, K2 v2 H2/H2b/H3 hypothesis family + outcome dispositions, K2 v2 H3 falsification-intervention design, K3 v2 sharpening criteria, pre-reg ablation panel contingency, multiple-testing-correction family, power-analysis floor, and the machine-verifiable YAML format.

This document supersedes any later informal hypothesis revisions for Phase 2 work. Post-hoc reframing of any criterion below is treated as exploratory and reported as such, not as a primary finding. The machine-readable companion `preregistration_v2.yaml` is the *authoritative* specification; this .md is the human-readable rendering of the same content.

---

## 0. Pre-registration discipline restated

The CLI verifier `confpert prereg verify` parses `preregistration_v2.yaml`, loads `baselines/results.json`, and runs each declared test against the actual data. The verifier refuses to run if `results.json` was modified after the YAML's git-commit timestamp + SHA-256 hash. This makes the pre-registration *machine-enforced*, not just rhetorical.

If a hypothesis below fails the pre-registered threshold, the disposition committed here is applied verbatim. No post-hoc reframing of a failed hypothesis as a pass at a relaxed threshold is permitted.

---

## 1. K1 v2: expanded benchmark conditions

### 1.1 Predictors (14 total)

Phase 1 (9, carried forward verbatim):
1. Mean (4 noise variants A/B/C/D)
2. Bilinear ridge
3. NoisyMean (control)
4. scGen
5. CPA
6. biolord
7. sVAE+ / SAMS-VAE
8. GEARS-uncertainty
9. STATE-SE-600M (ST-SE-Replogle + ST-SE-Tahoe)

Phase 2 additions (5):
10. **scGPT** — Cui 2024, transformer foundation model, ~100M params, bin-tokenised
11. **scFoundation** — Hao 2023 / xTrimo, transformer foundation, ~100M params
12. **Geneformer** — Theodoris 2023, transformer with rank-value encoding, ~10M params
13. **CellFM** — Zeng 2024, RetNet architecture, ~800M params
14. **GET** — Fu 2025 *Nature*, chromatin-conditioned, distinct input modality (ATAC + sequence)

Hyperparameters per predictor are pinned in `confpert.predictors.HPARAMS` table as of git commit at Phase 2A lock-in.

### 1.2 Datasets (10 total)

Phase 1 (5, carried forward verbatim):
1. Norman 2019 (K562 CRISPRa)
2. Replogle 2022 K562 essential
3. Replogle 2022 RPE1 essential
4. Adamson 2016 (K562 UPR)
5. Tahoe-100M cross-cell-line drug-perturbation subset

Phase 2 additions (5):
6. **Frangieh 2021** — Perturb-CITE-seq melanoma IFN response (multi-modal RNA + protein)
7. **Schmidt 2022** — Perturb-seq primary human T cells (primary-cell context)
8. **Datlinger 2017** — Jurkat T cells + primary cell lines (ECCITE-seq)
9. **McFaline-Figueroa 2024** — sci-Plex extended chemical perturbation
10. **Walker 2022** — Perturb-seq mouse BMDC (cross-organism)

Dataset SHA-256 fingerprints recorded at Phase 2C lock-in in `baselines/dataset_hashes_v2.json`. Loader code git-commit hash recorded in same file.

### 1.3 Splits (3 + 3 new)

Phase 1 (3, carried forward):
- (a) Within-perturbation 70/30 cell-level
- (b) Held-out-perturbation 50/50 calibration/test
- (c) Cross-cell-line K562 → RPE1

Phase 2 additions (3):
- (d) **Held-out organism**: mouse Walker 2022 → human Norman 2019 on gene-intersected space
- (e) **Held-out modality**: Frangieh 2021 RNA train → protein test
- (f) **Held-out perturbation-type**: Replogle 2022 (genetic) train → Tahoe (chemical) test on gene-symbol-overlapping perturbations

### 1.4 Discrepancies (unchanged — same 6 from Phase 1)

KS, Wasserstein-1, energy distance, MMD-RBF, bimodality coefficient match, variance ratio.

### 1.5 Coverage targets (unchanged)

α ∈ {0.05, 0.10, 0.20}.

### 1.6 Heavyweight-predictor dataset restriction

To stay within the $3,000 Modal budget (see `PHASE_2_PLAN.md` §14.5), the 5 heavyweight Phase 2 predictors (scGPT, scFoundation, Geneformer, STATE, CellFM) are run only on a 6-dataset subset: Norman, Replogle K562, Replogle RPE1, Tahoe, Frangieh, Schmidt. GET is demoted to optional extension and is not in the must-run set.

Lightweight + mid predictors (Mean, Bilinear, NoisyMean, scGen, CPA, biolord, sVAE+, GEARS) run on all 10 datasets where data-loading succeeds.

This restriction is **pre-registered**. Reporting MUST disclose missing (heavyweight × dataset) cells explicitly in K1 results table. The cell is reported as `N/A — heavyweight restriction per prereg v2 §1.6`.

---

## 2. K2 v2: hypothesis family

### 2.1 Multiple-testing family

The K2 v2 hypothesis family is {H2, H2b, H3}. Family-wise α controlled at 0.05 via Bonferroni correction → per-hypothesis α = 0.0167.

H1 and H1b from Phase 1 prereg are **separate** from the K2 v2 family. They are re-tested on the expanded grid (14 predictors × 10 datasets) for reporting purposes but use the original Phase 1 thresholds (α = 0.05, ≥3 of 4 main pre-reg datasets). H1/H1b are not part of the K2 v2 family-wise correction.

### 2.2 H2 — cell-line / data-corpus context covariate (cellular)

**Claim:** Calibration deviation differs significantly between K562-derived contexts and non-K562 contexts, controlling for predictor capacity.

**Operationalisation:**
- Dependent variable: `calibration_deviation = |target_coverage - achieved_coverage|`, averaged over (α ∈ {0.05, 0.10, 0.20}) × (6 discrepancies)
- Independent factors:
  - `cell_line_context`: binary, {K562-derived, non-K562}. K562-derived = {Norman, Replogle K562, Adamson}; non-K562 = {Replogle RPE1, Frangieh (melanoma), Schmidt (T cells), Datlinger (Jurkat partial), McFaline-Figueroa (multi), Walker (BMDC), Tahoe (multi)}. Borderline assignments documented in YAML.
  - `predictor_capacity_log`: continuous, log10(parameter_count), per predictor
- Test: two-way ANOVA with type-II sums of squares on the cell `(predictor × dataset)` × `cell_line_context × predictor_capacity_log`.
- Cell counts: 14 predictors × 10 datasets = 140 cells, minus heavyweight restrictions (5 × 4 = 20 missing) = ~120 cells.

**Success criterion (both must hold):**
- F-test p-value on the `cell_line_context` main effect < 0.0167 (Bonferroni-corrected)
- Effect size η² > 0.10 (medium-to-large by Cohen's convention)

**Failure to meet either condition** = H2 FAILS. No partial-pass framing.

### 2.3 H2b — data-corpus diversity covariate (training-data)

**Claim:** Calibration deviation decreases with training-corpus diversity, where diversity is the count of distinct (cell-line × perturbation) contexts in the predictor's training data.

**Operationalisation:**
- Compute `corpus_diversity_index = log10(N_distinct_train_contexts)` per predictor, declared in `confpert.predictors.HPARAMS.corpus_diversity_index`
- For each dataset, compute Spearman rank correlation between `corpus_diversity_index` and `calibration_deviation`
- Aggregate: number of datasets where ρ < -0.5 AND p < 0.0167 (permutation test, ≥10,000 permutations)

**Success criterion:** ρ < -0.5 (p < 0.0167) on ≥ 7 of 10 datasets.

**Failure to meet condition** = H2b FAILS.

### 2.4 H3 — architecture-family covariate

**Claim:** Architecture families have systematically different calibration profiles.

**Operationalisation:**
- Family assignment (locked here, see Phase 2 plan §14.12 + §3.1):
  - **F1: Deterministic + noise baselines** — {Mean(A), Mean(B), Mean(C), Mean(D), Bilinear, NoisyMean} (n=6)
  - **F2: Latent-delta / disentangled VAE** — {scGen, CPA, biolord, sVAE+} (n=4)
  - **F3: Transformer foundation** — {scGPT, scFoundation, Geneformer, STATE, CellFM} (n=5)
- Excluded from H3 with documented reason:
  - GEARS-uncertainty: GNN architecture, distinct from F1/F2/F3, n=1 → untestable as own family
  - GET: chromatin-conditioned distinct-input-modality, n=1 → untestable as own family
  Both still run in K1 and tested in H1/H1b/H2/H2b.
- Test: Kruskal-Wallis H-test across F1/F2/F3 means of `calibration_deviation`, averaged over datasets.

**Success criterion (both must hold):**
- Kruskal-Wallis p-value < 0.0167 (Bonferroni-corrected)
- Cliff's δ between best-family-mean and worst-family-mean > 0.50 (large effect by Cohen's convention)

**Aggregation:** unit is per-(predictor, dataset) cell (`mean_over_alphas_and_discrepancies`). Per-family n ≈ 40 (F1: 60, F2: 40, F3: 50) given the heavyweight-restriction exclusions.

**Failure to meet either condition** = H3 FAILS.

**Threshold rationale:** raised from δ > 0.30 (RC1) to δ > 0.50 (final lock) after `confpert.power.k2_v2_power_report()` revealed δ = 0.30 was substantially underpowered at the per-(predictor, dataset) aggregation sample size (0.36 power for n_per_group ≈ 40, α = 0.0167). At δ = 0.50 with the same sample size, simulated power is ≥ 0.80. Trade-off: large-effect threshold loses sensitivity to medium effects, but is the statistically-defensible choice at the data we will have. Documented in §6.

### 2.5 H3 falsification-intervention design (conditional on H3 result)

**Pre-committed design (run only if H3 PASSES, otherwise reported as "intervention skipped per prereg v2 §2.5"):**

```
IF best_family == F1 (baselines):
    intervention = "scale F2 latent-delta predictor (e.g., sVAE+) capacity 10x and retrain"
    target_predictor = "sVAE+ scaled to ~100M params"
ELIF best_family == F2 (latent-delta) AND worst_family == F3 (transformer):
    intervention = "add latent-delta encoder-decoder split head to a transformer baseline"
    target_predictor = "scGPT + latent-delta output head"
ELIF best_family == F3 (transformer) AND worst_family == F2 (latent-delta):
    intervention = "add a frozen pre-trained transformer encoder to a latent-delta predictor"
    target_predictor = "scGen + frozen scGPT encoder"
ELIF best_family == F3 AND worst_family == F1 (baselines):
    intervention = "augment a F1 baseline with transformer-derived gene embeddings as features"
    target_predictor = "Bilinear + scGPT gene embeddings"
ELSE:
    intervention = "skip — no clean architectural axis identified"
```

**Intervention success criterion:** retrained-modified predictor achieves calibration deviation reduction of ≥ 30% (relative) vs. the baseline-family member on ≥ 3 of 5 datasets in {Norman, Replogle K562, Replogle RPE1, Tahoe, Frangieh}.

**Intervention cost cap:** $300 Modal. If cap exceeded during intervention experiment, report partial results explicitly. Cap is enforced.

### 2.6 Outcome dispositions (2³ = 8 outcomes)

| H2 | H2b | H3 | Disposition | Headline framing |
|---|---|---|---|---|
| PASS | PASS | PASS | Triple-pass | "Three pre-registered covariates of calibration are confirmed: cell-line context, training-corpus diversity, and architecture family." |
| PASS | PASS | FAIL | Double-pass (no architecture story) | "Cell-line and corpus-diversity covariates pass; architecture-family does not." |
| PASS | FAIL | PASS | Double-pass (no corpus story) | "Cell-line and architecture-family covariates pass; corpus diversity does not." |
| FAIL | PASS | PASS | Double-pass (no cell-line story) | "Corpus diversity and architecture family pass; cell-line context does not." |
| PASS | FAIL | FAIL | H2-only | "Cell-line context dominates capacity; other axes are null." |
| FAIL | PASS | FAIL | H2b-only | "Training-corpus diversity dominates capacity; other axes are null." |
| FAIL | FAIL | PASS | H3-only | "Architecture family dominates capacity; other axes are null." |
| FAIL | FAIL | FAIL | Triple-null | "K2 v2 hypotheses do not survive pre-registered thresholds. Pre-registration-ablation-panel methodological contribution becomes the central paper finding (see §4)." |

**Triple-null contingency:** the pre-reg ablation panel (§4) is committed as Tier-1 contribution. Paper headline pivots to "Pre-registration prevents threshold-shopping bias in ML benchmarks" with the K2 v2 null result as the *worked example*.

### 2.7 Power analysis

Computed at lock-in time, before any Phase 2 run:

| Test | n_cells | Min detectable effect at 80% power, α=0.0167 |
|---|---|---|
| H2 ANOVA (cell_line main effect) | 120 | η² ≈ 0.07 (sufficient to detect medium effect 0.10) |
| H2b Spearman (per-dataset) | n=14 predictors (or 13 excluding NoisyMean control if specified) per dataset | min ρ ≈ ±0.62 — *borderline*; pre-reg threshold of ρ < -0.5 is below minimum detectable. Acknowledge: H2b is **underpowered** at per-dataset level; the success criterion requires ≥ 7 of 10 datasets to compensate, but each individual dataset test is underpowered. *Decision: accept this underpowering at lock-in.* Document in §6. |
| H3 Kruskal-Wallis (3 families, df=2) | per-(predictor, dataset) cells; ~40 per family (F1=60, F2=40, F3=50) | At δ = 0.30: 0.36 power (UNDERPOWERED). At δ = 0.50: 0.92 power (OK). Threshold locked at δ > 0.50 |

Power computations done with Python `statsmodels.stats.power` and verified with G*Power 3.1. Code snippet preserved in `paper_neurips_dnb_2026/power_analysis_v2.py` at lock-in time.

---

## 3. K3 v2: PRISM downstream sharpening

K3 Phase 1 result: 445 ConfPert signatures pass BH-FDR vs. 363 uncalibrated misses (ratio 1.23×). Sharpen with three additional layers:

### 3.1 K3a — biological interpretability layer

For each of the 445 ConfPert-passing signatures, query DrugBank + ChEMBL for known mechanism-of-action / target list. Compute the fraction of signatures whose top-5 Hallmark-pathway hits include a pathway containing the drug's known target.

**Pre-registered success criterion:** ≥ 30% of the 445 ConfPert-passing signatures match known target (a "strong" baseline; current literature reports ~20% match rate for unsupervised pathway-enrichment on PRISM, so 30% is a meaningful uplift).

### 3.2 K3b — calibration-vs-L2-rank ablation

Re-run K3 pipeline with L2-error-ranked signatures (top-50 genes by smallest L2 residual) instead of calibrated bimodal-mode signatures. Report BH-FDR pass count.

**Pre-registered success criterion:** ConfPert-calibrated signatures yield ≥ 1.30× the BH-FDR pass count of L2-ranked signatures. If ratio < 1.30 → K3 calibration value is weak; report honestly.

### 3.3 K3c — adversarial worked-example layer

Pick 3 specific compounds where ConfPert-calibrated signature passes BH-FDR AND L2-ranked signature does NOT. Manually inspect the signature gene list + pathway hits for biological plausibility (mechanism consistency with the compound's known activity).

**Pre-registered success criterion:** at least 2 of the 3 worked examples have biologically plausible mechanism explanations as judged by independent literature search. Report all 3 honestly regardless of judgement.

### 3.4 K3 v2 disposition

- K3 + K3a + K3b + K3c all pass → K3 promoted to second-headline finding
- Only K3 + K3b pass → K3 stays as third-headline, supporting role (status quo)
- Only K3 passes → K3 reported as in v1, downstream sharpenings reported as exploratory
- K3 fails on the expanded predictor set → K3 v1 result is *retracted* as exploratory and not load-bearing for v2 paper

---

## 4. Pre-registration ablation panel (Tier-1 if K2 v2 §2.6 triple-null)

**Committed at lock-in time as the Tier-1 contingency contribution if and only if K2 v2 lands triple-null.**

### 4.1 Design

Two parallel analysis pipelines run on the same Phase 2 results.json:

**Pipeline A — pre-registered:**
- Hypotheses, thresholds, and dispositions = exactly as specified in this document
- Run by `confpert prereg verify` against the locked YAML
- Output: PASS / FAIL per hypothesis with verbatim threshold

**Pipeline B — threshold-shopped (counterfactual):**
- Analyst (Aayan Alwani) sees results.json
- Picks threshold for each test that *maximises* the pass count for that hypothesis
- Documents the picked threshold explicitly (e.g., "ρ > 0.4 instead of 0.5", "p < 0.10 instead of 0.05", "≥ 4 of 10 datasets instead of ≥ 7")
- Output: PASS / FAIL per hypothesis with shopped threshold

### 4.2 Reporting

Side-by-side table:

| Hypothesis | Pre-reg threshold | Pre-reg result | Shopped threshold | Shopped result | Threshold-shopping Δ |
|---|---|---|---|---|---|
| H1 | ρ > 0.5, p < 0.05, ≥ 3 of 4 | (FAIL per v1 result, 2 of 5) | ρ > 0.4, p < 0.10, ≥ 3 of 5 | (likely PASS) | Δ = headline flip |
| H2 | F p < 0.0167, η² > 0.10 | ... | ... | ... | ... |
| H2b | ρ < -0.5, p < 0.0167, ≥ 7 of 10 | ... | ... | ... | ... |
| H3 | KW p < 0.0167, δ > 0.30 | ... | ... | ... | ... |

### 4.3 Argument

The threshold-shopping Δ is the *empirical demonstration* of pre-registration's methodological value. The paper reports it as the central methodological finding.

### 4.4 Disclosure

Pipeline B is run *after* Pipeline A and only for the ablation panel. The shopped thresholds and the analyst's reasoning are documented explicitly. No claim is made that Pipeline B "would have been chosen" — only that Pipeline B is a *plausible counterfactual* if pre-registration were absent.

---

## 5. K2 v2 verification protocol

The `confpert prereg verify` CLI executes the following on each invocation:

1. Parse `preregistration_v2.yaml`
2. Verify `sha256(preregistration_v2.yaml) == self.declared_sha256`
3. Verify `git log --oneline preregistration_v2.yaml | head -1` matches declared commit hash
4. Load `baselines/results.json`
5. Verify `results.json` modification timestamp > git commit timestamp of pre-reg v2 (i.e., results were generated AFTER pre-reg was locked)
6. For each declared hypothesis (H1, H1b, H2, H2b, H3):
   a. Filter results.json rows matching the hypothesis's predictor + dataset + α scope
   b. Compute the declared test statistic
   c. Compare to the declared threshold
   d. Emit PASS / FAIL with both observed value and threshold
7. Apply the §2.6 disposition table
8. Output the disposition headline
9. Hash the entire verification output + sign with git commit SHA

Output format: JSON + human-readable text. Cached as `paper_neurips_dnb_2026/prereg_v2_verification_<TIMESTAMP>.json`.

---

## 6. Acknowledged limitations of this pre-registration

- **H2b underpowered at per-dataset level** (see §2.7). The ≥ 7 of 10 rule partially compensates. Reviewers may flag this. Honest disclosure here.
- **Family assignment for H3 is researcher-chosen** at lock-in time, not unambiguously biology-derived. Alternative family-collapses (e.g., 4 families instead of 3) could yield different H3 outcomes. Documented as a robustness check in the §7 sensitivity panel.
- **Cell-line dichotomy for H2 (K562 vs non-K562) is researcher-chosen**. Alternative dichotomies (cancer vs primary, human vs mouse, etc.) tested as sensitivity checks but not pre-registered as primary.
- **Heavyweight × dataset restriction (§1.6)** means H1 / H2 / H2b cells are partially missing for the 5 Phase 2 heavyweight predictors. Reported as `N/A` with explicit count in K1 table.
- **Modal compute cap of $3,000** may force truncated runs. If cap is hit mid-experiment, the disposition is "K2 v2 results are partial; explicitly report which (predictor × dataset) cells were not run."

---

## 7. Sensitivity / robustness checks (committed but not part of primary hypotheses)

The following are **reported alongside primary K2 v2 results** but are not pre-registered as hypothesis tests:

1. H3 robustness under alternative family-collapses: {2 families (baselines vs deep-learning), 4 families (split F3 into transformer / large-transformer)}, {5 families (un-collapsed)}
2. H2 robustness under alternative cell-line dichotomies: {cancer vs primary, human vs mouse, K562 vs non-K562 (primary), CRISPRi vs CRISPRa vs chemical}
3. H1 / H1b re-run on expanded grid using the Phase 1 thresholds (≥ 3 of 4 main datasets, ρ > 0.5)
4. Calibration-deviation Pareto vs. mean L2 error (orthogonality claim)
5. Per-discrepancy hardness ranking
6. Bootstrap 95% CIs on every coverage number (n = 1000 stratified-by-perturbation resamples)
7. Sensitivity-to-test-subsampling: n_eval_perts ∈ {15, 30, 60}
8. Within-family vs between-family variance (sharper test of H3)

---

## 8. Hash record

The following must be recorded at Phase 2A lock-in (computed by `confpert prereg verify --emit-hashes`):

- Git commit hash of `preregistration_v2.md` at lock-in: `<TO_BE_FILLED_AT_LOCK_TIME>`
- Git commit hash of `preregistration_v2.yaml` at lock-in: `<TO_BE_FILLED_AT_LOCK_TIME>`
- SHA-256 of `preregistration_v2.md`: `<TO_BE_FILLED_AT_LOCK_TIME>`
- SHA-256 of `preregistration_v2.yaml`: `<TO_BE_FILLED_AT_LOCK_TIME>`
- OSF DOI: `<TO_BE_FILLED_AFTER_OSF_UPLOAD>`
- Modal image hashes (recorded per predictor wrapper run, stored in results.json row)
- Dataset SHA-256 fingerprints (recorded at Phase 2C lock-in)
- Random seed: 42 (unchanged from Phase 1)

---

## 9. Signed

Pre-registration v2 committed by Aayan Alwani via git timestamp + OSF mirror. Any reviewer or co-author requiring verification can:

1. Run `git log --follow preregistration_v2.{md,yaml}` to confirm the commit predates Phase 2 first run
2. Verify OSF DOI permanent timestamp predates Phase 2 first run
3. Run `confpert prereg verify --prereg preregistration_v2.yaml --results baselines/results.json` to reproduce the disposition computation
4. Diff against this `.md` file to confirm the YAML faithfully encodes the human-readable spec

Phase 2 lock-in target: 2026-06-08 (Phase 2A end per `PHASE_2_PLAN.md` §16).

---

## 10. Diff vs. Phase 1 preregistration.md

- Phase 1 K1: 8 predictors × 4 datasets × 3 splits. Phase 2 K1 v2: 14 predictors × 10 datasets × 6 splits.
- Phase 1 K2: H1 (capacity) + H1b (data). Phase 2 K2 v2: H1/H1b (carried, re-tested) + new H2/H2b/H3 family with Bonferroni-controlled α.
- Phase 1 had no falsification-intervention. Phase 2 adds conditional intervention design (§2.5).
- Phase 1 K3 had single criterion. Phase 2 K3 v2 adds K3a (target match) + K3b (L2 ablation) + K3c (worked examples).
- Phase 1 had no pre-reg ablation panel. Phase 2 adds it as Tier-1 contingency (§4).
- Phase 1 used markdown only. Phase 2 adds machine-verifiable YAML + CLI verifier.
- Phase 1 used git stamp. Phase 2 adds OSF mirror with permanent DOI.
