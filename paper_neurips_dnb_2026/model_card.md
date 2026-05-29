# Model Card for ConfPert Library

Following the Model Cards template (Mitchell et al. 2019, *FAT* '19).

**NOTE:** ConfPert is a wrapper + evaluation library, not a trained model. This model card describes the *library + benchmark* as the deliverable, not predictor weights. Each wrapped predictor (scGen, CPA, STATE, scGPT, etc.) ships with its own model card per its authors.

---

## Model details

- **Person/org developing model:** Aayan Alwani (Independent Researcher)
- **Library version:** `confpert` v0.2.0 (Phase 2 target release)
- **Type:** model-agnostic conformal-prediction wrapper for sample-producing predictors
- **Type of training:** none — ConfPert wraps already-trained predictors at inference / calibration time
- **License:** MIT (library code); pre-reg CC0; results CC-BY
- **Contact:** alwaniaayan6@gmail.com

## Intended use

### Primary intended uses
- Calibration evaluation of single-cell perturbation predictors at the per-population discrepancy level (KS, W1, energy, MMD, bimodality, variance-ratio)
- Pre-registered benchmarking under multiple-testing-corrected hypotheses
- Downstream PRISM-style enrichment analyses on calibrated bimodal subpopulation signatures
- Pre-registration template for ML benchmark protocols in general

### Primary intended users
- ML researchers working on single-cell predictor evaluation
- Computational biologists evaluating foundation cell model deployment readiness
- Methodologists studying pre-registration for ML benchmarks

### Out-of-scope use cases
- Clinical decision-making — ConfPert is not validated for patient-level decisions
- Drug repurposing — K3 PRISM signatures require independent biological validation before clinical use
- Predictor selection without license considerations — STATE has non-commercial restrictions

## Factors

### Relevant factors
The benchmark's findings are conditional on:
- Cell-line context (K562-derived vs. non-K562)
- Perturbation type (genetic CRISPRi/a vs. chemical small-molecule)
- Predictor architecture family (deterministic baselines, latent-delta VAE, transformer foundation)
- Training corpus diversity (number of distinct cell-line × perturbation contexts)
- Modality (RNA-only vs. RNA + protein for Perturb-CITE-seq)
- Organism (human vs. mouse)

### Evaluation factors
The K2 v2 hypothesis family explicitly tests interactions between these factors:
- H2: cell-line context
- H2b: corpus diversity
- H3: architecture family
With Bonferroni-corrected α = 0.0167 across the family.

## Metrics

### Performance measures
- Calibration deviation: `|target_coverage - achieved_coverage|`, lower is better
- BH-FDR-corrected pathway enrichment ratio (K3): higher is better
- Bootstrap 95% confidence intervals on every coverage value (n=1000 stratified resamples)

### Decision thresholds
Pre-registered in `preregistration_v2.yaml`:
- H1/H1b: ρ > 0.5, p < 0.05 (permutation), ≥3 of 4 main datasets
- H2: F-test p < 0.0167, η² > 0.10
- H2b: ρ < -0.5, p < 0.0167, ≥7 of 10 datasets
- H3: Kruskal-Wallis p < 0.0167, Cliff's δ > 0.30

### Variation approaches
Bootstrap CIs reported for all coverage values. Sensitivity to test-subsampling reported at n_eval_perts ∈ {15, 30, 60}. K2 v2 hypothesis tests are pre-registered before any new data is collected, eliminating analyst-side variability in threshold choice.

## Evaluation data

### Datasets
10 publicly available single-cell perturbation datasets (see Datasheet §3.1):
- Norman 2019 (K562)
- Replogle 2022 K562 + RPE1
- Adamson 2016 (K562)
- Tahoe-100M subset
- Frangieh 2021 (melanoma + Perturb-CITE-seq)
- Schmidt 2022 (primary T cells)
- Datlinger 2017 (Jurkat / primary)
- McFaline-Figueroa 2024 (chemical sci-Plex)
- Walker 2022 (mouse BMDC)

### Motivation
Datasets span: cell-line / primary cells, K562 / non-K562 contexts, genetic / chemical perturbation, human / mouse, RNA / RNA+protein. Selected for H2/H2b coverage.

### Preprocessing
Per pre-reg v2 §1.4: standard HVG selection (n_top_genes=512), log1p normalisation, filter cells (min_counts=200) and genes (min_cells=5).

## Training data

### N/A — wrapper, not a trained model
ConfPert does not train. Each wrapped predictor (scGen, CPA, biolord, sVAE+, GEARS, STATE, scGPT, scFoundation, Geneformer, CellFM, GET) is trained per its authors' default hyperparameters, pinned in `confpert.predictors.HPARAMS`.

## Quantitative analyses

### Unitary results (Phase 1, 2026-05)
- K1: 9 predictors × 5 datasets × 6 discrepancies × 3 α-levels = 2,028 calibration cells
- K2 v1 H1: 2 of 5 datasets pass (FAILS pre-reg overall threshold of ≥3 of 4)
- K3: 445 ConfPert signatures pass BH-FDR vs. 363 uncalibrated misses (1.23× ratio)

### Phase 2 target results
- K1 v2: ~12,960 calibration cells across expanded 14 × 10 grid
- K2 v2: H2/H2b/H3 hypothesis family with Bonferroni-corrected α; disposition per §2.6 of pre-reg v2
- K3 v2: K3 + K3a (DrugBank target match) + K3b (L2 ablation) + K3c (worked examples)

### Intersectional results
- Per-(cell-line-context × architecture-family) calibration heatmap
- Per-(perturbation-type × predictor-family) coverage scatter
- Pre-reg ablation: pre-registered thresholds vs. analyst-shopped thresholds (Tier-1 contribution if K2 v2 lands triple-null)

## Ethical considerations

### Sensitive use
Single-cell transcriptomic data; no patient identifiers. Cancer-cell-line dominance in the dataset list means findings may not generalise to primary-cell or healthy-tissue contexts. Schmidt 2022 (primary T cells) partially mitigates.

### Human life
None directly. K3 PRISM enrichment is hypothesis-generating for drug-target identification, not clinically actionable without orthogonal validation.

### Mitigations
- Reproducibility infrastructure (Modal image hashes, dataset SHA-256, random seed=42, semver, one-command rerun)
- Pre-registration prevents post-hoc threshold-shopping
- Honest reporting of failed hypotheses + null results

### Risks and harms
- Predictor mis-selection: if practitioners use a transformer predictor in a cell-line context where it fails calibration, downstream biological conclusions may be wrong. Mitigated by per-(predictor × dataset) reporting + explicit calibration deviation metrics.
- Inflated expectations of "foundation cell models": ConfPert's honest negative results may be misread as universal underperformance; we caveat per-context.

### Use cases (out of scope)
Same as Datasheet §5.3.

## Caveats and recommendations

- **Marginal, not conditional coverage:** ConfPert guarantees finite-sample marginal coverage; conditional coverage at every x is impossible distribution-free (Vovk 2012).
- **Cross-cell-line CP is undercoverage by design:** Tibshirani 2019 weighted CP is exposed but enabling per-domain weights is future work.
- **Cancer-corpus dominance:** Tahoe + most Perturb-seq datasets are cancer-cell-line. Frangieh + Schmidt partially mitigate. Conclusions may not transfer to primary cells without validation.
- **Wet-lab validation excluded:** per pre-reg v2 §6; deferred to Track C domain-venue paper.
- **Heavyweight-predictor restriction:** scGPT/scFoundation/Geneformer/STATE/CellFM run on 6 of 10 datasets due to Modal budget cap. Documented per cell.
