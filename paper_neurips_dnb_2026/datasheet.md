# Datasheet for ConfPert Benchmark

Following the Datasheets for Datasets template (Gebru et al. 2021, *Commun. ACM*).

**Benchmark name:** ConfPert v0.2
**Authors:** Aayan Alwani (Independent Researcher)
**Date:** Phase 2 draft 2026-05-25 (will be revised at Phase 2E)
**Companion paper:** "Capacity is not Calibration: A Pre-Registered Distribution-Aware Benchmark of Single-Cell Perturbation Predictors" (in prep, target ICLR 2027 D&B)

---

## 1. Motivation

### 1.1 Why was this benchmark created?
Single-cell perturbation prediction has crossed a saturation regime: deep learning foundation models (scGPT, scFoundation, scBERT, Geneformer, UCE) systematically underperform simple bilinear / mean baselines on standard mean-prediction metrics (Ahlmann-Eltze et al. 2025; Csendes et al. 2025; Wong et al. 2025; Bendidi et al. 2024). However, mean prediction is not the regime where these models will be used in practice — single-cell biology is intrinsically distributional. ConfPert provides the first pre-registered, distribution-aware benchmark with finite-sample conformal coverage guarantees on six per-population discrepancies for any sample-producing predictor.

### 1.2 Who created the benchmark, on whose behalf?
Created by Aayan Alwani (Independent Researcher). No funding source. No corporate sponsor.

### 1.3 Who funded the creation?
Self-funded. Modal compute spend ~$2,100 across Phase 1 + Phase 2 (target).

---

## 2. Composition

### 2.1 What do instances represent?
Each instance is a `(predictor, dataset, split, discrepancy, α-level)` cell. For each cell, ConfPert reports: target coverage (1 - α), achieved coverage from the conformal head, calibration deviation, conformal threshold τ, and per-cell metadata (Modal image hash, dataset SHA-256, predictor parameter count, runtime, random seed = 42).

### 2.2 How many instances total?
Phase 2 target: 14 predictors × 10 datasets × 6 splits × 6 discrepancies × 3 α-levels = ~15,120 cells, minus heavyweight-predictor × dataset restriction per pre-reg v2 §1.6 (~5 × 4 × 6 × 6 × 3 = 2,160 missing) → ~12,960 valid cells. Phase 1 had ~2,028 cells.

### 2.3 Is each instance a single data point?
Each instance is the *output* of one calibration + evaluation run. The *input* per instance is: (i) predictor configuration, (ii) dataset h5ad, (iii) split spec, (iv) discrepancy score function, (v) α target.

### 2.4 What data is included?
- Numerical coverage / deviation values (10-row × 14-column header in `results.json`)
- Pre-computed bootstrap CIs (Phase 2D output)
- Provenance metadata: random seed, Modal image hash, dataset SHA-256, git commit at run time

### 2.5 Is there a label or target?
The "label" for each cell is the *achieved coverage*. The "target" is `1 - α` (nominal coverage). The benchmark task is to minimise `|target - achieved|` (calibration deviation).

### 2.6 Is any information missing?
- Heavyweight predictors run only on 6 of 10 datasets (Phase 2 pre-reg §1.6) — 4 datasets per heavyweight predictor are documented N/A
- GEARS-uncertainty does not run on Tahoe (cell-gears does not natively support chemical perturbation)
- GET predictor is demoted to optional extension

### 2.7 Are relationships between instances made explicit?
Yes. The `results.json` row schema includes all (predictor, dataset, split, discrepancy, α-level) keys explicitly. Cells can be filtered to construct (predictor × dataset) × (discrepancy × α) tables.

### 2.8 Are there recommended splits?
The benchmark itself defines 6 splits (pre-reg v2 §1.3). Phase 2 K2 hypotheses operate on (predictor × dataset) aggregations across split-within-perturbation cells. Cross-cell-line and held-out-organism are reported separately.

### 2.9 Are there known errors / noise?
Phase 1 results.json has known systematic missingness for STATE × Tahoe (resolved in Phase 1 with pre-staged checkpoint workaround). All numerical coverage values include SHA-256 fingerprints to detect post-hoc modification.

### 2.10 Is the benchmark self-contained?
Mostly. The benchmark depends on:
- Source datasets (Norman 2019, Replogle 2022, etc.) — all publicly available, see §3
- Predictor weights — all publicly available except STATE (Arc Institute non-commercial license click-through)

### 2.11 Does the benchmark contain confidential / personal data?
No. All datasets are cell-line / primary-cell transcriptomic data without subject identifiers.

### 2.12 Could it be offensive / harmful?
No.

---

## 3. Collection process

### 3.1 How was data acquired?
Source datasets downloaded from public repositories:
- Norman 2019: GEO GSE133344
- Replogle 2022: GEO GSE164378 (K562 + RPE1)
- Adamson 2016: GEO GSE90546
- Tahoe-100M: cellxgene + Vevo (Phase 2 access required for full corpus)
- Frangieh 2021: GEO GSE166989
- Schmidt 2022: GEO GSE190604
- Datlinger 2017: GEO GSE92872
- McFaline-Figueroa 2024: GEO TBD
- Walker 2022: GEO GSE189574

### 3.2 What mechanisms were used to collect?
Each loader (`confpert.data.<dataset>`) downloads via standard HTTPS, computes SHA-256, applies pre-registered preprocessing (HVG selection, log1p, etc.), returns a `PerturbDataset`.

### 3.3 Is the benchmark a sample of a larger set?
No. All 10 datasets included in full (with HVG selection per pre-reg v2).

### 3.4 Who was involved in collection?
Aayan Alwani (sole). No human subjects.

### 3.5 Time frame?
Phase 1 collection: 2026-04 to 2026-05. Phase 2 collection: 2026-06 to 2026-08.

### 3.6 Ethical review processes?
N/A — public data, no subjects.

---

## 4. Preprocessing / labeling / cleaning

### 4.1 Was preprocessing done?
Yes. Pinned in pre-reg v2:
- `sc.pp.filter_cells(min_counts=200)`
- `sc.pp.filter_genes(min_cells=5)`
- `sc.pp.normalize_total(target_sum=1e4)`
- `sc.pp.log1p()`
- `sc.pp.highly_variable_genes(n_top_genes=512)`
- Subset to HVGs

### 4.2 Is raw data saved?
Yes. Raw h5ad files cached on Modal volume `causeflow-artifacts` at `/data/`.

### 4.3 Is preprocessing code public?
Yes. `confpert.data.<dataset>` modules. Git-stamped at Phase 2C lock-in.

---

## 5. Uses

### 5.1 Has the benchmark been used yet?
Phase 1 results (9 predictors × 5 datasets) used in the AI4Research workshop submission 2026-05-14.

### 5.2 What tasks could it support?
- Calibration evaluation of single-cell perturbation predictors
- K3-style downstream pathway-enrichment evaluation
- Pre-registration template for ML benchmarks beyond single-cell

### 5.3 Are there uses we should not support?
- Clinical decision-making (no validated link to patient outcomes)
- Drug repurposing without orthogonal validation
- Predictor selection without considering license restrictions (STATE = Arc non-commercial)

---

## 6. Distribution

### 6.1 Will the benchmark be distributed?
Yes. `pip install confpert` on PyPI v0.2.0 release. Source on GitHub (TBD URL).

### 6.2 How will it be distributed?
- PyPI: `confpert` package (library + CLI)
- GitHub: source + pre-reg + results.json + Modal entrypoints
- Codabench: competition page with auto-scoring (Phase 2D deliverable)
- Hugging Face Hub: predictor wrappers (under Apache-2.0)

### 6.3 When?
Phase 2 v0.2.0 release: target 2026-09-30 alongside ICLR D&B 2027 submission.

### 6.4 License?
- Code: MIT (per `pyproject.toml`)
- Pre-registration: CC0
- Results: CC-BY
- Datasets: per individual source licenses (mostly CC-BY)

### 6.5 IP / copyright restrictions?
- Tahoe-100M: Vevo clearance required for redistribution
- STATE checkpoints: Arc Institute non-commercial license

### 6.6 Export controls or regulatory restrictions?
None.

---

## 7. Maintenance

### 7.1 Who maintains?
Aayan Alwani at alwaniaayan6@gmail.com.

### 7.2 Updates?
- Pre-reg v2 frozen at Phase 2A.4 lock-in.
- Library: semver. v0.2.x for patches; v0.3.x for non-breaking additions.
- New predictors / datasets: PRs welcome via Codabench submissions + GitHub.

### 7.3 How will users be notified?
- `CHANGELOG.md` in repo
- PyPI release notes
- Codabench announcement bar

### 7.4 Plan to retire?
None planned. Snapshot at v0.2.0 → maintained as archival artefact even if active development pivots.
