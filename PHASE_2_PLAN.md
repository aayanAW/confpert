# ConfPert Phase 2 — Multi-Venue Pivot Plan

**Date:** 2026-05-25
**Status:** active — supersedes Phase 1 framing for downstream venues
**Author:** Aayan Alwani
**Predecessors:** `PHASE_1_PLAN.md` (initial benchmark scope), `preregistration.md` (locked K1/K2/K3 criteria, git-stamped 2026-05-03 commit `c7046e4b`).

---

## 0. Why this document exists

The AI4Research workshop reviewer rejected the v1 submission with three concerns: (i) off-topic for a testing-methodology workshop, (ii) ICML style file misused, (iii) prose reads in "LLM style." The independent senior-reviewer audit (2026-05-25) confirmed all three and added a fourth: the novelty claims around "no published conformal framework wraps standard population-level discrepancies" are contradicted by 2025-2026 multivariate / generative-model conformal-prediction literature (Meyer 2026, Zheng 2024, Thurin 2025, Ndiaye 2025, Dheur 2025, Yang 2026, Braun 2025).

The Phase 1 paper as currently framed cannot be developed into a serious top-tier ML methods publication without genuinely new statistical theory. However, the *empirical work, the pre-registration artefact, the library, and the K2 cell-line covariate finding* are real benchmark+protocol contributions that fit a Datasets-and-Benchmarks venue cleanly after scope expansion and reframing.

The AI4Research submission stays as-is (already submitted, awaiting review). Phase 2 builds a parallel, larger, properly-framed paper for archival venues without depending on AI4Research outcome.

---

## 1. Target venue cascade

Decision locked 2026-05-25: **primary target is ICLR 2027 Datasets & Benchmarks** (Sept-Oct 2026 deadline). Full-scope expansion, lower schedule risk, comparable archival venue prestige to NeurIPS D&B. NeurIPS D&B 2026 (June) is too tight for the full 14-predictor × 10-dataset scope at the depth this audit calls for.

Primary and backup venues:

| Slot | Venue | Type | Deadline (approx.) | Fit | Notes |
|---|---|---|---|---|---|
| **Primary** | **ICLR 2027 Datasets & Benchmarks** | Archival benchmark track | **~Sept-Oct 2026** | Strong | Full 14 predictors × 10 datasets. Time for K2 v2 hypotheses to be properly tested at scale. Time for proper LLM-style edit pass + external review |
| **Domain venue A** | Bioinformatics (Oxford) / Nature Methods Brief Communication | Single-cell methods | Rolling | Medium-strong if K3 PRISM is sharpened | Run in parallel after ICLR submission. Needs domain motivation lead, drop CP framing |
| **Domain venue B** | Cell Systems | Single-cell systems-biology | Rolling | Medium if K2 covariate finding is the lead | Needs at least one biological validation experiment |
| **Workshop** | NeurIPS 2026 AI4Science | Domain-aware ML workshop | ~Sept 2026 | Strong fit for current scope | Workshop, not archival; can run in parallel as a teaser for the ICLR D&B archival paper |
| **Skipped** | NeurIPS D&B 2026 | Archival benchmark track | ~June 2026 | Strong but too tight | Cannot hit full scope in 4 weeks; ICLR D&B 2027 is the cleaner shot |
| **In flight** | ICML 2026 AI4Research | Testing-methodology workshop | Submitted 2026-05-13/14 | Weak (per reviewer) | Decision May 31, 2026. Treat as throwaway |

Drop testing-methodology framing from Phase 2 entirely.

---

## 2. Reframed thesis

**Working title (v2):** *Capacity is not Calibration: A Pre-Registered Distribution-Aware Benchmark of Single-Cell Perturbation Predictors*

Alternates considered: "Calibration Decouples from Capacity in Single-Cell Perturbation Predictors" (conservative); "Conformal Coverage Reveals Why Foundation Cell Models Underperform" (aggressive). Lead with "Capacity is not Calibration" — short, claim-forward, memorable, falsifiable.

**One-paragraph pitch:**
Single-cell perturbation predictors are increasingly evaluated by mean-prediction metrics, but mean prediction is not the regime where these models will be used. We introduce a pre-registered, distribution-aware benchmark that wraps six per-population discrepancy scores in finite-sample conformal heads at four head levels and evaluates 14 published predictors across 10 publicly-available Perturb-seq and chemical-perturbation datasets, including held-out-cell-line, held-out-organism, and held-out-modality splits. Hypotheses, permutation thresholds, BH-FDR corrections, and outcome dispositions are committed to a machine-verifiable git-stamped pre-registration mirrored on OSF before any first run. The headline finding: capacity scaling does not improve calibration; cell-line and data-corpus context dominate, and architecture family predicts calibration profile more reliably than parameter count. As a falsification test, we modify the worst-family architecture along the best-family's distinguishing axis and show calibration improves, supporting the H3 claim beyond observational correlation. We release a pip-installable evaluation library (`confpert` v0.2 on PyPI), a Cell-Eval plug-in covering 5 of 6 discrepancies absent from the upstream Arc Institute registry, a Codabench leaderboard for community submissions, and a reusable machine-verifiable pre-registration template (`ml-benchmark-prereg-template`) with a CLI verifier.

**What is in:** larger benchmark (14×10), the K2 v2 H2/H2b/H3 pre-registered hypothesis family, the H3 *falsification-intervention experiment* (moves K2 v2 from observational to interventional — single most important methodological strengthener), the pre-reg-ablation panel (load-bearing if hypotheses land null), machine-verifiable pre-reg format with CLI verifier, library, Cell-Eval plug-in, Codabench leaderboard, datasheet + model card, K3 PRISM downstream, OSF pre-reg mirror.

**What is out:** "we introduce conformal coverage on six discrepancies" as a methods claim (engineering composition, not theory), testing-workshop framing, ICML style file, wet-lab validation (deferred to Track C).

**What is new vs. v1:**
- Scope: 14 predictors × 10 datasets vs. 9 × 5
- Frame: benchmark + pre-reg-as-methodology vs. CP-methods novelty
- H3 falsification-intervention experiment: predicts which architectural change should improve calibration based on K2 v2 finding, then runs it and measures
- Pre-reg-ablation panel made the central methodological contribution (load-bearing if K2 v2 lands null)
- Machine-verifiable pre-reg format with CLI verifier
- Codabench leaderboard
- Honest citation of 2025-2026 multivariate CP cluster (Meyer, Zheng, Thurin, Ndiaye, Dheur, Yang, Braun)
- Independent OSF pre-registration mirror

---

## 3. Scope expansion targets

### 3.1 Predictors (8 → 14)

Existing wrapped (Phase 1, all stay):
1. Mean (4 noise variants: A/B/C/D)
2. Bilinear ridge (Ahlmann-Eltze 2025)
3. NoisyMean (control)
4. scGen
5. CPA
6. biolord
7. sVAE+ / SAMS-VAE
8. GEARS-uncertainty
9. STATE SE-600M (ST-SE-Replogle + ST-SE-Tahoe)

To add (priority order, with per-predictor justification):
10. **scGPT** (Cui 2024) — most-cited foundation model. Required for comparison to existing benchmarks (Li 2025, Csendes 2025, Wong 2025). Defensive inclusion: omitting it invites reviewer complaint.
11. **scFoundation** (Hao 2023, xTrimo) — second most-cited. 100M params + 50M cells. Same defensive logic.
12. **Geneformer** (Theodoris 2023) — rank-value encoding, distinct from binning-based scGPT/scFoundation. Tests whether tokenization strategy correlates with calibration profile (H3 architecture-family).
13. **CellFM** (Zeng 2024) — 800M-param, RetNet architecture. Largest open scFM. Stress-tests H1 capacity hypothesis at 4× STATE.
14. **GET** (Fu 2025, *Nature*) — chromatin-conditioned, distinct input modality (ATAC-seq + sequence). Tests cross-modality calibration. Particularly novel because no other benchmark has evaluated GET on perturbation prediction with calibration metrics.

Family balance for K2 v2 H3 (must be locked before pre-reg v2):
- **Deterministic-baselines:** Mean (A noise), Bilinear ridge (n=2)
- **Noise-augmented baselines:** Mean (B/C/D noise variants), NoisyMean (n=4)
- **Sparse-additive:** sVAE+/SAMS-VAE (n=1) — note small family, may need to merge with latent-delta for testing
- **Latent-delta / disentangled:** scGen, CPA, biolord (n=3)
- **GNN-uncertainty:** GEARS-uncertainty (n=1) — small family; consider merging with latent-delta or excluding from H3
- **Transformer:** scGPT, scFoundation, Geneformer, STATE, CellFM (n=5)
- **Sequence-grounded:** GET (n=1) — small; merge with transformer or exclude

Decision before pre-reg v2 lock: collapse to 3 families with n ≥ 3 each — {deterministic+noise-baselines, latent-delta-VAE-family, transformer-family} — to give K2 v2 H3 reasonable statistical power. Document this collapse in pre-reg v2 BEFORE first run; document why excluded predictors (e.g., GEARS standalone, GET standalone) don't fit family assignments cleanly.

Defer (low marginal value or licensing friction in 4-month window): scBERT, UCE, scMamba, xVERSE, BMFM-RNA. Add as opt-in extensions after submission.

H3 falsification-intervention design: pick the strongest-family-distinguishing axis revealed by the H3 result (e.g., if sparse-additive is best, the distinguishing axis is "explicit per-mode latent partition"). Modify the worst-family architecture along that axis (e.g., add a sparse-additive head to a transformer baseline). Retrain. Show calibration improves. Pre-register this intervention design (which family wins → which intervention runs) BEFORE H3 results are computed. Budget: 1 retrained predictor variant on 5 datasets ≈ $300 Modal.

### 3.2 Datasets (5 → 10)

Existing (Phase 1, all stay):
1. Norman 2019 (K562 CRISPRa, ~7 cell lines tested in benchmark)
2. Replogle 2022 K562 essential
3. Replogle 2022 RPE1 essential
4. Adamson 2016 (K562 UPR)
5. Tahoe-100M subset (10 cell lines × 50 PRISM drugs)

To add (with per-dataset justification for H2 + H2b coverage):
6. **Frangieh 2021** (Perturb-CITE-seq, melanoma, IFN response) — multi-modal RNA + protein. Tests held-out-modality calibration. Cancer context (still K562-adjacent) but multi-modal extends the benchmark surface
7. **Schmidt 2022** (Perturb-seq, primary T cells) — *primary-cell context*. Breaks the cancer-corpus limitation. Critical for H2 cell-line covariate — if H2 lands, primary cells should behave most differently from K562
8. **Datlinger 2017 ECCITE-seq** (Jurkat T cells + primary cell lines) — earlier-era pioneering Perturb-seq. Useful as a historical baseline + small-scale sanity check
9. **McFaline-Figueroa 2024** (sci-Plex extended) — chemical perturbation, complements Tahoe. Larger compound space (~5,000 small molecules)
10. **Walker 2022** (Perturb-seq, BMDC) — non-K562 immune-cell context. Critical for H2 — BMDC is murine, also gives us cross-organism if we mirror with a mouse dataset

Reserve list (in case any of 6-10 has data-access friction): Tian 2021 (Perturb-seq human iPSC neurons), Roohani 2023 (the GEARS-paper original Norman/Adamson processing), Dixit 2016 (pioneering 4-perturbation), Jin 2020 (autism Perturb-seq).

Held-out splits to add (already enabled by the dataset choice):
- **Held-out organism** (mouse Walker 2022 → human Norman 2019, on the gene-intersected space) — addresses scoop-risk on cell-line-context finding by extending to organism shift
- **Held-out modality** (Frangieh 2021 RNA train → protein test) — tests whether calibration on RNA transfers to protein
- **Held-out perturbation type** (genetic train → chemical test, e.g., Replogle train → Tahoe test on gene-symbol-overlapping perturbation set) — tests perturbation-type generalisation

Data-access verification (must complete before Phase 2C starts):
- All 5 new datasets confirmed public + downloadable from GEO / cellxgene / Figshare / EBI
- All datasets have permissive licenses (CC-BY, MIT, or equivalent)
- All datasets have sufficient cell-counts post-QC (≥100 cells per perturbation for ≥30 perturbations) — exclude any that don't pass
- Two of the 5 above are flagged as needing license confirmation: Tahoe-100M (Vevo proprietary clearance status), Frangieh 2021 (CITE-seq protein layer rights). Verify before commit

### 3.3 Analyses to add

**Bootstrap CIs** on every coverage number (n=1000 resamples). Procedure: *stratified bootstrap over perturbations* — sample perturbations with replacement within calibration and test arms separately, recompute coverage on resampled set, take 2.5 / 97.5 percentile as 95% CI. Run on cached prediction sets to avoid re-running Modal. Cost: negligible (~$10).

**Sensitivity-to-test-subsampling**: vary `n_eval_perts` ∈ {15, 30, 60} and show coverage stability. Reports whether the headline claims survive sub-sampling. Cost: ~$50.

**Pre-reg ablation panel** (load-bearing if K2 v2 hypotheses land null — see Section 9):
- Run 1: *with* pre-reg disposition enforcement. All thresholds (ρ > 0.5, p < 0.05, ≥3-of-4-datasets, BH at q=0.05) locked from preregistration_v2.md.
- Run 2: *without* enforcement — analyst sees results, picks the threshold that maximises "pass count." Document explicitly.
- Report the gap: "without pre-registration, the headline would have been X (e.g., '4 of 5 datasets pass H1 at p < 0.10'); with pre-registration, the headline is Y (e.g., '2 of 5 datasets pass at p < 0.05'). The Δ is the threshold-shopping bias."
- This demonstrates the methodological contribution of pre-reg *empirically*, not just rhetorically. Centerpiece of paper Section 7.

**Power analysis** (must be in pre-reg v2 BEFORE first run): for each K2 v2 hypothesis, compute the minimum detectable effect at α = 0.0167 (Bonferroni-corrected) given:
- H1 / H1b / H2 / H2b Spearman tests: n = number of predictors per dataset (8–14). Min detectable ρ ≈ 0.65 at n = 10 with 80% power.
- H2 ANOVA: 140-cell (predictor × dataset) grid. Min detectable η² ≈ 0.05 with 80% power.
- H3 Kruskal-Wallis (3 families post-collapse): Min detectable Cliff's δ ≈ 0.25 with 80% power.
Document these *in pre-reg v2*. If pre-reg sample size is below the power-required size for the locked threshold → reduce the threshold or relax to "exploratory."

**Calibration-vs-accuracy Pareto**: scatter (mean L2 error) × (mean calibration deviation) for each (predictor, dataset). The argument: calibration is an orthogonal axis to accuracy. Some predictors dominate on calibration but not accuracy, and vice-versa. If no predictor dominates on both → confirms benchmark adds value beyond existing L2-only benchmarks.

**Per-discrepancy hardness ranking**: rank (KS, W1, energy, MMD, bimodality, variance-ratio) by how much they discriminate between predictors. Discriminating discrepancies are valuable; non-discriminating ones can be pruned in v3.

**Within-architecture-family variance**: for K2 v2 H3, report not just between-family means but within-family standard deviation. If within-family variance > between-family variance → H3 is null even if means look different. This is a sharper test than the headline F-statistic.

**H3 falsification intervention** (committed in pre-reg v2): take the worst-performing-family architecture, modify it along the best-performing-family's distinguishing axis, retrain on Norman + Replogle K562 + RPE1, show calibration deviation drops. The pre-reg specifies which family-pair the intervention will use *conditional on which families come out best/worst*, so the intervention is itself pre-registered and not data-snooped.

### 3.4 K2 v2 pre-registration
Decision locked 2026-05-25: **full H2 / H2b / H3 family**, pre-registered before any new run. The current K2 H1 ("capacity hurts calibration") fails at the pre-reg overall threshold (2 of 4 main datasets pass, need ≥3). The post-hoc cell-line covariate observation is honest but currently labelled "not pre-registered." A clean Phase 2 step: write a **K2 v2 pre-registration** that specifically tests the cell-line context, data-corpus, and architecture-family hypotheses on the expanded dataset list, git-stamp it before adding the new datasets. This converts the post-hoc observation into properly pre-registered hypotheses for the larger benchmark, and locks the architecture-family analysis as a primary finding rather than exploratory observation.

Full K2 v2 hypothesis family:
- **H2 (cell-line covariate):** Calibration deviation differs significantly between K562-derived and non-K562 contexts, controlling for predictor capacity. Test: two-way ANOVA-style F-test on the per-(predictor, dataset, cell-line-context) cell with α = 0.05 after Bonferroni correction across the 3 K2 v2 hypotheses (effective α = 0.0167). Success: F-test p < 0.0167 AND effect size η² > 0.10.
- **H2b (data-corpus covariate):** Calibration improves with corpus diversity (number of distinct cell-line × perturbation contexts in training). Test: Spearman ρ between log(corpus diversity index) and mean calibration deviation across predictors. Success: ρ < -0.5, p < 0.0167.
- **H3 (architecture-family):** Architecture families (Mean / Bilinear / sparse-additive {sVAE+} / latent-delta {scGen, CPA, biolord} / GNN-uncertainty {GEARS} / transformer {scGPT, scFoundation, Geneformer, STATE, CellFM}) have systematically different calibration profiles. Test: Kruskal-Wallis H-test across the 5 family means with α = 0.0167. Success: p < 0.0167 AND Cliff's delta > 0.30 between best-family and worst-family.

Outcome dispositions (must be locked in pre-reg v2 before first run):
- All 3 land → primary 3-hypothesis-confirmed story
- 2 of 3 land → publish with explicit "2 of 3" honest framing
- 1 of 3 lands → publish as single-hypothesis with the H1 honest-null carried forward
- 0 of 3 land → publish as null result, lean on K1 expansion + pre-reg ablation panel as the central contribution

Pre-register H2/H2b/H3 + dispositions BEFORE running the new datasets. Git-stamp the commit hash + SHA-256 of the new `preregistration_v2.md` before first new run. Mirror the pre-reg to a public OSF registration before the same first run for the strongest reviewer-defensible position.

---

## 4. Code / library / infrastructure changes

### 4.1 Library refactor (`confpert/src/confpert/`)
- Bump version 0.1 → 0.2 ; document API stability promise
- Move `cell_eval_plugin` to `confpert.celleval` with proper entry-point group registration
- Add `confpert.benchmark` module: declarative YAML/TOML config for (predictor × dataset × split) sweeps
- Add `confpert.report` module: generates the per-(predictor, dataset, score) figure + tables directly from `results.json`
- Add `confpert.power` module: power analysis for H1-style Spearman correlation tests
- Add bootstrap-CI utility under `confpert.metrics`

### 4.2 Predictor wrappers (`confpert/src/confpert/predictors/`)
Each new predictor gets its own module implementing the universal interface:
```
def predict_samples(X_ctrl_test, perturbation, n_cells) -> R^{n_cells x d}: ...
```
- `confpert/predictors/scgpt.py` — wrap scGPT via official API
- `confpert/predictors/scfoundation.py`
- `confpert/predictors/geneformer.py`
- `confpert/predictors/cellfm.py`
- `confpert/predictors/get.py`

Each wrapper must declare:
- Parameter count (for K2 H1)
- Training data size (for H1b)
- License terms (for redistribution)
- Heavyweight (requires Modal) vs lightweight (local) flag

### 4.3 Dataset loaders (`confpert/src/confpert/data/`)
Five new loaders, mirror the existing `replogle.py` / `norman.py` interface:
- `confpert/data/frangieh.py`
- `confpert/data/schmidt.py`
- `confpert/data/datlinger.py`
- `confpert/data/mcfaline_figueroa.py`
- `confpert/data/walker.py`

Each loader downloads from public GEO/cellxgene/Figshare, applies preprocessing pinned in `preregistration_v2.md`, returns a `PerturbDataset`.

### 4.4 Modal entrypoints (`confpert/scripts/modal_launch.py`)
- Add a new STATE-style pre-stage entrypoint per predictor where applicable
- Add `multi_calibrate(dataset_list, predictor_list)` mass-launch entrypoint
- Add a Modal-volume health-check function (verifies no stale partials)
- Reuse the v5-tested `.spawn()`-under-`--detach` pattern from this session

### 4.5 Test infrastructure
- Existing: 33/33 unit tests pass on local
- Add: integration test that runs the full K1 sweep on a synthetic 100-pert / 200-cell dataset in <60 s, asserts pre-reg thresholds compute correctly
- Add: CI workflow (GitHub Actions) that runs unit + integration tests on every PR

### 4.6 Compute estimate
Per `compute_estimate.md` v1 was ~$180 on Modal. v2 expansion estimate:

| Component | Predictors | Datasets | Modal $ |
|---|---|---|---|
| Existing 9 predictors × 5 new datasets | 9 | 5 | ~$300 |
| 5 new predictors × 10 datasets | 5 | 10 | ~$800 |
| Cross-cell-line + cross-organism splits | 14 | 4 splits | ~$200 |
| K2 v2 power analysis + ablations | n/a | n/a | ~$100 |
| Re-runs for bootstrap CI | all | all | ~$200 |
| Buffer (failures, debugging) | | | ~$300 |
| **Total estimated cumulative** | | | **~$1,900** |

Add to v1's $180 → ~$2,100 cumulative Modal spend for v2 paper.

---

## 5. Paper restructure (target: 9 pages NeurIPS D&B format)

### 5.1 Section plan
1. **Introduction** (~1 page) — lead with the K2 cell-line covariate finding, frame as distribution-aware calibration benchmark. Drop the "saturation" motivation. Drop "we introduce a model-agnostic conformal framework" framing.
2. **Related work** (~½ page) — honest acknowledgment of (a) multivariate / generative-model CP cluster (Meyer 2026, Zheng 2024, Thurin 2025, Ndiaye 2025, Dheur 2025), (b) parallel single-cell perturbation benchmarks (Li 2025, Bendidi 2024, Wong 2025), (c) ML pre-registration prior art (Hofman 2023).
3. **Benchmark design** (~1½ pages) — datasets, predictors, splits, discrepancies, pre-reg protocol. This is the contribution chapter.
4. **K1 calibration sweep** (~1 page) — main result table, calibration-vs-accuracy Pareto figure, per-discrepancy hardness ranking.
5. **K2 cell-line covariate analysis** (~1 page) — H1/H1b/H2/H2b/H3 results, the honest-null story, power analysis. THIS IS THE EMPIRICAL HEADLINE.
6. **K3 downstream PRISM utility** (~½ page) — calibrated signatures recover more BH-FDR pathways, robustness sweep.
7. **Pre-registration as transferable methodology** (~1 page) — the template, the disposition rules, the ablation panel showing how analyses would have differed without pre-reg. NEW SECTION not in v1.
8. **Library, Cell-Eval plug-in, reproducibility** (~½ page).
9. **Limitations and discussion** (~1 page).
10. **Datasheet / model card** (~1 page in appendix per NeurIPS D&B requirement).
11. **Appendix** (~10–15 pages, no page limit at D&B).

### 5.2 Figure plan (target: 4 main-text figures)
- **Fig 1:** Schematic of benchmark architecture (predictors × datasets × splits × discrepancies × heads × α levels) — replace any text-heavy "Method" prose
- **Fig 2:** K1 calibration deviation heatmap, all (predictor × dataset × discrepancy) cells with bootstrap CIs
- **Fig 3:** K2 cell-line covariate panel — per-dataset Spearman scatter, with K562/non-K562 colour, expanded to 10 datasets
- **Fig 4:** Calibration-vs-accuracy Pareto — argues calibration is orthogonal axis to L2/Pearson-delta

### 5.3 Table plan (target: 2 main-text tables)
- **Table 1:** Predictor characteristics (name, param count, training-cell count, license, citation, wrapped-in-ConfPert flag)
- **Table 2:** K1 mean calibration deviation across (predictor × discrepancy), averaged across datasets, with bootstrap 95% CIs

### 5.4 Datasheet & model card
NeurIPS D&B requires a datasheet for the benchmark (1 appendix section, follow Gebru et al. 2021 template) and a model card for the released library (Mitchell et al. 2019 template). Draft these from the existing `README.md` + `preregistration.md`.

---

## 6. Writing repair (the "LLM style" fix)

One pass with strict knife. Concrete rules:

- **Em-dash budget:** ≤ 2 per page. Replace excess with period + new sentence OR comma.
- **Drop rule-of-three lists:** "Findings: (i)… (ii)… (iii)…" → three separate sentences in prose.
- **Drop bold-shouting:** No `\textbf{PASS}` / `\textbf{FAIL}` / `\textbf{EXACT}` in prose. Reserve bold for table cells only.
- **Drop `\paragraph{}` in short sections:** No paragraph headers in any section under ½ page.
- **Replace inflated phrases:**
  - "load-bearing covariate" → "the main covariate" / "the strongest covariate"
  - "strictly better" → "lower" / "smaller deviation"
  - "explicitly not pre-registered" → "post-hoc"
  - "we find three things" → delete; restructure as prose
  - "X meets Y in our framework" → drop
- **Drop the "we introduce" repetition:** Use it ONCE in the abstract. Method section uses "ConfPert reports…" / "Each head calibrates…"
- **Drop self-praising adjectives:** "rigorous", "principled", "comprehensive", "extensive" — remove unless quantified.
- **Replace "saturation" motivation paragraph:** down to one sentence + multi-citation. Saturation is community consensus now (Ahlmann-Eltze 2025, Csendes 2025, Wong 2025, Bendidi 2024).
- **One declarative sentence per major claim:** no parenthetical asides for "what is" / "what is rare" constructions.
- **Read aloud test:** if a paragraph sounds like a textbook trying too hard, rewrite.

## 6.1 Style-file fix
- Delete the entire `icml2026.*` family from `paper_v2/` directory
- Pull official `neurips_data_2026.sty` + `neurips_data_2026.bst` when CFP opens
- If template not yet released, use `neurips_data_2025.sty` placeholder
- Strip the "Anonymous Authors¹" + "Title Suppressed" artifacts; the D&B template handles blind review natively
- Re-derive author block, affiliation block, abstract length for D&B requirements

---

## 7. Timeline & milestones

Locked decision: target ICLR D&B 2027 (Sept-Oct 2026 deadline). Roughly 4 months of runway. Full scope feasible.

### Track A — ICLR D&B 2027 (primary)
| Phase | Dates | Milestone |
|---|---|---|
| **Phase 2A: lock-in** | 2026-05-25 → 06-08 | Write `preregistration_v2.md` locking K2 v2 H2/H2b/H3 + dispositions. Git-stamp + SHA-256 + OSF mirror BEFORE any new run. Set up `confpert/paper_neurips_dnb_2026/` directory (rename later if needed). Choose NeurIPS / ICLR template once ICLR 2027 CFP opens |
| **Phase 2B: predictor expansion** | 06-09 → 07-06 | Wrap scGPT, scFoundation, Geneformer, CellFM, GET. Each wrapper passes integration test before Modal launch. Run all 14 predictors on existing 5 datasets. Bootstrap CI utility implemented. ~$700 Modal |
| **Phase 2C: dataset expansion** | 07-07 → 08-03 | Implement Frangieh, Schmidt, Datlinger, McFaline-Figueroa, Walker loaders. Run all 14 predictors on the 5 new datasets. Held-out organism + held-out modality splits. ~$900 Modal |
| **Phase 2D: K2 v2 analyses** | 08-04 → 08-24 | H2/H2b/H3 testing on full grid. Power analysis. Sensitivity-to-test-subsampling. Calibration-vs-accuracy Pareto. Per-discrepancy hardness ranking. Pre-reg ablation panel (with vs without disposition enforcement). ~$300 Modal |
| **Phase 2E: paper drafting** | 08-25 → 09-14 | Draft full 9-page paper + ~15-page appendix. Datasheet for the benchmark. Model card for the library. All four main figures generated from results.json. Reference 2025-2026 multivariate CP cluster honestly |
| **Phase 2F: editorial pass** | 09-15 → 09-28 | Mandatory LLM-style audit: em-dash budget, bold count, rule-of-three count, paragraph-header count, read-aloud test. External review by one domain expert. Final novelty re-check on arXiv / bioRxiv / OpenReview for Q3 2026 work that would scoop |
| **Submission** | ~09-30 | Submit to ICLR 2027 D&B (or earliest available archival D&B with open deadline) |

### Track B — NeurIPS AI4Science 2026 workshop (parallel teaser)
| Window | Goal |
|---|---|
| Aug 2026 | Compress Phase 2A-2C results into 4-page workshop submission to NeurIPS AI4Science. Frame as preview of full archival paper. Use existing K1 + partial K2 v2 |
| Sept 2026 | Submit AI4Science. Cite forthcoming ICLR D&B paper as full version |

### Track C — domain venue (Cell Systems / Bioinformatics / Nature Methods Brief)
Run after Track A submission. Lead with the cell-line covariate biology + K3 PRISM downstream, not the conformal methodology. Decision point: October 2026 after ICLR D&B submitted, depending on Track A scope and K3 strength.

### Track D — pre-reg template open-source (post-submission)
After ICLR D&B submission: extract `preregistration_v2.md` template into a standalone GitHub repo (`ml-benchmark-prereg-template`) + OSF mirror. Add a README, usage guide, and one applied example (ConfPert as the worked case). Aim for citable standalone artefact. ~½ day of work. Schedule: October 2026.

---

## 8. Risk register

| Risk | Likelihood | Impact | Mitigation |
|---|---|---|---|
| NeurIPS D&B June deadline impossible at full scope | High | Medium | Cut to 10 predictors × 7 datasets for Track A. Defer rest to Track B |
| scGPT / scFoundation wrapping harder than estimated | Medium | High | Allocate 2× buffer in W1-W2. Have CellFM as alternative if both fail |
| K2 v2 cell-line hypothesis fails on expanded data | Medium | Medium | Story still survives as honest-null; the pre-reg ablation panel becomes the headline instead. Don't post-hoc reframe |
| Scoop by Li 2025 follow-up or another benchmark | Medium-High | High | Pre-register publicly via OSF or arXiv preprint by 06-15. Stake the conformal-coverage + pre-reg angle |
| Conformal-prediction parallel work (Meyer 2026, Zheng 2024, etc.) overlaps further | Low-Medium | Low | Already partial scoop on methods novelty; v2 frames as benchmark contribution, so this is mitigated by reframing |
| Modal compute budget overrun | Low | Low | Hard-cap at $3,000; monitor weekly |
| Reviewer flags LLM-style prose again | Low if W4 editing is rigorous | High | Mandatory edit pass + external read. No `\textbf{}` outside tables |
| AI4Research acceptance forces awkward double-venue handling | Low (was 4/10 rated) | Low | If accepted: present workshop version, link to NeurIPS D&B v2 in workshop slides |

---

## 9. Decisions locked (2026-05-25)

| Decision | Choice | Rationale |
|---|---|---|
| Primary venue | **ICLR 2027 D&B** (Sept-Oct 2026) | NeurIPS D&B June 2026 too tight for full 14×10 scope. ICLR gives 4-month runway |
| K2 v2 scope | **Full H2 + H2b + H3 family** | "Best possible paper" — locks the architecture-family + cell-line + data-corpus story as primary pre-registered findings rather than post-hoc observations |
| Modal budget | **$3,000 hard cap** | Sufficient for full 14×10 scope + bootstrap + buffer per Section 4.6 |
| Pre-reg template open-source | **Yes, post-submission** (Track D) | Extract `preregistration_v2.md` into standalone GitHub + OSF mirror after ICLR D&B submission. Adds reviewer-pleasing transferable-methodology angle without blocking the submission |
| Wet-lab validation | **Out of scope for v2** | Per `preregistration.md` original exclusion. Re-evaluate for Track C domain submission |
| AI4Research v1 submission | **Leave as-is, await result** | Already submitted. Not load-bearing for Phase 2 |

Remaining open items not blocking Phase 2A start:
- External reviewer / co-author for the final editorial pass (Phase 2F) — adding one materially reduces LLM-style risk. To be confirmed by ~September 2026.
- Domain venue parallel track (C) — decision in October after Track A submission.

---

## 10. What stays from Phase 1

The following Phase 1 artifacts carry forward to Phase 2 unchanged:

- `preregistration.md` (git-stamped 2026-05-03 commit `c7046e4b`) — historical record of K1/K2/K3 v1 lock-in. Cited verbatim. Do NOT modify.
- `baselines/results.json` (SHA-256 fingerprinted rows) — every Phase 1 row stays, Phase 2 appends new rows
- The library API surface (`SCORES`, `PerturbationConformal`, `CDSplitConformal`, `CauchoisSubgroupConformal`, `RCPSCalibrator`) — append-only changes, no breaking renames
- The 5 Phase 1 datasets — all stay in the Phase 2 benchmark
- The 9 wrapped predictors — all stay
- The Cell-Eval plug-in — extends, doesn't reset
- The Modal app `confpert` + `causeflow-artifacts` volume — reused

The Phase 1 paper sources under `confpert/paper_ai4science/` stay frozen for the workshop submission. Phase 2 paper lives in a new directory `confpert/paper_neurips_dnb_2026/` with its own `main.tex`, `refs.bib`, `figures/`, and template.

---

## 11. Verification / self-check protocol

Per the project research-grade standards (novelty verification protocol), every Phase 2 milestone gets self-check:

- After K2 v2 pre-reg drafted: cite-verify any new citations against primary sources via Consensus / Semantic Scholar before commit
- After each new predictor wrapper: integration test on synthetic dataset before Modal launch
- After each new dataset loader: schema validation + cell-count sanity check
- After paper v2 draft: independent novelty verification — re-search arXiv / bioRxiv / OpenReview for any 2026-Q2 / Q3 work that would scoop the K2 v2 angle
- Before submission: full LLM-style audit pass (em-dash count, bold count, rule-of-three count, paragraph-header count); read aloud test
- Before submission: compile-check (0 unresolved `[?]` markers), anonymization-check, page-limit-check

---

## 12. References (verified via Consensus 2026-05-25)

Citations used in this plan, all verified against primary sources during the senior-reviewer audit. Full bibliographic detail in the cited audit document.

- Ahlmann-Eltze et al. 2025, *Nature Methods* — foundation model saturation
- Csendes et al. 2025, *BMC Genomics* — Train-Mean saturation
- Wong et al. 2025, *Bioinformatics* — Pfizer perturbation benchmark
- Bendidi et al. 2024 — perturbation foundation model benchmark
- Li et al. 2025 (bioRxiv) — 9 models × 17 datasets perturbation benchmark
- Hofman et al. 2023 — pre-registration for predictive modelling
- Nosek et al. 2018, *PNAS* — pre-registration foundations
- Pineau et al. 2020, *JMLR* — NeurIPS 2019 reproducibility report
- Klein et al. 2025 (otcp2025) — OT-based multivariate CP
- Thurin et al. 2025 — OT-based CP (parallel to Klein)
- Ndiaye 2025 — multivariate CPDs
- Dheur et al. 2025 — generalized multi-output CP
- Meyer et al. 2026 — kernel nonconformity (MMD) CP
- Zheng et al. 2024 — generative CP with vectorized scores
- Yang et al. 2026 — CP4Gen for generative models
- Braun et al. 2025 — min-volume conformal sets
- Asiaee et al. 2026 — selective CP for perturbation/intervention data
- Tibshirani et al. 2019 — weighted CP under covariate shift
- Cauchois et al. 2021 — subgroup-conditional conformal sets
- Bates et al. 2021 — RCPS
- Romano et al. 2019 — CQR
- Izbicki et al. 2022 — CD-split / HPD-split
- Vovk 2012 — distribution-free impossibility (conditional)
- Wu et al. 2025, *Genome Biology* — scFM benchmark
- Xu et al. 2021, *Patterns* — Codabench
- Gebru et al. 2021 — Datasheets for Datasets
- Mitchell et al. 2019 — Model Cards

---

## 13. Sign-off

This plan supersedes Phase 1's workshop-paper framing for archival venues. Phase 1 deliverables remain valid for their original target (AI4Research workshop) and as the foundation for Phase 2. Phase 2 commits to ICLR 2027 D&B (or earliest open archival D&B with Sept-Oct 2026 deadline) as the primary venue. The pre-registered scope expansion, the K2 v2 H2/H2b/H3 hypothesis family, the H3 falsification-intervention, the pre-reg-ablation panel as load-bearing methodological contribution, and the open-source pre-registration template are the five core new commitments.

Next concrete action upon user approval: open `confpert/paper_neurips_dnb_2026/` directory, draft `preregistration_v2.md`, git-stamp it.

---

## 14. Refinements (2026-05-25 ultrathink pass)

The v1.0 plan was an outline. The v1.1 refinements below close 9 gaps identified in deep review. They are normative additions, not optional. Read them as load-bearing extensions to Sections 1-13.

### 14.1 Venue cascade — honest uncertainty about ICLR D&B
ICLR's Datasets & Benchmarks track exists as of ICLR 2024 onward but is less established than NeurIPS D&B. Confirm exact deadline + CFP terms when ICLR 2027 CFP opens (typically June-July 2026). Backup primary if ICLR D&B 2027 is discontinued or has incompatible scope: **NeurIPS D&B 2027** (~June 2027 deadline, gives 12 months runway — even safer). Secondary backup: **ICML D&B 2027** if it exists (~Jan-Feb 2027 deadline). In all three cases, the same Phase 2 work product applies; only the template + deadline shift.

Action: monitor ICLR CFP weekly starting July 2026. Have the paper draft polish-ready by Aug 25 so it can be retargeted within 2 weeks of any venue change.

### 14.2 Scoop-defense strategy
Independent researcher targeting a hot subfield needs explicit defensive moves:

1. **Public pre-registration via OSF** at Phase 2A end (2026-06-08). OSF assigns a permanent timestamp + DOI; this stakes priority on the K2 v2 hypotheses.
2. **Pre-arXiv pre-print at Phase 2D end** (~2026-08-24). Even if paper not done, post a 4-page "extended abstract" with the K2 v2 headline finding to stake priority on the empirical result.
3. **Weekly literature scan** during Phase 2: Mon mornings, search arXiv / bioRxiv / OpenReview for any work that would scoop the cell-line covariate finding, the multivariate-CP-for-perturbation framing, or pre-registration-for-ML-benchmarks. Maintain `confpert/lit_notes/scoop_watch.md` log.
4. **Public benchmark URL** via Codabench (Phase 2D), so the artefact is community-citable before the paper is reviewed.
5. **Twitter / X / Bluesky / Mastodon thread** at OSF pre-reg post + at arXiv preprint — short technical summary with figures. Tag the relevant communities (#singlecell, #conformalprediction, #benchmarks).
6. **Email outreach** to specific PIs at Phase 2D: Aviv Regev (Genentech), Mohammad Lotfollahi (Theislab/Sanger), Tom Norman (MSKCC), Joseph Replogle (Weissman lab), Anastasios Angelopoulos (Berkeley CP), Hannah Bost (Theislab). Offer pre-print + benchmark URL; ask for feedback. The point isn't co-authorship; it's making sure they're aware before review.
7. **Submit a talk** to a relevant workshop (NeurIPS AI4Science, NeurIPS ML4H, RECOMB 2026 satellite). A talk is a public record that pre-dates the paper.

### 14.3 External engagement / collaborator outreach
Independent-researcher route is harder. Engagement strategy:

| Target | Why | Ask |
|---|---|---|
| Arc Institute (STATE authors: Adduri, Patapoutian, Davis et al.) | They built Cell-Eval; the plug-in is for their backbone | Share Cell-Eval plug-in pre-release. Ask if they want to merge it upstream. Possible co-author or acknowledgment |
| Lotfollahi / Theislab | Built scGen, biolord, CPA — direct stakeholder in calibration findings | Share K2 v2 results pre-submission. Ask about interpretation of architecture-family findings |
| Replogle / Norman (data authors) | Their datasets are central; calibration findings affect interpretation of their data | Notify of result, ask for biological interpretation review |
| Angelopoulos / Bates (CP) | Methodological expertise; can validate "off-the-shelf conformal" framing | Share pre-print, ask for review of CP-method positioning |
| Hofman (ML pre-reg) | Direct prior art on ML pre-reg | Share pre-reg template, ask for endorsement or co-authorship on a future pre-reg-template paper |

Action: schedule outreach emails by 2026-08-15 (Phase 2D mid-point) so responses come back before 09-15 editorial pass.

### 14.4 Reproducibility infrastructure
Concrete commitments:

- **Random seeds:** seed 42 (matches Phase 1). Document in every results.json row.
- **Modal image hashes:** record each Modal image hash in results.json row.
- **Dataset SHA-256:** every input h5ad gets SHA-256 fingerprint computed on first load, stored in results.json row.
- **Container hash:** Phase 2 paper repo built with Docker; image hash recorded per submission.
- **One-command reproduction:** `python -m confpert.reproduce --paper paper_neurips_dnb_2026 --predictor all --dataset all` runs the full Phase 2 benchmark from scratch on Modal. Cost-estimate flag prints expected $ before launch.
- **Semantic versioning:** `confpert` v0.2.0 = Phase 2 release. v0.2.x for patches. v0.3.x for non-breaking additions. Deprecation warnings for renamed methods over 3-version horizon.
- **Public changelog:** `confpert/CHANGELOG.md` per Keep-A-Changelog format.
- **Pinned dependency versions:** `requirements.txt` with exact `==` pins; `pyproject.toml` with version ranges for downstream installs.

### 14.5 Per-predictor compute budget caps
Hard caps prevent any single predictor from eating the budget:

| Predictor | Per-(predictor, dataset) cap | Total cap across 10 datasets | Rationale |
|---|---|---|---|
| Mean / Bilinear / NoisyMean | $5 | $50 | Local-compute mostly free |
| scGen | $20 | $200 | A100, ~15 min/dataset |
| CPA | $30 | $300 | A100, ~20 min/dataset |
| biolord | $20 | $200 | A100, ~15 min/dataset |
| sVAE+ / SAMS-VAE | $20 | $200 | A100, ~15 min/dataset |
| GEARS-uncertainty | $50 | $500 | Custom PertData prep is slow |
| STATE | $50 | $500 | Inference fast but pre-stage download can recur |
| scGPT | $50 | $500 | Fine-tune step on each dataset |
| scFoundation | $50 | $500 | Same |
| Geneformer | $40 | $400 | Smaller than scGPT |
| CellFM | $80 | $800 | 800M-param, larger |
| GET | $40 | $400 | Chromatin-conditioned, sequence input |
| **Subtotal** | | **$4,550** | |

Subtotal exceeds the $3,000 cap. Therefore: drop CellFM or GET from the must-have list, demote to optional extension. Phase 2 v1.1 decision: **drop GET from Phase 2A-2F**, keep CellFM. New subtotal: $4,150. Still over. Additional decision: **run heavyweight predictors (scGPT, scFoundation, Geneformer, STATE, CellFM) only on 6 of 10 datasets** (Norman, Replogle K562, Replogle RPE1, Tahoe, Frangieh, Schmidt — the 6 most important). New subtotal: $3,200. Acceptable with $3K cap + small overage tolerance.

If any predictor hits its per-predictor cap: stop, document partial-results-only status, do not exceed cap on retries.

### 14.6 Machine-verifiable pre-registration format (the methodological artefact)
The pre-reg template's novel contribution is being *machine-verifiable*. Concrete spec:

`preregistration_v2.md` is the human-readable spec. Companion file `preregistration_v2.yaml` is the machine-readable schema:

```yaml
version: "1.0"
project: confpert
git_commit: <hash>
sha256_of_yaml: <hash>
osf_doi: <doi>
hypotheses:
  H2:
    statement: "Calibration deviation differs..."
    test: "two_way_anova"
    metric: "calibration_deviation"
    factors: ["predictor_family", "cell_line_context"]
    alpha: 0.0167
    bonferroni_family: ["H2", "H2b", "H3"]
    success_criteria:
      - test: "f_test_pvalue"
        operator: "<"
        threshold: 0.0167
      - test: "eta_squared"
        operator: ">"
        threshold: 0.10
    disposition:
      pass: "primary_finding_cell_line_covariate"
      fail: "null_result_publish_honestly"
```

The CLI verifier `confpert prereg verify` (a) parses the YAML, (b) loads `results.json`, (c) runs each declared test against the actual data, (d) reports PASS/FAIL with the pre-registered criteria, (e) refuses to be run if `results.json` was modified after the YAML's `git_commit` timestamp. Outputs a signed report.

This is the genuine methods-contribution thread that survives independent of K2 v2 outcome. Could be its own paper at NeurIPS Workshop on Reproducibility in ML.

### 14.7 Codabench leaderboard
Public leaderboard at `codabench.org/competitions/confpert-2026`:
- Datasets: all 10 Phase 2 datasets, pre-split per pre-reg
- Metrics: 6 discrepancies × 4 heads × 3 α-levels (configurable)
- Submission: predictor implements the `predict_samples` interface, packaged as Docker image
- Auto-scoring runs the conformal pipeline + reports calibration deviation
- Public leaderboard with bootstrap CI bars

Action: open Codabench competition by Phase 2D end (~2026-08-24). Estimated effort: 3-5 days work to package + submit.

### 14.8 Datasheet & Model Card concrete spec
Per NeurIPS D&B requirements:

**Datasheet for the benchmark** (~3 pages appendix, per Gebru 2021):
1. Motivation (why this benchmark exists)
2. Composition (10 datasets, 14 predictors, 6 discrepancies)
3. Collection process (per-dataset provenance, downloads, preprocessing)
4. Preprocessing/labeling (HVG selection, normalisation, log1p, etc.)
5. Uses (calibration evaluation, downstream PRISM analysis, predictor development)
6. Distribution (PyPI + Codabench + GitHub)
7. Maintenance (changelog, deprecation policy, contact)

**Model card for `confpert` library** (~1 page appendix, per Mitchell 2019):
1. Model details: pip-installable Python library, v0.2.0
2. Intended use: calibration evaluation of single-cell perturbation predictors
3. Factors: cell-line context, perturbation type, model architecture family
4. Metrics: 6 discrepancies × 4 heads × 3 α-levels
5. Evaluation data: 10 datasets per Phase 2 manifest
6. Training data: N/A — wrapper, not a trained model
7. Ethical considerations: cancer-corpus dominance, primary-cell scarcity, no clinical decision-making use
8. Caveats: marginal-not-conditional coverage, exchangeability assumption

### 14.9 Contingency: K2 v2 lands null
If H2 + H2b + H3 all fail their pre-registered thresholds, the **pre-reg-ablation panel becomes the central paper contribution**. Restructure paper Section 7 (currently planned as a sub-section) to be Section 4 — the primary methodological contribution. Headline becomes:

*"Pre-registration protocols prevent threshold-shopping bias in ML benchmarks: we demonstrate that without enforcement, the K2 H1 capacity hypothesis would have been reported as a positive finding (4 of 5 datasets pass at threshold-shopped p < 0.10); with enforcement, the headline is the honest null. A machine-verifiable pre-registration format with a CLI verifier turns this enforcement into an auditable artefact."*

The benchmark scope expansion becomes the supporting empirical contribution. This pivot is pre-committed in this plan so it isn't a post-hoc save.

### 14.10 K3 PRISM downstream sharpening
Current K3 result: 445 ConfPert signatures pass BH-FDR vs. 363 uncalibrated misses (1.23× ratio). Modest. For Phase 2:

- Add a *biological interpretability* layer: of the 445 ConfPert-pass signatures, how many match known drug-target literature (DrugBank, ChEMBL)? If ≥50% match → much stronger result.
- Add a *gene-set-enrichment ablation*: would the same 445 signatures pass if we used L2-error-ranked signatures instead of calibrated signatures? If only ConfPert-calibrated signatures pass → calibration is causally important.
- Add an *adversarial* layer: pick 3 compounds where ConfPert-calibrated signature passes but uncalibrated misses; manually inspect for plausibility. Document worked examples.

If all three sharpenings land → K3 becomes a strong second-headline finding. If only one or two → K3 stays as third-headline, supporting role.

### 14.11 H3 falsification-intervention design (committed in pre-reg v2)
The single most important methodological strengthener.

Pre-register the intervention design *conditional* on H3 outcome (not data-snooped):

```
IF best_family == "sparse_additive" AND worst_family == "transformer":
    intervention = "add sparse-additive prior head to a transformer baseline"
    target_predictor = "Geneformer + sparse-additive output head"
ELIF best_family == "latent_delta" AND worst_family == "transformer":
    intervention = "add latent-delta encoder-decoder split to a transformer baseline"
    target_predictor = "scGPT + latent-delta head"
ELIF best_family == "transformer" AND worst_family == "sparse_additive":
    intervention = "scale sparse-additive predictor to transformer-comparable capacity"
    target_predictor = "sVAE+ scaled to 100M params"
ELSE:
    intervention = "skip — H3 ordering doesn't suggest a clear architectural axis"
    report_in_paper_as: "intervention skipped per pre-reg L<X>"
```

Run intervention on 5 datasets. Compare calibration deviation: intervention-modified vs. baseline. Pre-registered success: calibration deviation drops by ≥ 30% on ≥ 3 of 5 datasets.

Cost: ~$300 Modal (1 retrained predictor on 5 datasets).

### 14.12 Family-balance rebalancing for K2 v2 H3
Per refinement in Section 3.1: collapse to 3 families with n ≥ 3 to enable statistical testing. Decision locked here:

| Family | Members | Count |
|---|---|---|
| **Deterministic + noise baselines** | Mean(A), Mean(B), Mean(C), Mean(D), Bilinear, NoisyMean | 6 |
| **Latent-delta / disentangled VAE** | scGen, CPA, biolord, sVAE+/SAMS-VAE | 4 |
| **Transformer foundation** | scGPT, scFoundation, Geneformer, STATE, CellFM | 5 |

Excluded from H3 (clean reasons documented in pre-reg v2):
- GEARS-uncertainty: GNN architecture, distinct from all 3 families; n=1 makes it untestable
- GET: chromatin-conditioned, distinct input modality (sequence + ATAC, not RNA only); n=1
- Both run in K1 + K2 H1 + H2 + H2b; only excluded from H3

This locks total H3 sample to 15 predictors-with-noise-variants × 10 datasets = 150 cells. With 3 families, Kruskal-Wallis power is reasonable at 80% for medium effect.

---

## 15. Updated decision log (post-refinement)

| Item | Decision | Status |
|---|---|---|
| Primary venue | ICLR 2027 D&B with NeurIPS 2027 D&B as automatic fallback if ICLR D&B discontinued | Locked |
| K2 v2 hypothesis family | Full H2 + H2b + H3 | Locked |
| H3 falsification intervention | Pre-registered conditional design (Section 14.11) | Locked |
| H3 family collapse | 3 families {baselines, latent-delta, transformer}; GEARS + GET excluded from H3 with documented reason | Locked |
| Modal budget | $3,000 hard cap | Locked, but heavy predictors restricted to 6 of 10 datasets (Section 14.5) |
| GET predictor | Demoted to optional extension; not in Phase 2 must-have | Locked |
| Pre-reg ablation panel | Tier-1 contribution; load-bearing if K2 v2 lands null | Locked |
| Codabench leaderboard | Required deliverable by Phase 2D end | Locked |
| Datasheet + Model Card | Required appendix deliverables | Locked |
| Machine-verifiable pre-reg YAML + CLI verifier | Required artefact | Locked |
| OSF pre-reg mirror | Required at Phase 2A end | Locked |
| arXiv pre-print | Required at Phase 2D end | Locked |
| Pre-reg template open-source | Post-submission (Track D, Oct 2026) | Locked |

---

## 16. Updated next-action queue (start order)

### Phase 2A — DONE (2026-05-25)

| Milestone | Deliverable | Path |
|---|---|---|
| 2A.1 | Paper directory + README + structure | `paper_neurips_dnb_2026/{README.md,figures/,tex/,bib/}` |
| 2A.2 | `preregistration_v2.md` (10 sections) | `paper_neurips_dnb_2026/preregistration_v2.md` |
| 2A.3 | `preregistration_v2.yaml` + `confpert.prereg` module + CLI | `paper_neurips_dnb_2026/preregistration_v2.yaml`, `src/confpert/prereg.py`, `src/confpert/cli.py` |
| 2A.4 | SHA-256 lock + OSF upload manifest | `paper_neurips_dnb_2026/preregistration_v2.sha256`, `osf_upload_manifest.md` |
| 2A.5 | Dry-run verify passes against Phase 1 `results.json` (Lock check passed: True) | `paper_neurips_dnb_2026/prereg_v2_dryrun_20260525.json` |
| Bonus | Datasheet (Gebru 2021) + Model Card (Mitchell 2019) | `paper_neurips_dnb_2026/datasheet.md`, `model_card.md` |
| Bonus | Phase 2 predictor metadata stubs + Phase 2 dataset metadata stubs | `src/confpert/predictors_v2_stubs.py`, `src/confpert/data/_phase2_stubs.py` |
| Bonus | Unit tests: 25 new tests, all pass (58/58 total) | `tests/test_prereg.py`, `tests/test_predictors_v2_stubs.py`, `tests/test_phase2_dataset_stubs.py` |
| Bonus | Scoop-watch log + CHANGELOG + version bump to 0.2.0rc1 | `lit_notes/scoop_watch.md`, `CHANGELOG.md`, `pyproject.toml`, `src/confpert/__init__.py` |
| Bonus | README updated with Phase 2 callout + verifier usage | `README.md` |

### Phase 2A — remaining for user (requires explicit human action)

| Item | Why | Action |
|---|---|---|
| **Local git commit of Phase 2A deliverables** | Pre-reg YAML's `lock.git_commit_yaml` field needs a real commit SHA | `git add` + `git commit -m "Phase 2A lock: pre-reg v2 + verifier + stubs"` ; then re-run `python -m confpert.cli prereg emit-hashes --prereg paper_neurips_dnb_2026/preregistration_v2.yaml` |
| **OSF mirror upload** | Permanent timestamped registration (scoop defense + reviewer-defensible) | Follow `paper_neurips_dnb_2026/osf_upload_manifest.md` steps |
| **GitHub repo push (private)** | Off-site backup + collaborator access | User-authorized only |

### Phase 2B — predictor expansion (user-supervised Modal work)

Per pre-reg v2 §1.1 priority order: scGPT → scFoundation → Geneformer → CellFM. (GET demoted to optional per §14.5.)

| Phase | Dates | Milestone |
|---|---|---|
| **Phase 2B begins** | once user authorizes Modal launches | scGPT wrapper development + Modal entrypoint + integration tests |

See `PHASE_2_PLAN.md` §7 Track A row 2 for full Phase 2B-2F schedule.

End of Phase 2A: pre-reg v2 locked, machine-verifiable, dry-run verified. Phase 2B onwards may not modify pre-reg v2.
