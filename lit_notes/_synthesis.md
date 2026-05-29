# ConfPert Phase 0 Synthesis (revised, empirical-only)

**Status:** Phase 0 stop point. No project code until this document is approved.
**Scope:** 34 papers read in full across four tiers, plus Cauchois 2024 (Tier 1.5 supplement). All notes in `lit_notes/`.
**Revision:** 2026-05-03. Empirical-only pivot. No new theorems.

---

## (a) Gap proof

The single-cell perturbation prediction field has reached a documented crisis, and no published method provides finite-sample coverage guarantees over distributional discrepancies on cell populations. Three groups of evidence:

**Group 1, mean-only metrics fail.** Ahlmann-Eltze, Huber and Anders 2025 (`ahlmann_2025.md`) show that scGPT, scFoundation, scBERT, Geneformer, UCE, GEARS, and CPA all underperform a low-rank bilinear ridge baseline ($Y_{\text{train}} = G W P^T + b$, $K=10$ PCs of pseudo-bulk, $\lambda = 0.1$) on Norman, Replogle K562, Replogle RPE1, and Adamson under L2 and Pearson-delta. Csendes et al. 2025 BMC Genomics (`csendes_2025.md`) replicate the headline: Train Mean beats scGPT and scFoundation across all four Perturb-seq datasets. Ramakrishnan et al. 2025 (`ramakrishnan_2025.md`) extend the critique by training a per-gene marginal histogram model with Wasserstein-1 loss and outperforming GEARS on distributional structure beyond the mean.

**Group 2, calibration is the missing axis.** The Virtual Cell Challenge 2025 wrap-up (`vcc_2025.md`) flagged metric design and biological-relevance criteria as next-year priorities, observed that "almost all models performed worse than baseline on MAE," and awarded the Generalist Prize to a flow-matching generative model for heterogeneous responses. Cell-Eval 2025 (`cell_eval_2025.md`), the Arc Institute evaluation backbone built into STATE, has 17 registered metrics. Five of ConfPert's six target distributional discrepancies are missing from Cell-Eval (Wasserstein-1, Kolmogorov-Smirnov, MMD, bimodality coefficient match, variance ratio). Cell-Eval reports zero conformal or coverage diagnostics. PerturBench 2024 (`perturbench_2024.md`) explicitly admits the population-level gap and points to MMD/Wasserstein as future work.

**Group 3, the data-side endorses distributional framing.** Replogle et al. 2022 Cell (`replogle_2022.md`) chose nonparametric distributional tests (Anderson-Darling on z-scored single genes, permuted energy distance on top PCs) explicitly because of "incomplete penetrance" and per-cell heterogeneity. Norman et al. 2019 Science (`norman_2019.md`) document principal-curve heterogeneity within identical perturbation conditions. Unlasting / Doloris 2025 (`unlasting_2025.md`) states the bimodality phenomenon directly: "gene expression under the same perturbation often varies significantly across cells, frequently exhibiting a bimodal distribution that reflects intrinsic heterogeneity. This renders expectation-based metrics unreliable."

**The gap.** No published method provides finite-sample, distribution-free coverage guarantees on perturbation cell-population discrepancies. Conformal prediction is the obvious tool, but the conformal-prediction community and the perturb-seq community have not overlapped at the level of finite-sample-coverage methods. ConfPert closes this gap as a benchmark, library, and discovery pipeline.

## (b) Predictor selection

ConfPert wraps eight predictors stratified by parameter count for the K2 capacity-vs-calibration analysis. Selection criterion: native sample-producing API (or a documented sample-augmentation path), public code or trivial reimplementation, and span of parameter counts from $0$ to $10^9$.

| # | Predictor | Param count (approx) | Source | Predict() output |
|---|-----------|----------------------|--------|------------------|
| 1 | Mean predictor | 0 (lookup table) | trivial | per-pert observed mean |
| 2 | Bilinear ridge | $10^4$ | Ahlmann-Eltze 2025 eq. 1, 3 | per-pert mean |
| 3 | scGen | $\sim 10^6$ | `theislab/scgen` | per-cell sample (encode-add-decode) |
| 4 | CPA | $\sim 5 \times 10^6$ | `theislab/cpa` (`cpa-tools`) | parameterized distribution, sample explicitly |
| 5 | biolord | $\sim 2 \times 10^7$ | `nitzanlab/biolord` | per-cell mean and variance, sample with calibrated noise |
| 6 | sVAE+ / SAMS-VAE | $\sim 10^7$ | `insitro/sams-vae` | full posterior samples (NB decoder) |
| 7 | GEARS (uncertainty mode) | $\sim 5 \times 10^7$ | `cell-gears` v0.1.2 | (mean, sigma_hat^2) tuple, marginal Gaussian sample |
| 8 | STATE | $\sim 6 \times 10^8$ | `ArcInstitute/state` | sample-level cell-set output (B, S, output_dim) |

This eight-predictor stratification is the K2 instrument. We test the pre-registered hypothesis that calibration error increases with parameter count.

**scDFM is not in the K1 sweep.** Treat as related-work citation only. Saturated subfield, low marginal value to the sweep, and complicates the calibration-vs-capacity story. CellOT, CellFlow, Squidiff, PerturbDiff are excluded for the same reason plus per-dataset retrain compute cost.

**Adamson 2016 is not in the K1 main sweep.** Sparse, appendix only.

## (c) Conformal method selection (off-the-shelf only)

ConfPert applies four published conformal procedures to the perturbation cell-population setting. No new theorem.

**Per-gene marginal head.** Romano-Patterson-Candes 2019 CQR conformity score $E_i = \max\{q_{\alpha_{lo}}(X_i) - Y_i,\ Y_i - q_{\alpha_{hi}}(X_i)\}$ for quantile-band coverage on each gene marginally; Izbicki-Shimizu-Stern 2022 CD-split / HPD-split for density-level conformal regions on bimodal genes. Off-the-shelf.

**Per-perturbation marginal head.** Energy-distance and MMD-RBF as conformity scores against the held-out observed perturbed-cell-population. Calibration follows split conformal: compute scores on calibration perturbations, take the $\lceil(1-\alpha)(n+1)\rceil$ empirical quantile, accept new perturbation predictions whose score is below the threshold. Off-the-shelf.

**Per-population joint head.** OT-CP 2025 (`otcp_2025.md`, Klein-Bethune-Ndiaye-Cuturi) Brenier-merge to univariate norm with Lemma 3.5 finite-sample empirical-radius. Off-the-shelf, applied to perturbation cell populations for the first time.

**Subgroup-conditional head.** Cauchois-Gupta-Duchi 2024 JMLR (`cauchois_2024_kywyk.md`) for subgroup-conditional confidence sets. Subgroups are the predicted bimodal modes from the predictor's mixture density (responder vs non-responder). Off-the-shelf.

**Risk control (K2 procedure).** Bates et al. 2021 RCPS Theorem 1 with the Waudby-Smith-Ramdas variance-adaptive UCB on bounded losses. Off-the-shelf.

**Distribution shift.** Tibshirani-FoygelBarber-Candes-Ramdas 2019 NeurIPS weighted CP for the K562 calibration / RPE1 test split, with $w(x)$ estimated by a domain classifier. Wasserstein-Regularized CP 2025 ICLR (Xu et al.) as complement for general shift. Off-the-shelf.

**No proofs. No new score classes. No new coverage notions.** All conformal machinery is cited verbatim from prior art. The novelty is in the application, the benchmark, the empirical findings, and the discovery pipeline.

## (d) Three knockouts, all empirical

### K1 (benchmark + library)

First systematic finite-sample conformal-coverage benchmark for single-cell perturbation prediction. Eight wrappable predictors stratified by parameter count, four datasets (Norman 2019, Replogle K562, Replogle RPE1, Tahoe-100M cross-cell-line subset), three split types (within-perturbation, held-out perturbation, cross-cell-line), six distributional discrepancies (KS, W1, energy, MMD-RBF, bimodality coefficient match, variance ratio).

**Deliverables:** pip-installable `confpert` library; Cell-Eval plug-in registering the discrepancies that are missing from Cell-Eval at Phase 1 lock-in; public leaderboard table; reproducibility scripts on Modal.

**Pre-Phase-1 verification (engineering, not lit-reading):** before Phase 1 begins, run a Cell-Eval current-state audit. Concrete steps:

1. Clone `ArcInstitute/cell-eval` at HEAD on the day Phase 1 starts.
2. Open `src/cell_eval/metrics/_impl.py` and grep for the strings `wasserstein`, `kolmogorov`, `mmd`, `bimodal`, `variance_ratio`, `energy`. Record each as present or absent.
3. Document the commit hash, the date, and the present-vs-missing list in `baselines/cell_eval_audit.md`.
4. The plug-in registers only the metrics that are actually missing at lock-in time.

Phase 0 reading found energy distance present as a Pearson correlation across perturbations and the other five missing. Cell-Eval is actively updated; if Cell-Eval has added W1 or MMD between Phase 0 and Phase 1 lock-in, the metric contribution shrinks but the conformal calibration layer (the load-bearing piece) is unaffected. Do not panic if Cell-Eval pre-empts metrics. Do panic if we ship duplicates.

**Sample-generation specification for point-estimate predictors.** Mean predictor and Bilinear ridge produce per-pert mean expression vectors with no native cell-level variance. To compute distributional discrepancies (KS, W1, energy, MMD, bimodality, variance ratio) against them, we must add noise. The noise model is part of the predictor specification, not a confound. We pre-register four noise variants and report all four:

- **Variant A, no noise (point-mass at the mean).** Variance ratio $= 0$ by construction. KS, W1, energy degenerate. Reported as the floor.
- **Variant B, isotropic Gaussian.** $\hat{X} = \mu + \sigma_{\text{global}} \epsilon$ with $\sigma_{\text{global}}$ the pooled per-pert std on the calibration arm, $\epsilon \sim \mathcal{N}(0, I)$.
- **Variant C, per-gene marginal Gaussian.** $\hat{X}_j = \mu_j + \sigma_j \epsilon_j$ with per-gene calibration std. Equivalent to the HetPert noisy-mean baseline.
- **Variant D, full empirical-covariance Gaussian.** $\hat{X} = \mu + \Sigma^{1/2} \epsilon$ with $\Sigma$ the calibration covariance. Hits variance ratio plus second-moment correlations correctly; misses bimodality.

Reviewer-defensible because the K1 sweep then compares "predictor + noise model V" against native sample-producing predictors on equal footing. The noise model is a documented choice with a sensitivity analysis, not a hidden hyperparameter.

**Pre-registered:** report calibration error (target minus achieved coverage) and effective set size for every (predictor, dataset, split, discrepancy, $\alpha$, noise-variant for predictors 1-2) cell. No moving goalposts. Lock baselines in `baselines/results.json` before any tuning.

### K2 (calibration-vs-capacity AND calibration-vs-data, pre-registered)

**Two pre-registered hypotheses, both run, whichever lands becomes the headline.** Single-axis pre-registration (parameter count alone) is fragile because the broader calibration literature (Guo 2017 temperature scaling, deep ensembles, LLM calibration) documents that larger models trained on more data often have BETTER calibration after correction. STATE has 600M parameters AND 270M training cells AND a different objective function. Three confounds, one axis, locks us into a hypothesis the literature already refutes for some settings. We pre-register both directions and let the data decide.

**H1 (capacity hypothesis):** calibration error increases with parameter count at fixed nominal coverage. Operationalized as Spearman correlation between $\log(\text{params})$ and $|1 - \alpha - \text{achieved coverage}|$ across the eight predictors at $\alpha \in \{0.05, 0.10, 0.20\}$. Success: Spearman $\rho > 0.5$ with $p < 0.05$ on at least three of the four datasets, permutation test on predictor labels.

**H1b (data hypothesis):** calibration error decreases with training-set size at fixed nominal coverage. Operationalized as Spearman correlation between $\log(N_{\text{train cells}})$ and $|1 - \alpha - \text{achieved coverage}|$, expected sign negative. Success: Spearman $\rho < -0.5$ with $p < 0.05$ on at least three of the four datasets.

**Pre-registered outcome dispositions:**
- Both H1 and H1b land $\to$ richer story: capacity hurts AND data helps. Foundation-model designers should prioritize data scale over parameter scale for distributional fidelity.
- Only H1 lands $\to$ Ahlmann-Eltze-style headline: bigger is worse for calibration.
- Only H1b lands $\to$ data-scale-driven calibration improvement; pairs with VCC's Tahoe-100M push.
- Neither lands $\to$ null result. Publish honestly. K1 + K3 carry the paper.

**Pre-registered failure plan:** lock both hypotheses with timestamps before Phase 1 first run. Do not reframe post-hoc. If neither $\rho$ exceeds threshold but a third-axis pattern emerges (e.g., training-objective family clusters), report as exploratory in the discussion section, not as a pre-registered finding.

**Why this is a citable finding under any of the three positive scenarios:** pairs naturally with Ahlmann-Eltze 2025 ("DL doesn't beat linear baselines on means") and with the VCC 2025 wrap-up ("flow-matching for heterogeneous responses is the future"). Two of the three positive scenarios provide a clean axis-of-improvement directive for the field. The double-pre-registration also doubles the probability of a publishable empirical finding without sacrificing rigor.

### K3 (discovery K3 with orthogonal computational validation, no wet lab)

Reframe K3 from AUROC delta to a discovery contribution.

**Pipeline:**

1. Train ConfPert wrappers on Tahoe-100M observational arm. Predict cell-population responses for the 49 selective-and-predictable non-oncology compounds in Corsello 2020 PRISM (Pearson $r > 0.4$ in their primary screen, per `prism_2020.md`).
2. For each compound, extract calibrated bimodal subpopulations using Cauchois 2024 subgroup-conditional sets on the predicted-bimodal-mode partition. Tag the resistant mode as the predicted subpopulation maintaining baseline transcriptomic state under drug exposure.
3. Cross-reference resistant-subpopulation gene signatures against three orthogonal computational sources:
   - **DepMap CRISPR essentiality (`depmap_2017.md`):** do resistant-subpopulation marker genes show selective dependency in non-resistant cell lines?
   - **PRISM secondary screen viability AUC:** does the resistant-subpopulation transcriptomic signature align with the cell lines that retain viability under the drug?
   - **Hallmark MSigDB pathways:** does the signature enrich for known resistance pathways (mTOR, MYC, EMT, oxidative phosphorylation, hypoxia, p53)?
4. Repeat for uncalibrated baselines (raw GEARS predict, raw scGPT-architecture, raw STATE without conformal calibration). Compare hit rates.

**Multiple-comparisons correction (load-bearing).** The naive criterion "at least three signatures with Hallmark $p < 0.01$" without correction is loose enough to hit by chance over the full grid (49 compounds $\times$ ~50 Hallmark pathways $\times$ ~500 PRISM cell lines $\times$ ~17K DepMap genes). We pre-register Benjamini-Hochberg FDR correction at $q = 0.05$ over the joint compound $\times$ Hallmark-pathway grid for the enrichment test, plus a separate BH correction at $q = 0.05$ over the compound $\times$ DepMap-gene grid for the dependency-alignment test. A signature counts only if both corrected tests pass for that compound.

**Pre-registered success criterion (BH-corrected):** at least three drug-cell-line pairs where ConfPert's calibrated bimodal subpopulation produces a resistance signature with Hallmark enrichment $q < 0.05$ (BH over the full compound $\times$ pathway grid) AND DepMap selective-dependency alignment $q < 0.05$ (BH over the compound $\times$ gene grid), AND uncalibrated baselines fail at least two of those three at the same corrected thresholds.

**Pre-registered failure plan:** if K3 fails the BH-corrected criterion, report honestly. The K1 + K2 contributions still publish as a Datasets and Benchmarks track paper. Do not switch to looser uncorrected thresholds post-hoc.

**Why this is a discovery contribution:** the deliverable is "we surfaced N specific drug-resistance signatures that the field's current uncalibrated approach systematically misses, validated by two orthogonal computational sources at FDR-corrected significance." This is a discovery + methods paper, not just a methods paper. Cell Systems and Nat Methods Brief Comms accept this format. NeurIPS main track values it.

## (e) Narrative arc

**Section 1: Field in crisis.** Open with Ahlmann-Eltze 2025, Csendes 2025, Ramakrishnan 2025, Unlasting / Doloris 2025, and the VCC 2025 wrap-up. Cite Replogle 2022 and Norman 2019 as the data-side endorsement of the distributional framing.

**Section 2: ConfPert framework.** Define the six distributional discrepancies, the four conformal head levels (per-gene CD-split / HPD-split, per-perturbation energy / MMD, per-population OT-CP, subgroup-conditional Cauchois 2024), and the risk-control extension (Bates 2021). All procedures are off-the-shelf with verbatim citations. State the application context: cell-population perturbation prediction.

**Section 3 (K1, benchmark + library):** the 8-predictor x 4-dataset x 3-split sweep with calibration-error tables and effective-set-size tables. Pip library and Cell-Eval plug-in shipping with the paper.

**Section 4 (K2, calibration-vs-capacity):** the pre-registered empirical claim. Report Spearman correlations with permutation-test p-values. Plot calibration error vs $\log(\text{params})$ for each $\alpha$ and each dataset. Discussion: implications for the foundation-model-scaling argument in single-cell biology.

**Section 5 (K3, discovery):** the PRISM drug-resistance discovery pipeline. Report the surfaced drug-resistance signatures with their orthogonal-validation evidence. Show specific cases where uncalibrated baselines miss the signature.

**Section 6: Limitations.** Cell-line transcriptomics-to-bulk-PRISM bridge constraint. CPA / biolord / sVAE+ / scGen retrain per dataset. Cross-cell-line miscoverage observed empirically as a binding constraint. Marginal coverage is not conditional coverage.

**Section 7: Discussion.** Connection to VCC ecosystem (Cell-Eval likely adopts ConfPert metrics in v2). Forward-looking: active virtual cell experimentation via prediction-set-width acquisition (deferred future work).

### Explicit positioning vs PerturbDiff, scDFM, CellFlow, Squidiff, Unlasting

Single-paragraph differentiation, repeated in Section 1 and Section 7 to prevent reviewer conflation.

**PerturbDiff (Yuan et al. 2026, ICML 2026 sub), scDFM (Yu et al. 2026 ICLR), CellFlow (Klein et al. 2025 bioRxiv), Squidiff (He et al. 2025 Nat Methods), Unlasting / Doloris (Chi et al. 2025) all argue for distributional fidelity and produce sample distributions for cell populations.** None of them provides finite-sample coverage guarantees on those samples. Their distributional losses are training-time regularizers (MMD on flow endpoints in scDFM; squared RKHS norm = MMD$^2$ in PerturbDiff; energy-distance training in STATE) with no finite-sample-bounded calibration guarantee at evaluation time. ConfPert is orthogonal: we provide a model-agnostic conformal wrapper that produces finite-sample coverage on any sample-producing predictor, including these. The 2025-2026 distributional-perturbation-model wave produces calibration targets, not calibration guarantees; ConfPert provides the latter on top of the former.

**Why these are excluded from the K1 sweep.** Compute realism. CellOT, CellFlow, Squidiff, PerturbDiff lack public per-dataset checkpoints; per-dataset retraining for four datasets times each model exceeds the $5-8K compute budget. scDFM has public Norman + ComboSciPlex checkpoints but adds saturation noise to the K2 calibration-vs-capacity story (it is a 4-layer 512-hidden Transformer with no clean parameter-count placement on the K2 axis). The K1 sweep prioritizes parameter-count stratification (Mean $\to$ Bilinear ridge $\to$ scGen $\to$ CPA $\to$ biolord $\to$ sVAE+ $\to$ GEARS $\to$ STATE) over predictor-zoo breadth. Reviewer-defensible because the noise-model sensitivity analysis in K1 covers the "what if your wrappers aren't representative" concern, and the explicit positioning above covers the "why not these specifically" concern.

### Calibration-aware foundation-model design directive (Section 7 paragraph)

The K2 outcome implies a design directive for the next wave of perturbation foundation models. State it explicitly in the discussion regardless of which way K2 lands:

- **If H1 lands (capacity hurts calibration):** "Foundation-model designers should prioritize calibration-aware training, including conformal-aware losses and post-hoc temperature scaling at scale, before further parameter scaling. The Ahlmann-Eltze 2025 critique of mean accuracy generalizes to calibration."
- **If H1b lands (data helps calibration):** "Data-scale investments such as Tahoe-100M and the Arc Virtual Cell Atlas translate directly to better calibration. The field's 100M-cell push is justified for distributional fidelity, not just mean accuracy."
- **If both land:** combine: "Calibration improves with data, degrades with parameter count at fixed data; the design lever is data scale, not capacity."
- **If neither lands:** "We report a clean null on the parameter-count and training-set-size axes; training-objective-family analysis (MSE vs flow-matching vs energy-distance vs diffusion) is the next research direction."

Forward-looking design directives are what reviewers remember after the calibration tables blur. One paragraph, costs nothing, raises citation potential by giving the next foundation-model paper a clean reference target.

## (f) Top three risks

**Risk 1: Both K2 hypotheses fail.** Two pre-registered hypotheses (H1 capacity-hurts and H1b data-helps) increases the probability of a publishable empirical finding. If both fail, the headline evaporates. **Mitigation:** pre-register both hypotheses with timestamps before Phase 1 first run. Report whichever lands. If neither lands, report as a null result and lean on K1 (benchmark) plus K3 (discovery) plus the small-proposition aggregation result to carry the paper. Do not reframe to a third axis post-hoc.

**Risk 2: K3 finds zero BH-corrected drug-resistance signatures.** PRISM is bulk barcode viability; the bridge from per-cell transcriptomic prediction to PRISM resistance is non-trivial. The BH-corrected criterion is appropriately strict and may produce zero hits. **Mitigation:** anchor to Corsello 2020's own bimodality-coefficient selectivity index. Pre-register the BH-corrected success criterion at $N \geq 3$ signatures with Hallmark $q < 0.05$ and DepMap dependency $q < 0.05$. If we get zero, K3 publishes as a falsifiable null result. The framework + benchmark + small proposition still carries the paper.

**Risk 3: Compute and engineering realism on eight predictors.** STATE 600M-parameter inference plus seven smaller predictors across four datasets and three splits is roughly 96 distinct calibration runs, multiplied by four noise variants for the two point-estimate predictors. Per-dataset retraining for scGen, CPA, biolord, sVAE+ adds ~50 training jobs. **Mitigation:** lock baselines in `baselines/results.json` before any framework development. Use STATE SE-600M public checkpoint plus the lightweight-predictor stack. Total compute estimate $3 - 5K Modal/Lambda H100 spot, within the $5 - 8K budget with margin. Cell-Eval audit step is engineering, $0 cost.

**Risk 4 (added): Cell-Eval already adopts our discrepancies before we ship.** Cell-Eval is actively updated by Arc Institute. If they add MMD or W1 in a release between now and our Phase 1 lock-in, our metric contribution shrinks. **Mitigation:** the Cell-Eval audit step is part of the K1 deliverables. The conformal calibration layer (per-gene, per-perturbation, per-population, subgroup-conditional) is the load-bearing contribution, NOT the metric set. Five new metrics shipping with the library is a useful bonus, not the headline.

## Honest novelty grade with this revision

**7.5 - 8.25.**

- 7 baseline because the calibration framework is a domain port + combination of off-the-shelf conformal procedures.
- +0.5 from K2 if either H1 (capacity hurts) or H1b (data helps) lands. The double-hypothesis pre-registration roughly doubles the probability of a publishable empirical finding without sacrificing rigor.
- +0.5 if K3 surfaces $\geq 3$ BH-FDR-corrected drug-resistance signatures that uncalibrated baselines miss. Discovery contribution on top of methods contribution.
- +0.25 from explicit positioning vs PerturbDiff / scDFM / CellFlow / Squidiff / Unlasting. Reframes ConfPert as the calibration wrapper for the 2025-2026 distributional-perturbation-model wave, not a competitor to it.

The earlier draft included a six-discrepancy aggregation proposition for a +0.25 bump. Dropped per Phase 0 review: the project does not have an external theorist available within a one-week window for sanity-check, and reviewer-rejection risk on unverified math outweighs the marginal novelty. Paper-quality cost is negligible per the same review. If K1+K2+K3 land cleanly, a verified aggregation theorem is the natural V2 addition.

Not 9 - 10. Empirical-only papers in conformal prediction rarely hit paradigm-shift territory. But 8 is a strong main-track NeurIPS / ICML pitch and gives Cell Systems / Nat Methods Brief Comms a real shot if K3 lands. The 9-floor would require either a wet-lab validation arm (not in scope) or a genuinely new conformal theorem (deferred to V2).

## Real-impact ceiling

- **NeurIPS / ICML / ICLR main track:** medium-to-high probability. K2 + K3 do most of the work.
- **NeurIPS Datasets and Benchmarks track:** very high probability. The benchmark + library is the natural fit.
- **Cell Systems / Nat Methods Brief Comm:** medium probability with clean K3.
- **VCC ecosystem foundation:** very high. Cell-Eval likely integrates ConfPert calibration metrics within 6 - 12 months if it ships first.
- **Citation potential over 3 years:** medium-to-high. Benchmark + surprising empirical finding + discovery is a stronger citation triple than any single piece.

## Phase 0 stop point

This revised synthesis covers (a) gap proof, (b) predictor selection, (c) off-the-shelf conformal-method choice, (d) three empirical knockouts with pre-registered success criteria, (e) narrative arc, (f) top three risks. All 34 paper notes plus Cauchois 2024 supplement are at `lit_notes/<slug>.md`.

No project code. Awaiting approval before Phase 1.
