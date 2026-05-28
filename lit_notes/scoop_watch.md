# ConfPert Phase 2 Scoop Watch Log

Weekly literature scan for work that could scoop Phase 2's contribution. Started 2026-05-25.

**Per `PHASE_2_PLAN.md` §14.2:** Mon mornings, search arXiv / bioRxiv / OpenReview for any 2026-Q2 / Q3 work that would scoop:

1. The K2 v2 cell-line / architecture-family covariate finding
2. The multivariate-CP-for-perturbation framing
3. Pre-registration-as-transferable-methodology for ML benchmarks
4. The distribution-aware single-cell perturbation benchmark scope

**Search queries to run each Monday:**

- `arxiv.org/abs +"single-cell" +"conformal"` last 7 days
- `biorxiv.org +"distribution-aware" +"perturbation"` last 7 days
- `openreview.net +"single-cell" +"benchmark" +"calibration"` open submissions
- Google Scholar: `"perturb-seq" "benchmark" "foundation model"` last week
- Semantic Scholar: papers citing Ahlmann-Eltze 2025 / Csendes 2025 / Li 2025 last week
- bioRxiv: `single-cell perturbation prediction benchmark` last 7 days

---

## 2026-05-25 — baseline scan

Existing related work (already known at Phase 2A.4 lock-in):

| Paper | Authors | Year | Type | Threat level | Notes |
|---|---|---|---|---|---|
| Deep-learning-based gene perturbation effect prediction does not yet outperform simple linear baselines | Ahlmann-Eltze et al. | 2025 *Nature Methods* | Saturation report | LOW (saturation thesis is now community consensus, ConfPert's contribution is calibration angle) | 102 citations |
| Benchmarking foundation cell models for post-perturbation RNA-seq prediction | Csendes et al. | 2025 *BMC Genomics* | Saturation report | LOW | 40 citations |
| Simple controls exceed best deep learning algorithms... | Wong et al. (Pfizer) | 2025 *Bioinformatics* | Saturation report | LOW | 19 citations |
| Benchmarking Transcriptomics Foundation Models for Perturbation Analysis: one PCA still rules them all | Bendidi et al. | 2024 | Benchmark | MEDIUM (parallel benchmark, smaller scope) | 29 citations |
| A Systematic Comparison of Single-Cell Perturbation Response Prediction Models | Li et al. | 2025 bioRxiv | Benchmark | **HIGH** (9 models × 17 datasets × 24 metrics — larger scope) | 19 citations. ConfPert's edge: pre-reg + conformal coverage |
| Multivariate Conformal Prediction using Optimal Transport (otcp2025) | Klein et al. | 2025 ArXiv | CP method | LOW (cited in v1 paper) | 20 citations |
| Optimal Transport-based Conformal Prediction | Thurin et al. | 2025 ArXiv | CP method | LOW (parallel to otcp2025) | 18 citations |
| Minimum Volume Conformal Sets for Multivariate Regression | Braun et al. | 2025 ArXiv | CP method | LOW | 25 citations |
| Beyond Uncertainty Sets: Optimal Transport for Multivariate CPDs | Ndiaye | 2025 ArXiv | CP method | LOW | 1 citation |
| Generalized Conformity Scores for Multi-Output Conformal Regression | Dheur et al. | 2025 ArXiv | CP method | LOW | 20 citations |
| A Kernel Nonconformity Score for Multivariate Conformal Prediction | Meyer et al. | 2026 | CP method (MMD-based!) | **MEDIUM** (overlaps with one of ConfPert's 6 discrepancies) | 0 citations |
| Generative Conformal Prediction with Vectorized Non-Conformity Scores | Zheng et al. | 2024 | CP for generative | **MEDIUM** (overlaps with ConfPert's sample-from-generative-model angle) | 4 citations |
| CP for Generative Models via Adaptive Cluster-Based Density Estimation | Yang et al. | 2026 | CP for generative | MEDIUM | 0 citations |
| Partial Causal Structure Learning for Valid Selective Conformal Inference under Interventions | Asiaee et al. | 2026 | CP for perturbation data | **HIGH** (directly applies CP to Replogle K562 CRISPRi perturbation data) | 0 citations. Different angle (selective inference for descendant identification), but in the same space |
| Pre-registration for Predictive Modeling | Hofman et al. | 2023 ArXiv | ML pre-reg methodology | LOW (only 6 citations; ConfPert builds on this, not scooped) | The closest prior art for the "pre-reg-as-transferable-methodology" thread |

**Top scoop risks:**
1. **Li 2025** — bigger benchmark. Our edge: pre-reg + conformal + transferable pre-reg template.
2. **Asiaee 2026** — CP for Replogle perturbation data. Our edge: distinct angle (population-level calibration, not selective descendant inference) + larger predictor + dataset scope.
3. **Meyer 2026** — MMD-based nonconformity. Our edge: we use MMD as one of 6 discrepancies, not the central methodological claim.
4. **Zheng 2024** — generative CP with vectorized scores. Our edge: we apply to single-cell specifically + pre-reg + K2 v2 cell-line findings.

**Defensive positioning:** All four scoop-risk papers cited explicitly in Phase 2 paper. Frame ConfPert's contribution honestly as (a) the empirical K2 v2 finding (cell-line context dominates capacity), (b) the pre-reg + machine-verifiable artefact, (c) the falsification-intervention design, NOT methodological novelty in CP itself.

---

## Weekly scan template

```markdown
## YYYY-MM-DD — week N

### New papers (citation, threat level, action needed)

- [Citation 1] (LOW/MEDIUM/HIGH) — note + action

### Existing papers — citation count change
- Paper X went from N → M citations
- New citations of relevance: ...

### Action items
- [ ] Cite new work in Phase 2 paper section X
- [ ] Update OSF registration with new related work list
- [ ] Email author of paper Y to flag overlap / coordinate
```

---

## Citation-monitoring source set

Add new papers to the table above; do not delete old entries (preserved for historical defense).
