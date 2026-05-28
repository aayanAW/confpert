# Cell-Eval audit (Phase 1 lock-in)

**Audit date:** 2026-05-03 (overnight Phase 1).
**ConfPert preregistration commit:** `c7046e4b4c08106acd52c37e7ea8410e126d2ae4`.

## Cell-Eval HEAD at audit time

- **Repo:** `https://github.com/ArcInstitute/cell-eval`
- **Commit hash:** `2db6b9c65498c2e801cc1cf96f2e42e65e7b5016`
- **Commit timestamp:** 2026-03-24 11:11:42 -0700
- **Last commit message:** "Merge pull request #225 from ArcInstitute/cell-eval-0.7.0"

## Six-discrepancy presence audit

Per the Phase 0 reading and Phase 1 plan, ConfPert defines six distributional discrepancies. The audit greps each of {`wasserstein`, `kolmogorov`, `smirnov`, `ks_`, `mmd`, `bimodal`, `variance_ratio`, `energy`} in `src/cell_eval/metrics/`.

| # | Discrepancy | Cell-Eval HEAD | Notes |
|---|---|---|---|
| 1 | KS | **MISSING** | No `ks_2samp`, `kolmogorov_smirnov`, or `ks_` reference in any metric file |
| 2 | Wasserstein-1 | **MISSING** | No `wasserstein` reference |
| 3 | Energy distance | **PRESENT (different framing)** | `pearson_edistance` reports Pearson correlation on energy distance across perturbations, not the raw per-perturbation absolute values ConfPert uses |
| 4 | MMD-RBF | **MISSING** | No `mmd` reference |
| 5 | Bimodality coefficient match | **MISSING** | No `bimodal` reference |
| 6 | Variance ratio | **MISSING** | No `variance_ratio` reference |

## Cell-Eval registered metrics (17 total)

From `src/cell_eval/metrics/_impl.py`:

**ANNDATA_PAIR (expression-space, 9):**
1. `pearson_delta`
2. `mse`
3. `mae`
4. `mse_delta`
5. `mae_delta`
6. `discrimination_score_l1`
7. `discrimination_score_l2`
8. `discrimination_score_cosine`
9. `pearson_edistance`
10. `clustering_agreement`

**DE (differential-expression, 7+):**
1. `overlap_at_{50, 100, 200, 500, N}` (5 variants)
2. `precision_at_{50, 100, 200, 500, N}` (5 variants)
3. `de_spearman_sig`
4. `de_direction_match`
5. `de_spearman_lfc_sig`
6. `de_sig_genes_recall`
7. `de_nsig_counts`
8. `pr_auc`
9. `roc_auc`

## Conclusion

**5 of 6 ConfPert discrepancies are absent from Cell-Eval HEAD as of `2db6b9c6`.** ConfPert's plug-in adds: KS, Wasserstein-1, MMD-RBF, bimodality coefficient match, variance ratio. Energy distance is present in Cell-Eval but in a Pearson-correlation-across-perturbations framing rather than ConfPert's raw per-perturbation form; ConfPert exposes both but does not duplicate the existing Cell-Eval metric.

## Phase 1 plug-in shipping plan

The `confpert.cell_eval_plugin` module (Phase 1 deliverable) registers 5 new metrics with Cell-Eval's `metrics_registry` at import time, named:

- `confpert_ks`
- `confpert_w1`
- `confpert_mmd_rbf`
- `confpert_bimodality_match`
- `confpert_variance_ratio`

The `confpert_` prefix prevents naming collisions if Cell-Eval upstream adds a same-named metric in a future release. Energy distance is consumed via `cell_eval.metrics_registry.get("pearson_edistance")` rather than re-registered.

## Risk acknowledgment

Cell-Eval is actively maintained. If a release between this audit (2026-05-03) and ConfPert paper submission adds Wasserstein or MMD, the metric contribution shrinks but the conformal calibration layer (per-gene CD-split, per-perturbation marginal, per-population OT-CP, subgroup-conditional Cauchois) is unaffected. Per the synthesis Risk 4: "do not panic if Cell-Eval pre-empts metrics; do panic if we ship duplicates."

## Reproducibility

```bash
cd /Users/aayanalwani/FLM4S/confpert/external
git clone --depth 1 https://github.com/ArcInstitute/cell-eval.git
cd cell-eval && git rev-parse HEAD  # 2db6b9c65498c2e801cc1cf96f2e42e65e7b5016
grep -E "wasserstein|kolmogorov|mmd|bimodal|variance_ratio|energy" src/cell_eval/metrics/*.py
```
