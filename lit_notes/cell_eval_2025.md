# Cell-Eval 2025 (Arc Institute, in STATE paper)

**Venue:** Built into Adduri et al. 2025 STATE paper (bioRxiv 2025.06.26.661135). Released as a standalone package `cell-eval` (PyPI + GitHub `ArcInstitute/cell-eval`, latest v0.7.0). The STATE manuscript explicitly directs publications using cell-eval to cite the State paper.

## Claim

Cell-Eval is a standardized, registry-based evaluation harness for predictors of single-cell perturbation response. It compares two AnnData objects (predicted vs. real) and computes a fixed suite of metrics covering global expression error, mean-effect correlation, perturbation discrimination by rank, energy-distance correlation, leiden clustering agreement, and a battery of differential-expression metrics (top-K overlap/precision, Spearman LFC correlation, direction match, recall, ROC/PR AUC). It runs DE via the companion `pdex` package on both predicted and real anndatas, then aggregates per-perturbation scores. Cell-Eval is the evaluation backbone Arc uses for the Virtual Cell Challenge (VCC) and STATE benchmarks.

Verbatim from README: "This package provides a comprehensive suite of metrics for evaluating the performance of models that predict cellular responses to perturbations at the single-cell level."

## Method

Cell-Eval ingests `(adata_pred, adata_real)`, a control label, and a perturbation column, then runs differential expression separately on each anndata via `pdex`, and dispatches metrics through a global `MetricRegistry`. Metrics fall in two registered types: `ANNDATA_PAIR` (operate on raw expression / embedding matrices) and `DE` (operate on DE result tables). Results return per-perturbation dicts plus an aggregated `agg_results` table. Profiles select metric subsets: `full`, `minimal`, `vcc`, `de`, `anndata`, `pds`.

Verbatim: "It generally revolves around a *real* anndata and a *predicted* anndata where it measures the general differences between the two across a variety of metrics."

## Metrics (exact, from `_impl.py`, `_anndata.py`, `_de.py`)

ANNDATA_PAIR metrics:

1. `pearson_delta` — Pearson r between predicted and real mean perturbation effect deltas (mean(pert) - mean(ctrl)). Best = 1. Per-perturbation.
2. `mse` — sklearn `mean_squared_error` between predicted and real per-perturbation mean expression vectors. Best = 0. Per-perturbation.
3. `mae` — sklearn `mean_absolute_error` between per-perturbation mean expression vectors. Best = 0. Per-perturbation.
4. `mse_delta` — MSE on perturbation-control delta vectors (effect-space). Best = 0.
5. `mae_delta` — MAE on perturbation-control delta vectors. Best = 0.
6. `discrimination_score_l1`, `_l2`, `_cosine` (PDS) — for each perturbation p, compute pairwise distance (l1/l2/cosine) between pred-effect[p] and all real-effect[q]; rank where the matching index falls; report `1 - rank/N`. Best = 1. L1 forces `embed_key=None` (expression space). Excludes the target gene by default. This is the Virtual Cell Challenge primary metric.
7. `pearson_edistance` — energy distance (Szekely-Rizzo) computed per-perturbation as `2 * mean(d(pert, ctrl)) - mean(d(pert, pert)) - mean(d(ctrl, ctrl))` with euclidean default; then Pearson r between the real-edistance vector and predicted-edistance vector across perturbations. Best = 1.
8. `clustering_agreement` (`ClusteringAgreement` class) — compute per-perturbation centroids in real and pred, build kNN graphs (n_neighbors=15), Leiden cluster real at resolution 1.0 and predicted at resolutions {0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 2.0}, score AMI/NMI/ARI between cluster labelings, return best score. Default metric AMI. Best = 1.

DE metrics (FDR threshold default 0.05; sort by abs fold change):

9. `overlap_at_{N,50,100,200,500}` — Jaccard / set overlap of top-K significant DE genes per perturbation. Best = 1. (When `k=None` it uses the number of real significant genes, "top N".)
10. `precision_at_{N,50,100,200,500}` — precision at top-K. Best = 1.
11. `de_spearman_sig` (`DESpearmanSignificant`) — Spearman correlation between count of significant genes per perturbation in real vs predicted. Single scalar.
12. `de_direction_match` (`DEDirectionMatch`) — for genes significant in real, fraction whose log2FC sign matches in predicted (per perturbation).
13. `de_spearman_lfc_sig` (`DESpearmanLFC`) — Spearman correlation between real and pred log fold-changes restricted to real-significant genes (per perturbation).
14. `de_sig_genes_recall` (`DESigGenesRecall`) — recall = |sig_real intersect sig_pred| / |sig_real| per perturbation.
15. `de_nsig_counts` — diagnostic counts of significant genes (real, pred). No best value.
16. `pr_auc` — average-precision over all genes per perturbation, labels = (real FDR < 0.05), scores = -log10(pred FDR). Best = 1.
17. `roc_auc` — same labels/scores, ROC AUC. Best = 1.

Profiles (`_pipeline/_runner.py`):
- `MINIMAL_METRICS = [pearson_delta, mse, mae, discrimination_score_l1, overlap_at_N, precision_at_N, de_nsig_counts]`
- `VCC_METRICS = [mae, discrimination_score_l1, overlap_at_N]` (this is the Virtual Cell Challenge official trio)
- `pds` profile = `[discrimination_score_l1]`
- `full` = all DE + all ANNDATA_PAIR.

Aggregation: `MetricsEvaluator.compute()` returns `(results, agg_results)` as polars DataFrames. The `cell-eval score` CLI / `score_agg_metrics` normalizes user `agg_results.csv` against a baseline `agg_results.csv` to produce a leaderboard score.

## Datasets / experiments

The harness is dataset-agnostic; in the STATE paper and VCC it is applied to:
- Replogle / Replogle-Nadig (genome-scale Perturb-seq, K562 + Jurkat + RPE1)
- Tahoe-100M (chemical perturbation, ~100M cells)
- Parse / scPerturb-style consortium datasets
- Virtual Cell Challenge held-out test split (H1 hESC perturbation set)

Splits supported via the companion `cell-load` package: `zeroshot` (hold out entire cell types) and `fewshot` (hold out perturbations within cell types).

## What we steal for ConfPert

Cell-Eval is the leaderboard-grade reporting frontend; ConfPert wraps any sample-producing predictor with conformal calibration and reports through Cell-Eval to stay comparable.

Direct adoption:
- `discrimination_score_l1` (PDS) - we adopt as our identifiability check; it is the VCC primary.
- `overlap_at_N`, `precision_at_N` - top-DE recovery, we report alongside our distributional metrics.
- `de_spearman_lfc_sig`, `de_direction_match` - DE-side correctness.
- `pearson_edistance` - Cell-Eval already covers energy distance (correlated form). We retain our raw per-perturbation E-distance discrepancy because Cell-Eval reports a correlation across perturbations rather than per-perturbation energy distance values; ours is finer-grained.

Metric gap (additions ConfPert contributes on top of Cell-Eval):
- **Wasserstein-1 (W1) per gene** - Cell-Eval has no W1 metric. We add it (per-gene 1D W1, mean across genes). Justifies marginal-distribution fidelity beyond mean-vector MAE.
- **Kolmogorov-Smirnov (KS) per gene** - Cell-Eval has no KS. We add it as a non-parametric tail-sensitive marginal check.
- **MMD-RBF (multivariate)** - Cell-Eval has no MMD. We add it for joint-distribution discrepancy with kernel bandwidth set by median heuristic.
- **Bimodality coefficient match** - Cell-Eval has no shape statistic for non-Gaussian marginals. We add per-gene bimodality coefficient comparison (real vs pred) since perturbed expression often becomes bimodal (on/off knockdown) and MAE/MSE miss this entirely.
- **Variance ratio** - Cell-Eval has `mse_delta` on means but no second-moment fidelity check. We add per-gene Var(pred)/Var(real) so over-smoothed predictors (a known failure mode of point estimators) get penalized.

Net: ConfPert's six discrepancies (KS, W1, energy, MMD-RBF, bimodality match, variance ratio) intersect Cell-Eval at energy distance only. The other five are genuine additions.

## What we wrap

Not a model - an evaluation harness. ConfPert's pipeline contract: produce sample-level AnnData predictions per (cell_type, perturbation), run `cell-eval run --profile full` for the leaderboard view, and run our `confpert.eval` for the six discrepancies plus conformal coverage diagnostics (PI width, marginal coverage at level alpha, conditional coverage by cell type).

## Failure modes / caveats

- **No native distributional discrepancy on marginals.** Cell-Eval's `mse`, `mae`, `pearson_delta`, `mse_delta`, `mae_delta` all collapse cells to per-perturbation mean vectors before comparison. A predictor that perfectly matches the mean but produces zero-variance output (mode collapse) gets perfect MSE/MAE/pearson_delta. Only `pearson_edistance` and `discrimination_score` use cell-level information; both still aggregate.
- **Energy distance reported as correlation, not absolute.** `pearson_edistance` is a Pearson r between real-edist and pred-edist vectors across perturbations. A model that uniformly under-predicts effect magnitude can still correlate at 1.0. ConfPert reports raw per-perturbation energy distance.
- **DE metrics depend on `pdex` and chosen FDR threshold (0.05).** Different DE callers move scores meaningfully; the harness fixes the choice for comparability but it is not invariant.
- **`discrimination_score_l1` ignores embedding spaces by design** (forced `embed_key=None`); l2/cosine variants do honor embeddings.
- **Cell-Eval does not report uncertainty calibration, prediction-interval coverage, or any conformal diagnostic.** All point-prediction.
- **`clustering_agreement` is non-deterministic across runs** (Leiden, kNN graph) unless seeds are fixed.
- **No support for non-Gaussian shape statistics** (skewness, kurtosis, bimodality), which are exactly the regimes where mean-only metrics mislead.

## Code URLs (confirmed)

- `https://github.com/ArcInstitute/cell-eval` (v0.7.0, MIT, 127+ stars)
- `https://github.com/ArcInstitute/state` (model code; calls cell-eval via `state tx predict`)
- `https://github.com/ArcInstitute/cell-load` (dataloaders / split TOMLs)
- Paper: `https://www.biorxiv.org/content/10.1101/2025.06.26.661135v2`
- Manuscript hub: `https://arcinstitute.org/manuscripts/State`
- VCC: `https://virtualcellchallenge.org/`

## Verbatim quotes worth keeping

1. "This package provides a comprehensive suite of metrics for evaluating the performance of models that predict cellular responses to perturbations at the single-cell level."
2. "The metrics are built using the python registry pattern. This allows for easy extension for new metrics with a well-typed interface."
3. "It generally revolves around a *real* anndata and a *predicted* anndata where it measures the general differences between the two across a variety of metrics."
4. From `discrimination_score` docstring: "Best score is 1.0 - worst score is 0.0... Sort by distance (ascending - lower distance = better match)... Find rank of the correct perturbation... Normalize rank by total number of perturbations."
5. From `edistance` docstring and code: "Compute Euclidean distance of each perturbation-control delta... `2 * delta - sigma_x - sigma_y`" (the standard Szekely-Rizzo energy distance), then "pearsonr(d_real, d_pred).correlation" across perturbations.
6. "Any publication that uses this source code should cite the State paper."
