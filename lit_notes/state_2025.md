# Adduri et al. 2025 - STATE (Arc Institute)

**Venue:** bioRxiv preprint, doi:10.1101/2025.06.26.661135 (v1 posted 2025-06-27, v2 update). Arc Institute manuscript page: https://arcinstitute.org/manuscripts/State. Not yet peer reviewed. License CC-BY for the preprint, CC-BY-NC-SA 4.0 for code, custom Arc Research Institute State Model Non-Commercial License for weights. Corresponding author Yusuf H. Roohani; first author Abhinav K. Adduri.

## Claim

STATE is a transformer-based virtual-cell model that predicts cellular responses to genetic, signaling, and chemical perturbations across diverse cellular contexts. Two interlocking modules: State Embedding (SE) for observational cell representation and State Transition (ST) for perturbation-conditioned set-to-set prediction. The paper claims more than 30 percent improvement in perturbation effect discrimination on large datasets, roughly twice the accuracy in identifying differentially expressed genes versus prior baselines (CPA, scGPT, scVI, GEARS, pseudobulk, mean-baselines), and is reported as the first model to consistently beat simple linear baselines on these tasks.

## Method

**ST architecture (state_transition.py).** Llama-style transformer backbone with bidirectional attention. Default config: hidden_dim 768, 8 layers, 12 heads, head_dim 64, intermediate 3072, cell_set_len 512. Small: 672/4/8/128. Large: 1488/6/12/512. Rotary embeddings disabled. Input is a cell-set of length S; per cell, basal expression embedding (HVG vector or SE embedding) plus perturbation embedding are summed in input space, projected to hidden_dim, transformed across the cell-set, and decoded via an MLP. Optional residual prediction (default True): output = decoder(transformer_output + control_cells). Optional add-ons: batch token, confidence token (predicted scalar loss), LoRA, finetune VCI counts decoder.

**ST loss (verbatim from code).** Default `distributional_loss: energy` implemented as `geomloss.SamplesLoss(loss="energy", blur=0.05)`. Options: `sinkhorn`, `mse`, or `se` which is a `CombinedLoss` of `sinkhorn_weight * sinkhorn + energy_weight * energy` (defaults 0.001 and 1.0). Operates on full predicted vs target cell sets per perturbation, so the objective is distribution-matching at the set level rather than pointwise. Optional auxiliary CE losses for batch token and MSE for confidence token. Optional decoder loss when finetuning a counts decoder.

**SE architecture (emb/nn/model.py).** Lightning module `StateEmbeddingModel` with FlashTransformer encoder layers, learnable CLS token, gene-token encoder MLP (Linear, LayerNorm, SiLU). Public checkpoint SE-600M. Gene tokens initialized from ESM2 embeddings of the protein product (configs profile esm2-cellxgene at size 5120; alternatives: esm2_3B-scbasecamp size 2560, evo2-scbasecamp size 4096). Loss options in emb/nn/loss.py: WassersteinLoss, KLDivergenceLoss, MMDLoss, TabularLoss. SE produces per-cell embeddings used as the basal input for ST.

**Tokenization.** Genes are tokens with embedding from a protein language model (ESM2, optionally Evo2 for cross-species). Each cell is a bag of (gene_token, expression_value). Cells are gathered into cell-sets ("sentences") of length 128 to 512 and the transformer attends across cells in a set with bidirectional attention.

**Pretraining and fine-tuning.** SE pretrained on roughly 167M observational human cells (CellxGene + scBaseCamp + Tahoe; SE-167M profile). ST trained on over 100M perturbed cells across 70 cell contexts: Tahoe-100M (chemical), Replogle-Nadig (CRISPRi), Parse-PBMC. Hydra training; TOML splits for zeroshot (held-out cell types) and fewshot (held-out perturbations within a cell type). Preprocessing: normalize, log1p, top-2000 HVGs into `obsm["X_hvg"]`.

**Predict-time output (critical for ConfPert).** ST.predict produces cell-set outputs: given a cell-set of basal cells plus a perturbation, the model returns a same-shape set of predicted post-perturbation cells (B, S, output_dim) where output_dim is either the HVG count (gene space) or SE embedding dim. By default this is sample-level prediction, not a single mean; the energy-distance loss is computed against the empirical perturbed cell-set and so the predictions empirically match the target distribution. There is a `predict_mean` flag that switches to mean-only prediction. The `state tx infer` CLI writes predictions to an h5ad (or .npy of the predictions matrix). The ST output is typically returned in the same gene or HVG basis as the input.

## Datasets / experiments

Training: Tahoe-100M (chemical perturbations across roughly 60 cancer cell lines), Replogle-Nadig CRISPRi Perturb-seq (K562, RPE1, Jurkat, hepatocytes), Parse-PBMC, plus Arc Virtual Cell Atlas observational data totaling about 167M cells for SE. Stated coverage: more than 100M perturbed cells across 70 cellular contexts.

Evaluation: zeroshot held-out cell types and fewshot held-out perturbations within cell types. Reported headline numbers: more than 30 percent improvement in perturbation discrimination, roughly 2x improvement in DE gene identification, and identifies strong perturbations in novel cellular contexts not seen at training time.

## Metrics

Cell-Eval (github.com/ArcInstitute/cell-eval) is the official evaluation harness. Metrics enumerated in `src/cell_eval/metrics/_impl.py`:

Distributional / per-perturbation metrics on (pred, real) AnnData pairs:
- `pearson_delta`: Pearson correlation between predicted and real mean delta-from-control vectors.
- `mse`, `mae`: per-perturbation mean squared and mean absolute error vs control.
- `mse_delta`, `mae_delta`: same but on perturbation-control deltas.
- `discrimination_score_l1`, `discrimination_score_l2`, `discrimination_score_cosine`: normalized rank similarity of each predicted perturbation to its true counterpart under the named distance. This is the "perturbation discrimination" headline metric.
- `pearson_edistance`: Pearson correlation of predicted vs real energy-distance from control across perturbations.
- `clustering_agreement`: agreement between clustering of real vs predicted perturbation centroids.

Differential-expression metrics on a `DEComparison`:
- `overlap_at_50`, `overlap_at_100`, `overlap_at_200`, `overlap_at_500`, `overlap_at_N`: top-k DE gene overlap at fixed k or full ranked list.
- `precision_at_50`, `precision_at_100`, `precision_at_200`, `precision_at_500`, `precision_at_N`.
- `de_spearman_sig`: Spearman on counts of significant DE genes.
- `de_direction_match`: agreement on direction of DE changes.
- `de_spearman_lfc_sig`: Spearman on log-fold-changes of significant genes.
- `de_sig_genes_recall`: recall of significant genes.
- `de_nsig_counts`: counts of significant genes.
- `pr_auc`, `roc_auc`: precision-recall and ROC AUC for significant DE recovery.

DE is computed via Arc's `pdex` package on the AnnData inputs. The CLI exposes a `--profile` flag (e.g. `full`) to select metric suites, and `cell-eval prep` and `cell-eval score` for VCC-formatted submissions and baseline-normalized scoring.

## What we steal for ConfPert

Cell-Eval is adopted as the ConfPert evaluation backbone without modification. Import it, feed `(adata_pred, adata_real)` with `pert_col` and `control_pert`, take `agg_results`. The `discrimination_score_*`, `pearson_delta`, `pearson_edistance`, and `overlap_at_k`/`precision_at_k` metrics map onto ConfPert's distributional-coverage and DE-recovery axes. `clustering_agreement` and `de_direction_match` give complementary signal not redundant with raw distance.

## What we wrap

STATE is wrappable predictor #8. It is the most recent canonical scFM-style perturbation predictor, has public weights (SE-600M, ST checkpoints trained on Tahoe-100M and Replogle-Nadig from the Colab notebooks), and its native predict returns a sample-level cell-set, the exact input format ConfPert needs for distributional conformal calibration. Repo: https://github.com/ArcInstitute/state. CLI: `state tx infer --model-dir ... --adata ... --pert-col ...`.

## Failure modes / caveats

- License. Code is CC-BY-NC-SA 4.0; weights are under a non-commercial Arc license with an Acceptable Use Policy. ConfPert artifacts downstream of STATE inherit non-commercial restrictions.
- Tahoe-100M licensing. Hosted by Vevo Therapeutics; redistribution constraints must be checked before vendoring eval splits.
- Compute. SE-600M public checkpoint is 600M params; inference requires GPU. Official Singularity container assumes large-storage mounts.
- Predict-time outputs are conditional on a basal cell set, so calibration must be conditioned on both perturbation and basal context. The `predict_mean` flag discards the distributional structure that energy-loss training optimizes for.
- Energy-distance training matches distributions consistently but is lower-power than MMD with a tuned kernel for tail discrimination, which matters for ConfPert tail-coverage claims.
- The bioRxiv preprint is not peer reviewed; at access time the PDF and JATS XML are Cloudflare-gated. Methods here are reconciled from the abstract, the Arc news post, GitHub source, and model configs. Headline numbers are from v1; v2 may revise.
- Bidirectional Llama backbone with rotary embeddings disabled. The position-embedding choice for cell-set permutation invariance is worth checking against any equivariance claims.

## Code URLs

- Model and CLI: https://github.com/ArcInstitute/state
- Eval harness: https://github.com/ArcInstitute/cell-eval
- DE backend: https://github.com/ArcInstitute/pdex
- Dataloaders: https://github.com/ArcInstitute/cell-load
- Preprint: https://www.biorxiv.org/content/10.1101/2025.06.26.661135v2
- Manuscript landing: https://arcinstitute.org/manuscripts/State
- Inference Colab on Tahoe-100M ST: https://colab.research.google.com/drive/1bq5v7hixnM-tZHwNdgPiuuDo6kuiwLKJ
- SE embedding Colab: https://colab.research.google.com/drive/1uJinTJLSesJeot0mP254fQpSxGuDEsZt
- VCC training Colab: https://colab.research.google.com/drive/1QKOtYP7bMpdgDJEipDxaJqOchv7oQ-_l

## Verbatim quotes worth keeping

From the v1 abstract (bioRxiv 2025.06.26.661135): "Here, we introduce STATE, a transformer model that predicts perturbation effects while accounting for cellular heterogeneity within and across experiments. STATE predicts perturbation effects across sets of cells and is trained using gene expression data from over 100 million perturbed cells. STATE improved discrimination of effects on large datasets by more than 30% and identified differentially expressed genes across genetic, signaling and chemical perturbations with significantly improved accuracy. Using its cell embedding trained on observational data from 167 million cells, STATE identified strong perturbations in novel cellular contexts where no perturbations were observed during training. We further introduce Cell-Eval, a comprehensive evaluation framework that highlights STATE's ability to detect cell type-specific perturbation responses, such as cell survival."

From state_transition.py docstring: "This model: 1) Projects basal expression and perturbation encodings into a shared latent space. 2) Uses an OT-based distributional loss (energy, sinkhorn, etc.) from geomloss. 3) Enables cells to attend to one another, learning a set-to-set function rather than a sample-to-sample single-cell map."

From cell-eval README: "This package provides a comprehensive suite of metrics for evaluating the performance of models that predict cellular responses to perturbations at the single-cell level."

ConfPert context note: STATE is the most recent canonical scFM in the perturbation lineage and ships with Cell-Eval. Cell-Eval is the ConfPert evaluation backbone (per user spec: do not reinvent). STATE is wrappable predictor #8.
