# Cui et al. 2024 - scGPT

**Venue:** Nature Methods 2024 (published doi 10.1038/s41592-024-02201-0; bioRxiv v2 2023.04.30.538439, posted 2 Jul 2023, methods identical to published version).

## Claim

scGPT is a 12-layer transformer pretrained on 33M normal human scRNA-seq cells from CELLxGENE (51 organs/tissues, 441 studies) via masked iterative gene-expression generation with MSE on per-cell binned expression values (eq. 10). Fine-tuned variants address cell-type annotation, batch integration, scMultiomic integration, GRN inference, and perturbation prediction (Perturb-GEP). The paper reports highest Pearson and Pearson-delta on top-20 DE genes against GEARS and CPA on Adamson and Norman.

## Method

**Three-token input (Methods 4.1).** Each cell is a length-M sequence of (gene-token, expression-value-token, condition-token). Gene tokens are integer ids over the human gene vocabulary plus `<cls>` and `<pad>`. Expression values are binned per-cell into B equal-mass bins of non-zero counts ("a new set of bin edges is computed for each cell"). Condition tokens are position-wise integers per gene (eq. 4) used for, among other things, perturbation indicators. Embedding (eq. 5): h^(i) = emb_g(t_g) + emb_x(x) + emb_c(t_c).

**Architecture (Methods 4.8).** "The pre-trained foundation model has an embedding size of 512. It consists of 12 stacked transformer blocks with 8 attention heads each. The fully connected layer has a hidden size of 512." Max input length M = 1200 (only non-zero genes input at pretrain). FlashAttention used: "we leverage the accelerated self-attention implementation by Flash-Attention." Pretrain optimizer: Adam, lr 1e-4, weight decay 0.9 per epoch, batch 32, 6 epochs. 99.7% / 0.3% train/val split on 33M whole-human cells.

**Pretraining loss (eq. 10, Methods 4.3.3).** A specialized non-causal attention mask supports both gene-prompt (predict masked gene values from observed) and cell-prompt (predict whole-genome from `<cls>` cell-type embedding) generation in iterative steps (top-1/K confidence selection). Loss is MSE between an MLP head on the final representation and the integer bin index target at unknown (masked) positions:

> L = (1/|U_unk|) sum_{j in U_unk} (MLP(h_n^(i)) - x_j^(i))^2.   (eq. 10)

x_j is the binned expression value (integer 0..B). Gene-prompt and cell-prompt losses are summed within each step.

**Fine-tune objectives (Methods 4.4).** GEP repeats eq. 10 at masked positions over the fine-tune data:

> tilde-x^(i) = MLP(h_n^(i)),  L_GEP = (1/|M_mask|) sum_{j in M_mask} (tilde-x_j^(i) - x_j^(i))^2.   (eq. 11)

Perturb-GEP is defined as a variant of GEP: "We refer to this variation as perturb-GEP. We maintain the MLP estimator in equation 11, but utilize the post-perturbation gene expression as the target x_j^(i). In perturb-GEP, the model is supposed to predict the post-perturbation expression of all input genes." There is one scalar MLP head per gene; no variance head, no mixture, no NB likelihood at perturb fine-tune time.

**Perturbation pipeline (Methods 4.5, "Perturbation response prediction").** Two changes from the rest of the framework, quoted verbatim:

> "First, we used log1p transformed expression values as input and target values instead of binned values, to better predict the absolute post-perturbation expression for this task. Second, we appended a binary condition token at each input gene position to indicate whether the gene is perturbed."

Training pairs are (control_cell, perturbed_cell): "Instead of utilizing the masked and unmasked expression values of the same cell as the input and learning target, we employed a control cell as the input and the perturbed cell as the target. This was achieved by randomly pairing a non-perturbed control cell with each perturbed cell to construct the input-target pairs. The input values consisted of all the non-perturbed gene expression values in the control cells. Consequently, the model learned to predict the post-perturbation responses based on the control gene expression."

Pre-processing: HVG = 5000 plus union of any perturbed gene; for one-gene perts the test split has zero overlap with training perts; for Norman two-gene the splits are 0/2, 1/2, 2/2 unseen in training.

**Inference (Methods S.5 perturbation).** A single mean per perturbation is produced, not samples:

> "scGPT predicted the representative perturbation response for each perturbation condition from a single vector of sample control gene expression (i.e., of size 1 x M genes), obtained by averaging the gene expression of all control cells in the dataset."

So per-cell predictive variance under scGPT-Perturb is zero by architecture: MSE-trained scalar head + deterministic decoder + averaged control input collapse the per-condition output to one mean vector.

## Datasets / experiments

- **Pretrain.** 33M human cells from CELLxGENE Census (May 15 2023 release), 51 organs / 441 studies; Brain 13.2M, Blood 10.3M, Lung 2.1M, Heart 1.8M, Pancreas 210k, Kidney 814k, Intestine 94.5k.
- **Adamson Perturb-seq, K562 CRISPRi:** 87 one-gene perturbations.
- **Norman Perturb-seq, K562 CRISPRa:** 105 one-gene + 131 two-gene = 236 perturbations across 105 unique target genes; 5,565 possible 1-gene + 2-gene combinations of those 105 targets.
- Splits and metrics follow GEARS.

## Metrics

- **corr (Pearson) on top-20 DE genes per pert.** "We reported these Pearson metrics on the top 20 genes of most changed expression (DE genes)." DE selected via the dataset-defined DE list per perturbation (GEARS convention).
- **corr(delta).** Pearson on the change vector (predicted - control) vs (true - control), top-20 DE.
- Both are mean-only metrics; neither penalizes a degenerate point predictor that ignores cell-to-cell heterogeneity.
- Headline claim: scGPT achieves the highest correlation on both datasets across both metrics (paper Fig. 3A; "highest correlation 7/8" framing comes from the scGPT vs GEARS vs CPA bar chart across {Adamson, Norman 0/2, 1/2, 2/2} x {corr, corr-delta}).

## What we steal for ConfPert

scGPT is wrappable predictor #5. `predict(pert, control_cell)` returns one log1p mean expression vector. Under the paper's evaluation protocol it returns one mean per pert (averaged control). Like GEARS it is mean-only: per-cell variance = 0 by architecture (scalar MLP head, MSE training, deterministic forward). Clean K1 demonstrator: even a 33M-cell-pretrained foundation model gives zero distributional coverage out of the box. Stage-2 HetPert claim generalizes: MSE + per-pert conditioning forces variance recovery (VR) -> 0 regardless of pretraining scale. ConfPert layers calibrated coverage on top of any sample-augmented variant (dropout MC, residual quantile head).

## What we wrap

Two paths. (1) **Architecture probe.** Re-implement the gene-token + binary-perturbation-token + log1p transformer per Methods 4.5. Sufficient for K1 and the HetPert variance-collapse claim. (2) **Published checkpoint.** github.com/bowang-lab/scGPT, whole-human checkpoint plus Adamson/Norman fine-tune scripts; HuggingFace mirror. ~$10/run on Modal H100.

## Failure modes / caveats

- Architecturally precluded from variance recovery at perturb fine-tune time. No NB, no Gaussian variance head, no mixture, no flow.
- Internal pretrain/fine-tune inconsistency: pretraining target is integer per-cell bin index (eq. 10), perturb fine-tune target is log1p continuous expression (Methods 4.5). MSE on bin indices and MSE on log1p are not the same loss and the embedding layer for x is reused across both regimes; the paper does not characterize what this does to gene-token semantics.
- `scgpt` pip package is fragile: hard pin on `flash-attn` build, transitive `IPython` dep often breaks installs, CUDA version drift; expect to vendor and patch.
- Inference protocol uses one averaged control vector per dataset, not per-cell; this further removes any path to a sample-level distribution from the published pipeline.
- Evaluation is mean-only (corr, corr-delta on top-20 DE). Recent benchmarks (Ahlmann-Eltze 2025; PerturBench 2024) show these metrics are saturated by linear or "predict-control" baselines, so the 7/8 leadership claim does not transfer to distributional or causal evaluations.

## Code URLs

- https://github.com/bowang-lab/scGPT
- HuggingFace mirror for whole-human, brain, blood, heart, lung, kidney, pan-cancer checkpoints (linked from the README).
- Tutorial: tutorials/Tutorial_Perturbation.ipynb in the repo.

## Verbatim quotes worth keeping

1. Pretraining loss, eq. 10: "L = 1/|U_unk| sum_{j in U_unk} (MLP(h_n^(i)) - x_j^(i))^2 ... where U_unk denotes the set of the output positions for unknown genes, and the x_j^(i) is the actual gene expression value to be predicted."
2. Fine-tune loss, eq. 11: "tilde-x^(i) = MLP(h_n^(i)), L_GEP = 1/|M_mask| sum_{j in M_mask} (tilde-x_j^(i) - x_j^(i))^2."
3. Perturb-GEP definition: "We refer to this variation as perturb-GEP. We maintain the MLP estimator in equation 11, but utilize the post-perturbation gene expression as the target x_j^(i). In perturb-GEP, the model is supposed to predict the post-perturbation expression of all input genes."
4. Perturb input changes: "we used log1p transformed expression values as input and target values instead of binned values ... we appended a binary condition token at each input gene position to indicate whether the gene is perturbed."
5. Inference protocol: "scGPT predicted the representative perturbation response for each perturbation condition from a single vector of sample control gene expression (i.e., of size 1 x M genes), obtained by averaging the gene expression of all control cells in the dataset."
6. Architecture: "The pre-trained foundation model has an embedding size of 512. It consists of 12 stacked transformer blocks with 8 attention heads each."

## ConfPert role

Predictor #5 of 9. Demonstrates K1 (foundation-model mean-only -> zero distributional coverage out of the box). Confirms the Stage-2 HetPert claim that MSE + per-pert conditioning forces VR -> 0 independent of capacity or pretraining scale.
