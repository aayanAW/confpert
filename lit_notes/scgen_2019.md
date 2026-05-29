# Lotfollahi, Wolf, Theis 2019 - scGen

**Venue:** Nature Methods 2019. DOI: 10.1038/s41592-019-0494-8. Authors: Mohammad Lotfollahi, F. Alexander Wolf, Fabian J. Theis (Helmholtz Munich / TUM).

## Claim

scGen uses a variational autoencoder (VAE) plus latent vector arithmetic to predict perturbed single-cell expression from controls. The core operation is delta_p = mean(z_perturbed) - mean(z_control) in latent space; predicted perturbed cells are produced by encoding held-out controls, adding delta_p, and decoding back to gene space. The model generalizes across cell types, studies, and species, predicting perturbation responses for cell types observed in only the control condition.

## Method

**Architecture.** Standard Gaussian VAE. Defaults from the public implementation: `n_hidden=800`, `n_latent=100`, `n_layers=2`, `dropout_rate=0.2`. Encoder outputs latent mean qz_m and latent variance qz_v; decoder is a feedforward network with linear output head producing the reconstructed gene-expression vector. Input is log1p-normalized expression (typically restricted to ~7000 highly variable genes per the README guidance).

**Loss (verbatim from `_scgenvae.py`):**
```
kld = kl(Normal(qz_m, torch.sqrt(qz_v)), Normal(0, 1)).sum(dim=1)
rl  = ((x - px) ** 2).sum(dim=1)
loss = (0.5 * rl + 0.5 * (kld * self.kl_weight)).mean()
```
Reconstruction is squared error (Gaussian likelihood, identity link), KL is to a standard normal prior, weights are 0.5/0.5 with a tunable `kl_weight`.

**Prediction procedure.**
1. Train VAE on the joint dataset (control + perturbed cells, all observed cell types).
2. Compute condition centroids in latent space: z_ctrl_bar = mean over latent codes of control cells; z_stim_bar = mean over latent codes of perturbed cells (computed within matched cell types when available).
3. Compute the perturbation vector: `delta = z_stim_bar - z_ctrl_bar`.
4. For a held-out cell type observed only under control: encode each control cell c -> z_c, add delta -> z_c + delta, decode -> predicted perturbed expression for cell c.

**Per-cell sampling behavior.** The public `predict()` API returns an AnnData of per-cell predictions, not a single mean vector. From the source: `stim_pred = delta + latent_cd` where `latent_cd` is the per-cell control encoding. Heterogeneity comes from (a) the per-cell control input variance (each control cell has a distinct latent code), and (b) optionally sampling z ~ N(qz_m, qz_v) at encode time and/or sampling from the decoder's output distribution. The predict() default is deterministic on the encoder mean; stochastic samples are accessible via the `sample()` method on the underlying VAE.

**API signature (public repo, `scgen/_scgen.py`):**
```
def predict(
    self,
    ctrl_key=None,
    stim_key=None,
    adata_to_predict=None,
    celltype_to_predict=None,
    restrict_arithmetic_to="all",
) -> AnnData
```
Returns `(predicted_adata, delta)` where `predicted_adata` is one predicted cell per input control cell.

## Datasets / experiments

- **IFN-beta stimulation (Kang et al., 2018) PBMCs.** The headline benchmark: train on all cell types in control + a subset of cell types in IFN-beta condition; hold out one cell type (e.g., CD4 T cells) entirely from the stimulated set; predict its IFN-beta response. Compared against alternative approaches (CVAE, scVI baseline, linear models).
- **Intestinal epithelium (Haber et al., 2017).** Salmonella and H. polygyrus infection responses across enterocytes, goblet cells, etc.
- **Hematopoiesis / cross-study generalization.** Demonstrating transfer across studies.
- **Cross-species prediction.** Mouse <-> human / rat for LPS or related stimuli, demonstrating that the latent arithmetic transfers across species when conditioning on shared cell-type identity.
- **Out-of-sample cell types.** Key validation: a cell type whose perturbed state was never seen in training is predicted by applying delta learned from other cell types.

## Metrics

- **R^2 on mean expression.** Pearson R^2 between predicted mean expression vector (averaged over predicted cells of a held-out type) and observed real perturbed mean. Reported both on all genes and on the top differentially expressed (DE) genes. Typical headline numbers in the IFN-beta task: R^2 around 0.94-0.98 on all genes, lower on top-100 DE genes.
- **R^2 on variance.** Same regression but on per-gene variance across cells, capturing whether scGen reproduces cell-to-cell heterogeneity rather than only the population mean.
- **Top DE gene overlap / DE recovery.** Differential expression run on real perturbed vs control and on predicted perturbed vs control; overlap of top-K (typically top 100) DE genes is reported.
- **Qualitative.** UMAP / PCA embedding of real-stim, control, and predicted-stim cells; violin plots for marker genes (e.g., ISG15) showing distribution match.
- **Built-in plotting.** `model.reg_mean_plot()` and `model.reg_var_plot()` in the released package implement the mean and variance R^2 evaluations directly.

## What we steal for ConfPert

scGen is the simplest wrappable predictor: a clean encode -> add delta -> decode pipeline with a stable `predict()` API that returns per-cell AnnData, which is exactly the input shape ConfPert's calibration step expects. Useful as a strong-and-simple baseline in K1's nine-predictor sweep. The latent arithmetic gives a clean handle on per-cell heterogeneity: control input variance is preserved through the encode-decode, so the marginal predictive distribution over a held-out cell type is non-degenerate without extra effort. The released `reg_mean_plot` and `reg_var_plot` give us a free sanity-check on whether the wrapped predictor's mean/variance behavior matches the published scGen numbers before we drop it into the conformal pipeline.

## What we wrap

scGen is wrappable predictor #2. Repo: github.com/theislab/scgen. Built on scvi-tools (`scgen.SCGEN(adata)` -> `model.train()` -> `model.predict(ctrl_key, stim_key, celltype_to_predict)`). Public weights are not distributed as separate checkpoints; the standard practice (and what the tutorials assume) is to retrain on the target dataset for a few hundred epochs (`max_epochs=100` with early stopping in the official tutorial). For ConfPert this means scGen is wrapped as a re-trainable predictor per dataset, not a frozen pretrained model.

## Failure modes / caveats

- **Additive perturbation in latent space.** delta_p is a single mean vector; the method assumes the perturbation effect is approximately additive in latent space and roughly cell-type-invariant (or at least transferable). Combinatorial perturbations (drug A + drug B), dose-response, and strongly non-linear interactions are not captured by a single delta.
- **No explicit conditioning.** Unlike CVAE, scGen does not condition the encoder/decoder on perturbation labels; the latent space must spontaneously align conditions for arithmetic to work. Empirically this works on the benchmarked datasets but is fragile when batch effects dominate biological signal.
- **Single-condition VAE limits OOD generalization.** Distribution shifts at the cell-type level (e.g., a held-out type with very different baseline transcriptional state) push delta arithmetic out of the regime where it was estimated.
- **Variance prediction is weaker than mean.** Reported R^2 on variance is consistently lower than R^2 on mean, meaning the per-cell heterogeneity match is the weakest part of the model and is exactly where conformal calibration adds value.
- **Gaussian likelihood on log-normalized counts.** The squared-error reconstruction is convenient but mis-specified for raw counts; users must log1p-normalize first.
- **delta is computed from population means**, so cell-type-specific perturbation effects (effect modification by cell state) are averaged away unless `restrict_arithmetic_to` is used to compute delta within strata.

## Code URLs

- https://github.com/theislab/scgen
- https://scgen.readthedocs.io/en/stable/
- https://scgen.readthedocs.io/en/stable/tutorials/scgen_perturbation_prediction.html
- Paper: https://www.nature.com/articles/s41592-019-0494-8

## Verbatim quotes worth keeping

From the public repo and docs (load-bearing for the wrapper API):

1. `predict()` docstring: "Predicts the cell type provided by the user in stimulated condition."
2. README framing: scGen is "a generative model to predict single-cell perturbation response across cell types, studies and species."
3. README usage scope: "Train on a dataset with multiple cell types and conditions and predict the perturbation effect on the cell type which you only have in one condition."
4. `predict()` returns: a "numpy matrix of our predicted cells and the second one is the difference vector between our conditions" (delta).
5. Loss (verbatim): `loss = (0.5 * rl + 0.5 * (kld * self.kl_weight)).mean()` with `rl = ((x - px) ** 2).sum(dim=1)` and `kld = kl(Normal(qz_m, sqrt(qz_v)), Normal(0,1)).sum(dim=1)`.

## ConfPert integration note

Confirmed: `predict()` returns per-cell predictions (one predicted perturbed cell per input control cell), packaged as an AnnData. This is the right shape for ConfPert's K1 calibration step, which needs a sample-producing predictor rather than a mean-only predictor. For K1 we treat scGen as the lightweight reference baseline against scVI / scFoundation / GEARS / CPA in the nine-predictor sweep. Per-cell heterogeneity comes from the diversity of control inputs; if we want extra stochasticity we can additionally sample z ~ N(qz_m, qz_v) at encode time using the `sample()` method on the underlying VAE module.
