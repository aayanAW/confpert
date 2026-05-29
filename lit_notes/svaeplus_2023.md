# Bereket and Karaletsos 2023 - sVAE+ / SAMS-VAE

**Venue:** NeurIPS 2023. arXiv 2311.02794. Authors: Michael Bereket, Theofanis Karaletsos (Insitro).

## Claim

SAMS-VAE (sVAE+ in some figures) is a generative latent variable model for single-cell perturbation response that decomposes a perturbed cell's latent state into a sample-specific basal latent and a sparse, additive, perturbation-specific mechanism shift. Verbatim from the paper: "SAMS-VAE models the latent state of a perturbed sample as the sum of a local latent variable capturing sample-specific variation and sparse global variables of latent intervention effects." The model is pitched as combining "compositionality, disentanglement, and interpretability for perturbation models." Output is a full posterior over perturbed cell states, so it is a sample-producing distributional predictor, not a point estimator.

## Method

Latent decomposition. Each perturbed cell's latent z is the sum z = z_basal + sum_p (e_p \* mask_p) where z_basal is a per-cell local latent (analogous to CPA's basal state), e_p is a global per-perturbation mechanism-shift latent, and mask_p is a sparsity gate (Bernoulli or relaxed Bernoulli) over latent dimensions. The shift is added only along the dimensions selected by mask_p, which is the "sparse additive mechanism shift" in the title. Perturbation effects compose by summing the masked shifts of the active perturbations, giving combinatorial generalization with no extra parameters per pair.

Generative model. p(x | z) is a count-likelihood decoder (negative binomial in the released code) over genes; p(z_basal) is a standard Gaussian; e_p has a Gaussian prior; mask_p has a Beta-Bernoulli or concrete-relaxation sparsity prior with a tunable sparsity hyperparameter. Inference is amortised variational with an encoder q(z_basal | x) and global variational posteriors over (e_p, mask_p).

ELBO. Standard mean-field-style ELBO with three terms: reconstruction log p(x | z), KL on the local basal latent, and KL on the global mechanism-shift latents and their masks (the latter penalises non-sparse masks via the Beta-Bernoulli / concrete prior).

Predict-time. To predict a perturbed expression for a held-out perturbation set P, sample z_basal ~ q(z_basal | x_control), add sum_{p in P} e_p \* mask_p drawn from the global posterior, and decode through p(x | z). Because both the basal posterior and the global mechanism-shift posterior are stochastic, draws produce a full predictive distribution per (cell, perturbation set), not just a mean.

Identifiability. The sparse-additive structure plus the per-perturbation mechanism prior is what gives SAMS-VAE its claimed disentanglement and compositional generalisation; without sparsity the additive shift collapses into entangled global directions. The paper introduces "a framework for evaluation of perturbation models based on average treatment effects with links to posterior predictive checks" to assess this.

## Datasets / experiments

- Norman et al. 2019 Perturb-seq (CRISPRa, K562). The repo distributes a preprocessed Norman dataset (~1.6 GB). Used for combinatorial generalisation: held-out double-gene perturbations whose singletons are seen in training.
- Replogle et al. 2022 genome-scale Perturb-seq (essential genes, K562/RPE1). Repo distributes a preprocessed Replogle dataset (~550 MB). Used for held-out single perturbations and large-scale fits.
- sci-Plex (Srivatsan et al. 2020) is referenced as a small-molecule benchmark consistent with the CPA evaluation protocol.

The paper benchmarks against CPA and conditional VAE baselines and reports gains on in-distribution reconstruction, OOD perturbation prediction, and combinatorial generalisation, plus interpretability of the inferred mechanism-shift latents which "correlate strongly to known biological mechanisms."

## Metrics

- R^2 on predicted vs observed mean expression across all genes for each held-out perturbation.
- R^2 on top-DE genes (the same top-K DEG protocol used by CPA / scGen).
- Average treatment effect (ATE) framework introduced in this paper, with posterior predictive checks: predicted vs observed condition-level ATE on held-out perturbations.
- MMD-style discrepancy between predicted and observed perturbed-cell distributions for combinatorial / OOD splits.
- Sparsity statistics on the inferred mask_p (fraction of active latent dimensions per perturbation) as an interpretability metric.

## What we steal for ConfPert

SAMS-VAE is a clean fit for ConfPert as wrappable predictor #7. Three properties matter:

1. It is a genuine sample-producer. The full variational posterior over (z_basal, e_p, mask_p) plus the count-likelihood decoder gives Monte Carlo draws per (cell, perturbation set), which is exactly what ConfPert's conformal layer consumes. No need to retrofit a sampler the way we do for CPA's mean head.
2. The sparse-additive prior gives compositional OOD generalisation on combinations, which is the regime where conformal coverage under shift is most informative. ConfPert can report calibrated intervals on Norman double-gene held-outs where SAMS-VAE is already strong on point R^2.
3. The ATE-based evaluation in the paper aligns with ConfPert's distributional targets, so we can directly compare conformal coverage against the paper's posterior-predictive-check protocol.

For K1's predictor sweep we wrap SAMS-VAE alongside CPA, scGen, CellOT, and the count-based baselines, using identical Norman / Replogle splits.

## What we wrap

SAMS-VAE is wrappable predictor #7. The Insitro repo at https://github.com/insitro/sams-vae is the canonical implementation. Install via conda lockfile + `pip install -e .` (Linux env was used for paper results; macOS arm64 lock is also provided). Datasets are pulled by `python download_datasets.py --replogle --norman` and cached in `$SAMS_VAE_DATASET_DIR`. Training is config-driven: `python train.py --config tests/models/sams_vae_correlated.yaml`. Sweeps run via redun + wandb. Precomputed metrics and checkpoints are hosted on a public S3 bucket (AWS CLI required). Repo is mostly Jupyter notebooks (~98%) for paper replication, with the model code in `sams_vae/` and replication scripts under `paper/experiments/`.

## Failure modes / caveats

- Sparsity prior tuning. The Beta-Bernoulli / concrete sparsity hyperparameter (alpha-equivalent) directly controls how many latent dimensions a perturbation can touch. Too tight and the model underfits strong perturbations; too loose and disentanglement collapses to a CPA-like dense additive embedding.
- Posterior collapse on z_basal. Standard VAE failure mode: when the decoder is too expressive relative to the basal latent, q(z_basal | x) collapses to the prior and all variation is absorbed into the perturbation shift, breaking the cell-conditional counterfactual interpretation.
- Mechanism-shift identifiability without ground-truth pathway labels. The model recovers shift directions that are sparse and additive but the mapping from latent dimensions to biological pathways is post hoc and not guaranteed unique; rotation/permutation ambiguity persists across runs.
- Combinatorial generalisation is genuinely OOD only when the singletons are seen; predictions for double-gene sets where one constituent is unseen revert to weak priors.
- Count likelihood (NB) requires correct size-factor handling; mis-set library-size offsets bias mean predictions on lowly expressed genes, same caveat as CPA's NB head.
- No native calibration guarantee on the predictive distribution. The variational posterior is approximate; nominal credible intervals from posterior draws are not coverage-calibrated. This is exactly the gap ConfPert closes.

## Code URLs

- GitHub: https://github.com/insitro/sams-vae
- arXiv: https://arxiv.org/abs/2311.02794
- NeurIPS 2023 proceedings entry (search by title; OpenReview ID redirects from the proceedings hash page).

## Verbatim quotes worth keeping

1. "SAMS-VAE models the latent state of a perturbed sample as the sum of a local latent variable capturing sample-specific variation and sparse global variables of latent intervention effects."
2. "compositionality, disentanglement, and interpretability for perturbation models."
3. "a framework for evaluation of perturbation models based on average treatment effects with links to posterior predictive checks."
4. "Interpretable latent structures which correlate strongly to known biological mechanisms."
5. "The results in the paper were generated using the Linux environment." (repo README, relevant for reproducibility on our cluster.)

## Predict() output format - confirmed

The forward pass returns a full count-likelihood distribution per (cell, perturbation set). For ConfPert we draw S samples by jointly sampling (z_basal_s ~ q(z_basal | x_control), e_p_s, mask_p_s ~ q_global) and decoding to NB(mu_s, theta_s), then drawing one count vector per s. This gives S full transcriptome samples per query, which feed straight into the conformal calibration layer alongside CPA and scGen samples on identical Norman / Replogle splits.
