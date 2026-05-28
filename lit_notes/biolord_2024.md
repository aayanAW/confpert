# Piran et al. 2024 - biolord

**Venue:** Nature Biotechnology 2024. Piran Z, Cohen N, Hoshen Y, Nitzan M. "Disentanglement of single-cell data with biolord." Nat Biotechnol 42(11):1678-1683. doi:10.1038/s41587-023-02079-x. PMID 38225466. PMCID PMC11554562. Epub 2024-01-15.

## Claim
biolord is a deep generative framework that performs hard disentanglement of single-cell expression into known attributes (categorical or ordered: drug, dose, time, age, tissue, infection state) and unknown attributes (residual cell state). Each known attribute gets a dedicated encoder/embedding subnetwork; unknown attributes are stored as per-sample latent codes regularized toward zero. The decoder reconstructs expression from the concatenated latent. By holding the unknown-attribute embedding fixed and swapping known-attribute embeddings to target values, biolord produces counterfactual single-cell expression for unseen drug, dose, perturbation, time, or age combinations. Headline result: predicts cellular response to unseen drugs and unseen genetic perturbations, beating chemCPA-pre and GEARS on sci-Plex 3 and Norman 2019 Perturb-seq.

## Method
**Architecture.** Per-attribute encoder subnetworks. Categorical attributes (e.g., drug identity, cell line, donor) get shared embedding tables optimized via latent optimization. Ordered/continuous attributes (e.g., dose, age, time) get small MLPs (depth 2, width 256) mapping the scalar feature to a latent code. Unknown attributes are stored as per-cell latent vectors `z_u` with Gaussian noise regularization. The decoder G takes the concatenation of all known-attribute embeddings plus `z_u` and outputs reconstructed expression (NB or Gaussian likelihood).

**Hard disentanglement via constrained optimization.** Loss combines a completeness term `L_cmp` (NLL/MSE reconstruction) and a minimality term `L_min = λ ||z_u||^2` constraining the unknown-attribute latent. Quoted: "The model's loss function attempts to maximize the accuracy of the reconstruction (enforcing completeness) while minimizing the information encoded in the unknown attributes (limiting its capacity)." This is not adversarial in the GAN sense; the disentanglement is enforced by *latent optimization* on `z_u` plus capacity constraint, so that any information explained by known attributes flows through their dedicated encoders rather than leaking into `z_u`.

**Inductive bias.** Quoted: "Randomly initialized latent codes trivially do not contain any information on known or unknown attributes" whereas encoder networks initialized randomly start from a "perfectly entangled state," making subsequent disentanglement harder. This is the central architectural argument vs CPA.

**Difference vs CPA.** CPA uses a single autoencoder with adversarial classifiers removing covariate signal from a bottleneck `z_basal` and then adds back learned drug+dose embeddings. biolord (i) replaces the bottleneck encoder with per-cell latent optimization, (ii) replaces adversarial covariate removal with a capacity penalty plus dedicated per-attribute encoders, (iii) supports arbitrary mixes of categorical and ordered attributes uniformly. Result: more flexible covariate space, better OOD on multi-attribute regimes (drug + dose + cell line).

**Counterfactual prediction.** Quoted: "Given a control cell and unseen (target) labels as input, biolord can predict the gene expression of the unseen cellular states." Hold `z_u` of a reference cell fixed, swap the target-attribute embedding to the unseen value, decode. The `predict()` API returns AnnData with predicted "expression mean and variance"; per-cell samples are produced by passing the corresponding control AnnData through `predict()` so each control cell yields one counterfactual cell.

## Datasets / experiments
- **sci-Plex 3 (Srivatsan 2020):** 649,340 cells, 3 cancer cell lines, 188 compounds, 4 dosages. OOD task: predict response to held-out drugs at held-out doses.
- **Norman 2019 Perturb-seq (GI dataset):** ~89,357 cells, 131 two-gene + 105 one-gene CRISPRa perturbations. OOD splits: zero/one/two genes unseen.
- **Perturb-seq one-gene set:** 65,899 cells, 81 perturbations.
- **Plasmodium liver infection atlas:** spatiotemporal scRNA-seq, 5 infection timepoints + control. Used for infection-state counterfactuals and biolord-classify semi-supervised extension of an abortive hepatocyte population from 36 hpi back to 24 and 30 hpi.
- **Fetal chromatin atlas:** multi-omic disentanglement of tissue, age, cell-type.

Concrete numbers (best-of-class shown):
- sci-Plex 3 OOD r^2 at 10 uM dose: **biolord 0.76 +/- 0.0005 vs chemCPA-pre 0.51 +/- 0.0062**.
- sci-Plex 3 with 10% subsampling: biolord r^2 = 0.63 +/- 0.0003 (robustness).
- Drug-embedding uncertainty (lower better): biolord 0.19 +/- 0.005 vs RDKit 0.32 +/- 0.008.
- Norman two-gene Perturb-seq normalized MSE (biolord vs GEARS): two unseen 0.50 vs 0.53; one unseen 0.35 vs 0.39; zero unseen 0.20 vs 0.28.
- One-gene unseen Perturb-seq normalized MSE: biolord 0.37 vs GEARS 0.47.

Note: aging and IFN PBMC variants exist in companion analyses / preprint v1; the published Nature Biotech version foregrounds sci-Plex 3, Perturb-seq one/two-gene, and the Plasmodium atlas.

## Metrics
- r^2 on per-perturbation mean expression (primary).
- r^2 on top-100 differentially expressed genes (DE-gene focused).
- Normalized MSE for genetic perturbation OOD (Norman).
- Drug-embedding uncertainty via k-NN pathway-prediction on the latent graph.
- Retrieval accuracy via `evaluate_retrieval`.
- Reported as mean +/- SEM across seeds.

## What we steal for ConfPert
biolord is **wrappable predictor #6**. Like CPA and scGen it produces samples at predict time via the learned attribute decoder applied to a control input cell, but it is materially stronger on multi-attribute / multi-dose / unseen-combination regimes because of the hard disentanglement and per-attribute encoder architecture. ConfPert wraps biolord output samples with conformal calibration:
- For each control cell c and target perturbation t, draw `m` samples by perturbing `z_u(c)` with small Gaussian jitter and decoding under target attribute embedding e(t). This gives the per-(c,t) sample cloud needed for split-conformal radii on per-cell expression vectors or DE-gene subsets.
- Use `predict(adata, indices, nullify_attribute)` to obtain mean and variance heads; combine with sampled `z_u` perturbations to reconstruct a calibrated predictive distribution.
- biolord's strength on Norman two-gene OOD is exactly the K1 (combination) regime where ConfPert needs tight CP intervals.

## What we wrap
biolord = wrappable predictor #6. github.com/nitzanlab/biolord. Reproducibility scripts at github.com/nitzanlab/biolord_reproducibility. Built on scvi-tools, so the `Biolord` class subclasses scvi `BaseModelClass`. Public model weights are not advertised; the reproducibility repo provides training scripts for sci-Plex 3 and Perturb-seq that we can rerun to obtain checkpoints. API: `Biolord.setup_anndata(adata, ordered_attributes_keys=..., categorical_attributes_keys=..., categorical_attributes_missing=..., retrieval_attribute_key=..., layer=...)`, then `model.train(...)`, then `model.predict(adata=control_adata, indices=..., batch_size=512, nullify_attribute=None)` returning two AnnData (mean and variance). Per-cell counterfactuals: pass control cells through `predict()` after rewriting the attribute keys in `adata.obs` to target values, or use `compute_prediction_adata(adata, adata_source, target_attributes, add_attributes=None)` which generates predictions across attribute combinations. Latent for unknown-state sampling: `get_latent_representation_adata(adata, indices, batch_size, nullify_attribute)`.

## Failure modes / caveats
- Hard-disentanglement assumption strict; correlated attributes degrade OOD predictions. Quoted: "It is unclear what is the desired outcome when attributes are correlated... predictions over unseen combinations may yield unpredictable results, which is a known limitation of neural networks." Direct hit on K1 multi-perturbation combos.
- Lack of interpretability: "As with any deep generative framework, biolord suffers from the lack of direct interpretability. We overcome this by suggesting various downstream analysis tools." Conformal wrapping does not fix this but does add coverage guarantees.
- Latent optimization is per-sample, so adding a brand-new cell at test time requires either treating it as a control (its `z_u` is implicit in the decoded reference) or running a short fitting step. ConfPert calibration set must be constructed from cells the model has indexed.
- No native generative likelihood for unseen donors/cell-lines outside the training categorical vocabulary; biolord cannot embed an unseen categorical level without retraining or a similarity-based proxy. ConfPert under K3 PRISM should ensure the drug vocabulary at calibration time is a subset of training.
- Multi-attribute training stability: dependence on `lambda` (capacity penalty) and on choice of which attributes are "known"; biolord-classify is needed when labels are partial.
- Reported variances are very small (SEM ~1e-4); replicate seeding before trusting deltas.

## Code URLs
- https://github.com/nitzanlab/biolord
- https://biolord.readthedocs.io
- https://github.com/nitzanlab/biolord_reproducibility

## Verbatim quotes worth keeping
1. "Biolord is a deep generative method for disentangling single-cell multi-omic data to known and unknown attributes, including spatial, temporal and disease states."
2. "The model's loss function attempts to maximize the accuracy of the reconstruction (enforcing completeness) while minimizing the information encoded in the unknown attributes (limiting its capacity)."
3. "Randomly initialized latent codes trivially do not contain any information on known or unknown attributes."
4. "Given a control cell and unseen (target) labels as input, biolord can predict the gene expression of the unseen cellular states."
5. "It is unclear what is the desired outcome when attributes are correlated... predictions over unseen combinations may yield unpredictable results, which is a known limitation of neural networks."
6. "As with any deep generative framework, biolord suffers from the lack of direct interpretability. We overcome this by suggesting various downstream analysis tools."

---
ConfPert context: biolord is the most attribute-flexible wrappable predictor in the stack. K3 PRISM drug-resistance is a natural fit because PRISM has explicit drug + dose + cell-line attributes that map cleanly to biolord's categorical+ordered API. `predict()` returns mean and variance AnnData; per-cell counterfactual samples come from passing each control cell through `predict()` after rewriting `adata.obs` to target attribute values. Confirmed predict() API is per-cell (one output cell per input control cell), so split-conformal calibration on per-cell residuals or per-pert pseudobulk residuals is direct.
