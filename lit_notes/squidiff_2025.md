# He et al. 2025 - Squidiff

**Venue:** Nature Methods 2025 (article s41592-025-02877-y, vol. 23, pp. 65-77, 2026 issue). Preprint: bioRxiv 10.1101/2024.11.16.623974, Nov 2024.

## Claim

Squidiff is a conditional latent diffusion model for single-cell transcriptomic prediction. It predicts transcriptomic responses across cell types under environmental perturbations (cell differentiation, gene knockout/knockdown, drug response, irradiation, growth factor stimulation). The acronym stands for Single-cell QUantitative Inference of stimuli responses by a DIFFusion model. Conditioning is via a learned semantic latent variable z_sem that encodes cell type, environmental change, disease state, and optional drug compound features (rFCFP). The model generates new gene expression vectors representing distinct cellular states by manipulating z_sem (interpolation, addition) and decoding via DDIM reverse sampling.

## Method

Three components: (1) semantic encoder Enc producing z_sem from input expression x_0 plus optional auxiliary features (e.g., drug rFCFP fingerprint), (2) forward diffusion q(x_t | x_{t-1}) = N(sqrt(1-beta_t) x_{t-1}, beta_t I) with T=1000, beta_t linearly spaced 0.001 to 0.01, (3) DDIM-based reverse process with noise predictor epsilon_theta(x_t, t, z_sem).

Single-step forward: x_t = sqrt(alpha_bar_t) x_0 + sqrt(1 - alpha_bar_t) epsilon, epsilon ~ N(0, I).

Denoiser predicts clean signal: f_theta(x_t, t, z_sem) = (1/sqrt(alpha_bar_t)) (x_t - sqrt(1 - alpha_bar_t) epsilon_theta(x_t, t, z_sem)).

Conditional reverse step is deterministic DDIM: p_theta(x_{t-1} | x_t, z_sem) = N(f_theta(x_t, t, z_sem), 0) at t=1, otherwise q(x_{t-1} | x_t, x_0 = f_theta(x_t, t, z_sem)).

Training objective: L = sum_{t=1..T} E_{x_0, eps_t} || epsilon_theta(x_t, t, z_sem) - eps_t ||_2^2, eps_t ~ N(0, I). Optimizer Adam, learning rate 1e-4.

Architecture: semantic encoder is an MLP with linear layers, batch norm, ReLU. Denoiser is initial linear layer plus a sequence of MLP blocks with timestep embeddings (sinusoidal position encoding psi(t)) and z_sem injection; SiLU activations.

No classifier-free guidance is reported. The diffusion runs in the gene expression space (not in a separate compressed VAE latent). The "latent" referred to is the semantic conditioning vector z_sem, not a compressed pixel-style latent. So strictly this is a conditional DDIM with semantic encoder, not a Stable-Diffusion-style two-stage latent diffusion. For drug conditioning, a PRNet-style adapter ingests rFCFP fingerprints: z_sem' = Enc(x_0, rFCFP).

Latent manipulation modes: (a) interpolation between z_sem of two endpoints to predict intermediate states (e.g., differentiation timepoints), (b) addition of a perturbation vector (delta_z_sem) learned from one cell type and transferred to another for cross-cell-type generalization.

Hardware: NVIDIA H100 80GB HBM3. Training a dataset with 8000 cells and 500 genes takes ~15 minutes.

## Datasets / experiments

- Synthetic Splatter-simulated data: 3 cell types, 7,200 cells (2,400 train/test per type).
- iPSC differentiation (Zenodo 3625024): iPSC to mesendoderm to definitive endoderm; train on day 0 and day 3, predict day 1 and day 2.
- Gene perturbation (Norman et al., GSE133344): K562 cells; train on control, +PTPN12, +ZBTB25 singles; test on PTPN12+ZBTB25 double.
- 4i melanoma drug combination data (ETH Zurich collection 609681): ~50 features per condition.
- Glioblastoma drug response (GSE148842): myeloid + tumor + oligodendrocyte cells; train myeloid on 6 drugs and other types on etoposide; predict all drug-by-cell-type combinations.
- Sci-Plex3 (GSM4150378): random split and unseen-drug split.
- Blood vessel organoids (BVOs, figshare 27948633): endothelial, mural, fibroblasts, progenitors; full sequencing at day 11; predict days -1 to 17 plus 4Gy neutron irradiation and G-CSF stimulation.

## Metrics

Reported metrics: Pearson correlation R, R^2, Jaccard index over top differential expression gene sets, Spearman correlation of predicted pseudotime vs real timepoints, cell-type proportion comparisons. No MMD, no Wasserstein, no energy distance reported as primary metrics. No FID-equivalent. No coverage / calibration metrics.

Concrete numbers:

- iPSC day 1: Squidiff R = 0.96 +/- 0.02 vs scGen R = 0.91 +/- 0.03 (p = 1.97e-35).
- iPSC day 2: Squidiff R = 0.95 +/- 0.01 vs scGen R = 0.89 +/- 0.04 (p = 1.48e-9).
- iPSC day 1 R^2: 0.92 +/- 0.02 vs scGen 0.84 +/- 0.05 (p = 3.97e-32).
- iPSC day 3 R^2: 0.90 +/- 0.03 vs scGen 0.79 +/- 0.06 (p = 7.22e-35).
- ZBTB25+PTPN12 double: Squidiff R = 0.84 +/- 0.05 vs GEARS R = 0.71 +/- 0.08; visual R^2 ~0.75 (Squidiff) vs ~0.40 (GEARS).
- Sci-Plex3 random split: Squidiff Pearson p = 3.58e-2; R^2 p = 1.67e-2.
- Sci-Plex3 unseen drug: Pearson p = 9.96e-4; R^2 p = 7.46e-2 (not significant).
- BVO endothelial cells in mural progenitor clusters: Squidiff 1.8 +/- 0.3 vs scGen 0.5 +/- 0.2 (p = 1.02e-3).
- BVO irradiated fibroblasts: R ~ 0.85-0.90, R^2 ~ 0.72-0.81 (read off scatter plots).

Baselines: scGen (iPSC, BVO), GEARS (gene perturbation), PRNet (unseen drug), Monocle / PAGA / Scorpius (trajectory comparison). No ablations on architecture choices reported.

## What we steal for ConfPert

Squidiff is wrappable predictor in the diffusion family. ConfPert wraps Squidiff sample sets with conformal coverage layers. Concretely, for each (control x, perturbation p) query we draw N DDIM samples from Squidiff conditioned on z_sem(p), then build a conformal prediction set (per-gene CQR or split conformal on a held-out perturbation calibration fold), giving marginal or perturbation-conditional coverage guarantees over the predicted expression vector.

Squidiff is one of the 9-10 wrappable diffusion / generative perturbation predictors in our wrappable pool; it sits alongside CellOT, sVAE+, scGen, CPA, GEARS as a reference predictor. We do not compete with Squidiff on diffusion modeling capability. We only wrap it.

Useful properties for wrapping: (i) DDIM is deterministic given z_sem and noise seed, so calibration sets are reproducible; (ii) z_sem is exposed and editable, enabling perturbation-conditional calibration sets; (iii) per-cell sampling is cheap relative to training.

## What we wrap

GitHub: https://github.com/siyuh/Squidiff (main package), https://github.com/siyuh/Squidiff_reproducibility (reproducibility scripts). Zenodo archive of code: https://doi.org/10.5281/zenodo.15061773. Authors release training scripts and config; checkpoint availability for the BVO and iPSC models is via the reproducibility repo and figshare data deposit. Need to confirm at integration time which exact .pt checkpoints are public for the iPSC, GBM drug response, Sci-Plex3, and BVO experiments. If a given experiment does not ship a checkpoint, retraining on the listed public datasets is feasible (~15 min on H100 for 8000 cells x 500 genes).

## Failure modes / caveats

1. Sampling cost. Diffusion sampling is iterative; ConfPert needs N samples per query for conformal calibration, multiplying inference cost by N x DDIM_steps per cell. Plan to use reduced DDIM step counts (e.g., 50) and possibly distillation if tight latency targets emerge.
2. Linearity assumption on z_sem. Authors state interpolation in z_sem is approximate and produces averaged trajectories; "the current assumption of linearity in semantic variables may only provide approximate predictions in highly complex scenarios." Implication for ConfPert: conformal coverage may be loose for trajectory midpoints because the underlying predictive distribution has bias, and conformal corrects for miscoverage but not for systematic bias relative to true biology.
3. Unseen-drug OOD weakness. "When Squidiff encounters a completely unseen drug, prediction is limited." Sci-Plex3 unseen-drug R^2 was not statistically significant. ConfPert-wrapped Squidiff will need explicit OOD detection or shift-aware conformal (cf. Tibshirani 2019, WCP shift) when querying drugs outside training set. This is exactly the regime where ConfPert's coverage-under-shift machinery earns its keep.
4. No formal ablations on architecture / conditioning depth.
5. No distributional metric (MMD, Wasserstein) reported in the paper, so reproducing prior single-cell perturbation benchmarking conventions requires recomputing those ourselves on the released checkpoints.
6. Latent space averaging across timepoints with different growth factor schedules (the day 2 / day 3 z_sem mismatch) suggests heteroskedastic prediction error structure that CQR-style conformal handles better than fixed-width split conformal.
7. Diffusion training cost. "diffusion models generally require more computational resources compared to other generative frameworks, such as VAEs or generative adversarial networks."
8. Limited in vivo validation. Authors acknowledge "further validation using in vivo models would strengthen its translational potential."
9. 4i dataset has limited cellular features and impacts capacity to model complex drug interactions.

## Code URLs

- https://github.com/siyuh/Squidiff
- https://github.com/siyuh/Squidiff_reproducibility
- https://doi.org/10.5281/zenodo.15061773
- bioRxiv preprint: https://www.biorxiv.org/content/10.1101/2024.11.16.623974v1
- Nature Methods: https://www.nature.com/articles/s41592-025-02877-y
- PMC: https://pmc.ncbi.nlm.nih.gov/articles/PMC12407682/

## Verbatim quotes worth keeping

1. "Squidiff is a conditional denoising diffusion implicit model generating new transcriptomes representing distinct cellular states."
2. "To generate new single-cell gene expression data over time and in response to stimuli, Squidiff employs two methods of latent manipulation: interpolation and addition."
3. "biologically meaningful information, such as cell type, environmental changes, and disease states, are encoded via the semantic latent variable zsem."
4. "When Squidiff encounters a completely unseen drug, prediction is limited. Unlike cell type transfer learning, where Squidiff can generalize across different cell types, predicting the response to an entirely new drug is challenging, as the model has never been exposed to its molecular characteristics during training."
5. "the current assumption of linearity in semantic variables may only provide approximate predictions in highly complex scenarios, requiring future refinement."
