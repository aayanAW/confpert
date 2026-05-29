# Roohani, Huang, Leskovec 2024 - GEARS

**Venue:** Nature Biotechnology 2024 (42:927-935). DOI: 10.1038/s41587-023-01905-6. Preprint: bioRxiv 10.1101/2022.07.12.499735 (2022).

## Claim

GEARS predicts post-perturbation gene expression for single-gene and combinatorial perturbations from Perturb-seq data, including perturbation sets with genes never experimentally perturbed. A gene co-expression GNN over learnable gene embeddings is paired with a GO-derived perturbation similarity GNN over learnable perturbation embeddings; for multigene sets perturbation embeddings are summed, then conditioned on a cross-gene MLP and decoded per gene to a scalar effect added to a sampled control cell to give g-hat = z-hat + g_ctrl. Training uses an autofocus MSE with exponent (2+gamma) plus a sign-based direction loss; an optional Kendall and Gal heteroscedastic sigma^2 head exposes a per-gene log-variance trained with L_unc = exp(-s_u)(g_u - g_hat_u)^(2+gamma).

## Method

Architecture (Methods, Fig. 1b). Each gene u has a learnable embedding x_u^gene in R^d, refined by a GNN over a co-expression graph G_gene (top H_gene Pearson neighbors above threshold delta): h_u^gene = GNN_theta_g(x_u^gene, G_gene). Each candidate perturbation has a learnable embedding x_v^pert refined by a GNN over a GO-derived similarity graph G_pert built from the Jaccard index J_{u,v} = |N_u intersect N_v| / |N_u union N_v| over GO pathway sets: h_v^pert = GNN_theta_p(x_v^pert, G_pert). For a perturbation set P = (P_1,...,P_M) embeddings are summed: h_P = MLP_theta_c(sum_i h_{P_i}^pert). Per gene u this decodes to z_u; a cross-gene MLP h_cg = MLP_theta_cg(z) provides transcriptome-wide context, and a per-gene head gives z_u-hat = w_u_cg(z_u || h_cg) + b_u_cg. Final prediction g-hat = z-hat + g_ctrl, so "This allows GEARS to focus only on learning perturbation effects."

Loss (verbatim from "Autofocus direction-aware loss"): the autofocus loss is

L_autofocus = (1/T) sum_{k=1..T} (1/T_k) sum_{l=1..T_k} (1/K) sum_{u=1..K} (g_u - g-hat_u)^(2+gamma)

with elevated exponent giving "a higher weight to differentially expressed genes." The direction loss is

L_direction = (1/T) sum_k (1/T_k) sum_l (1/G) sum_u [sign(g_u - g_u^ctrl) - sign(g-hat_u - g_u^ctrl)]^2

with total prediction loss L = L_autofocus + lambda L_direction. The public utils.py has gamma = 2 hard-coded (i.e. exponent 4) and direction_lambda = 1e-3 default.

Heteroscedastic uncertainty head (Kendall and Gal 2017, ref. 52): "A Gaussian likelihood N(g-hat_u, sigma-hat_u^2) is used to model the postperturbation gene expression value for gene u under perturbation P, where g-hat_u is the predicted postperturbation scalar and sigma-hat_u^2 is the variance." A separate gene-specific MLP head predicts the log-variance s_u = log(sigma-hat_u^2) = w_u^unc h_u^post-pert + b_u^unc. The uncertainty loss is

L_unc = (1/T) sum_k (1/T_k) sum_l (1/G) sum_u exp(-s_u)(g_u - g-hat_u)^(2+gamma)

"By encouraging log variance to be large when the error is large, the log variance is learned to be a proxy of model uncertainty." Crucially: in gears.py the public predict() returns the dictionary of mean predictions when uncertainty is disabled, and a tuple (mean_dict, uncertainty_dict) when uncertainty=True. The sigma^2 head is exposed only as a confidence/regularization proxy and the default predict() path returns the mean.

## Datasets / experiments

Norman et al. 2019 (GSE133344, K562 CRISPRa, 131 two-gene perturbations, primary combinatorial benchmark); Replogle et al. 2022 K562 (1,092 perturbations); Replogle et al. 2022 RPE-1 (1,543 perturbations, >170k cells each via genome-wide Perturb-seq); Adamson et al. 2016 (GSE90546, UPR multiplexed CRISPRi); Dixit et al. 2016 Perturb-seq (GSE90063); Jost et al. 2020 titrating CRISPRi (GSE132080); Tian et al. 2019 multimodal CRISPRi (GSE124703); Replogle et al. 2020 combinatorial (GSE146194); Horlbeck et al. 2018 genetic landscape map (GSE116198). Train/test splits hold out perturbations entirely, with three two-gene generalization regimes (0/2 unseen, 1/2 unseen, 2/2 unseen) reflecting how many perturbed genes were experimentally observed during training.

## Metrics

(i) MSE on top-20 differentially expressed genes per perturbation, normalized to the no-perturbation baseline ("the vast majority of genes do not show substantial variation between unperturbed and perturbed states, we restricted our m.s.e. analysis to the harder task of only considering the top 20 most differentially expressed genes"); (ii) Pearson correlation on the change-in-expression delta over control across all genes; (iii) fraction of top-20 DE genes with predicted change in opposite direction (direction error); (iv) Jaccard similarity between predicted and true sets of top-20 DE genes; (v) precision@10 for genetic interaction subtype identification (synergy, suppression, neomorphism, redundancy, epistasis); (vi) top-10 accuracy and R^2 on GI scores. Baselines: CPA, GRN-based linear propagation adapted from CellOracle, and a no-perturbation null.

Headline results: 30-50 percent MSE improvement over CPA/GRN on top-20 DE genes for single-gene perturbations on Replogle K562 (-29.2 percent) and RPE-1 (-48.9 percent); on Norman two-gene perturbations the normalized MSE improvement reaches -53.8 percent in the 2/2-unseen regime, -47.2 percent at 1/2-unseen, and -32.4 percent at 0/2-unseen. Pearson correlation on change in expression is "more than two times better" than baselines (+382.9 percent on Replogle K562, +499.4 percent on Replogle RPE-1, +24.3 percent on Norman). Precision@10 for GI subtypes improves +39 percent (synergy), +14 percent (suppression), +56 percent (neomorphism), +252 percent (redundancy), +89 percent (epistasis), with +40 percent precision@10 averaged across four subtypes claimed in the abstract. GI score R^2 reaches "approximately 0.4 for synergy, neomorphism and redundancy, whereas it was only around 0.0 for the same interactions when predicted by CPA." When trained directly to predict cell fitness, R^2 between 0.64 and 0.93.

## What we steal for ConfPert

GEARS is wrappable predictor #4 and the most rigorous demonstration target for ConfPert. The Kendall-Gal sigma^2 head means GEARS can in principle produce a per-gene per-perturbation Gaussian, so ConfPert can wrap GEARS-with-uncertainty to give it the best possible distributional output and then layer conformal calibration on top to provide finite-sample COVERAGE that the sigma^2 head itself does not provide on its own (the head is a model-confidence regularizer, not a calibrated predictive distribution). Important framing point for ConfPert positioning: even with the sigma^2 head exposed, the marginal predictive is Gaussian, and a Gaussian has bimodality coefficient bounded above by 1/3 (well below the 5/9 threshold), so GEARS cannot recover bimodal cell-fate genes regardless of how well sigma^2 is calibrated. This becomes the canonical paragraph distinguishing ConfPert from "just calibrate GEARS' sigma^2." The autofocus exponent (2+gamma) with gamma=2 also tells us GEARS' MSE is dominated by large-deviation DE genes, so naively conformalizing GEARS residuals will over-pad easy genes; we want CQR-style or per-gene-quantile calibration.

## What we wrap

GEARS released as `cell-gears` on PyPI; latest version is 0.1.2 (released 2023-12-13). Install via `pip install cell-gears` after PyG. Use the official PertData pipeline: `pert_data = PertData('./data'); pert_data.load(data_name='norman'); pert_data.prepare_split(split='simulation', seed=1); pert_data.get_dataloader(...)`. The README references a "Tutorial on how to train an uncertainty-aware GEARS model"; we will train with uncertainty=True so we can access the (mean_dict, uncertainty_dict) tuple. Public Norman + Replogle pretrained checkpoints are not linked from the README; we either train from scratch on Modal H100 or load community checkpoints. The HetPert real-GEARS Modal H100 probe found VR=2.5e-11 because predict() returns the mean only - the fix is to enable the uncertainty branch.

## Failure modes / caveats

Knowledge-graph dependency: "it also limits the ability of GEARS to predict outcomes for perturbing previously unperturbed genes that are not well connected in this graph." Autofocus gamma is fixed (paper text writes general (2+gamma); shipped code hard-codes gamma=2 making the effective exponent 4) and is not tuned per-dataset. Direction loss weight default 1e-3. Training is CPU-prohibitive at full Replogle scale (~$5/scFM training run on Modal H100 in our setup). The sigma^2 head is a model-confidence proxy, not a calibrated predictive distribution: "the log variance is learned to be a proxy of model uncertainty." Marginal predictive is Gaussian per gene, so any bimodal/heavy-tailed cell-fate genes are out of reach by construction. The default predict() returns only the mean, which is the source of the VR=2.5e-11 collapse observed in the HetPert probe.

## Code URLs

- https://github.com/snap-stanford/GEARS, `pip install cell-gears` (PyPI 0.1.2, 2023-12-13).
- Loss: gears/utils.py `loss_fct` and `uncertainty_loss_fct`.
- Model: gears/model.py with `self.uncertainty_w = MLP(...)` for the sigma^2 head and `SGConv` layers for co-expression (`self.layers_emb_pos`) and GO (`self.sim_layers`) GNNs.

## Verbatim quotes worth keeping

1. "A Gaussian likelihood N(g-hat_u, sigma-hat_u^2) is used to model the postperturbation gene expression value for gene u under perturbation P, where g-hat_u is the predicted postperturbation scalar and sigma-hat_u^2 is the variance."

2. "We add an additional gene-specific layer to predict the log variance term s_u = log(sigma-hat_u^2) ... and learn it through a modified Bayesian neural network loss" with L_unc = (1/T) sum_k (1/T_k) sum_l (1/G) sum_u exp(-s_u)(g_u - g-hat_u)^(2+gamma).

3. "By encouraging log variance to be large when the error is large, the log variance is learned to be a proxy of model uncertainty."

4. Autofocus: L_autofocus = (1/T) sum_k (1/T_k) sum_l (1/K) sum_u (g_u - g-hat_u)^(2+gamma); "We designed an autofocus loss that automatically gives a higher weight to differentially expressed genes by elevating the exponent of the error."

5. "GEARS makes use of a Bayesian formulation to overcome this challenge by outputting an uncertainty metric that is inversely correlated with model performance."

6. From utils.py (verbatim): `losses += torch.sum((pred_p - y_p)**(2 + gamma) + reg * torch.exp(-logvar_p) * (pred_p - y_p)**(2 + gamma)) / pred_p.shape[0] / pred_p.shape[1]` with `gamma = 2` and default `reg = 0.1`.
