# Bunne et al. 2023 - CellOT: Neural Optimal Transport for Perturbation Response

**Venue:** Nature Methods 2023, vol. 20, no. 11, pp. 1759-1768. DOI: 10.1038/s41592-023-01969-x. PMC10630137.

## Claim

CellOT learns a deterministic neural optimal transport (OT) map T_k per perturbation k, mapping the unperturbed (control) single-cell distribution to the perturbed distribution. Given any control cell x, the model predicts its perturbed counterpart as T_k(x) = grad g*_{theta_k}(x), the gradient of an input-convex neural net (ICNN) Brenier potential. Because single-cell assays destroy cells, paired control/perturbed observations of the same cell do not exist; CellOT exploits OT duality to recover the per-cell map from unpaired distributional samples. The paper claims CellOT "outperforms current methods at predicting single-cell drug responses, as profiled by scRNA-seq and a multiplexed protein-imaging technology" and generalizes to held-out patients (lupus IFN-beta, glioblastoma panobinostat), across species (LPS innate immune response in pigs, rabbits, mice, rats), and to hematopoietic developmental trajectories.

## Method

CellOT instantiates the Makkuva et al. 2020 neural-OT formulation. Two ICNNs g_theta and f_phi parameterize Brenier dual potentials. Training solves a min-max dual: "(g_theta*, f_phi*) <- argmax_phi min_theta C_{rho_c, rho_k} - V_{rho_c, rho_k}(g_theta, f_phi), where T* = grad g*_theta". Each ICNN layer is "h_{i+1} = sigma_i(W_i^x x + W_i^z h_i + b_i) and f(x; theta) = h_L" with non-negativity constraints on the W^z weights to enforce convexity in x. Default architecture: 4 hidden layers, width 64; learning rate 1e-4; batch size 256. Trained per (cell-type, perturbation) pair: "CellOT learns the optimal pair of dual potentials (g_{theta_k}*, f_{phi_k}*) ... for each perturbation k in K". Prediction is deterministic and one-to-one: "CellOT then predicts the transformation of a control cell x_i^c upon perturbation k via x_hat_i^k = grad g*_{theta_k}(x_i^c)". Pushing N control cells through grad g* yields N predicted perturbed cells, so the empirical pushforward is a sample-producing distributional predictor (not a Gaussian mean predictor and not a generative model with explicit density). Operates either in raw feature space or on a learned latent (autoencoder) space depending on dataset.

## Datasets / experiments

- **4i (multiplexed iterative immunofluorescence imaging) on melanoma:** "97,748 cells ... two melanoma tumor cell lines (ratio 1:1) ... 34 different treatments". Per-marker protein readouts.
- **sci-Plex 3 (scRNA-seq):** drug screens with "34 and 9 different treatments" across cell lines.
- **Lupus IFN-beta scRNA-seq (Kang et al. 2018):** "eight patients with lupus ... response of eight patients with lupus to interferon (IFN)-beta". Used for held-out-patient (o.o.d.) generalization.
- **Glioblastoma sci-Plex (Zhao et al.):** "seven patients ... panobinostat treatment outcomes of glioblastoma patients". Held-out-patient evaluation.
- **Cross-species LPS:** innate immune stimulation in "pigs, rabbits, mice and rats ... stimulated using LPS"; CellOT reconstructs responses both i.i.d. and o.o.d.
- **Hematopoiesis:** "oligopotent and multipotent progenitor cell subpopulations" with timepoints "days 2, 4 and 6" used to model developmental flux as OT.

Comparisons are against identity (predict perturbed = control), scGen, CPA, and observed lower bounds (sample-vs-sample distances within true perturbed). Concretely: "CellOT outperforms the baselines in both metrics across all treatments, typically by one order of magnitude" on 4i. CellOT MMD curves "approach the lower bound (observed setting), whereas the baseline methods often do not improve much over the identity setting". On lupus IFN-beta o.o.d.: "All models ... show little performance drop when modeling the treatment outcome on a new patient". On cross-species: "CellOT accurately reconstructs the innate immune response in both mouse and rat in the i.i.d. and o.o.d. setting".

## Metrics

- **MMD (maximum mean discrepancy)**, multi-scale RBF kernel: "Low values of MMD imply that all moments of two distributions are matched". Primary distributional metric.
- **l2 feature means:** "the l2-distance between means of the observed and predicted distributions". First-moment metric.
- **R^2 on per-feature DE / perturbation-marker genes/proteins:** "average correlation coefficient r^2 between the predicted and observed data on all features", and on marker-gene subsets (DE genes).
- **Wasserstein-2** on subsampled distributions; sliced Wasserstein and energy distance reported in supplementary comparisons.
- **UMAP/2D projection visualizations** for qualitative inspection of pushforward overlap with true perturbed cells.

## What we steal for ConfPert

CellOT is a wrappable distributional predictor: given control cells {x_1, ..., x_N}, it returns a sample {grad g*(x_1), ..., grad g*(x_N)} that approximates the perturbed distribution. ConfPert wraps the CellOT pushforward as a base distributional predictor and adds split-conformal coverage on top of its outputs (e.g., conformal prediction sets for individual marker dimensions, MV/sliced-Wasserstein conformal balls around predicted distributions, or HPD-based conformal regions per cell). Crucially, CellOT's training objective is already a distributional (OT/Wasserstein) loss, so its uncalibrated multivariate distributional metrics (MMD, energy distance, sliced W2) should be substantially stronger than mean-only predictors like CPA's reconstruction or scGen's latent shift. ConfPert's incremental claim over CellOT is therefore not "better Wasserstein" but **finite-sample distribution-free coverage guarantees** (validity at level 1 - alpha) for held-out perturbations or held-out patients, with the K1 demo showing that conformalizing CellOT preserves CellOT's good distributional fidelity while delivering nominal coverage where vanilla CellOT has none.

## What we wrap

- Wrappable predictor #3 in the ConfPert wrapper roster.
- Repo: https://github.com/bunnech/cellot (BSD-3-Clause).
- Entry points: `scripts/train.py` (training) and `scripts/evaluate.py` (eval, supports `iid`/`ood` modes and `data_space`/`latent` spaces). Configs in `configs/tasks/` (e.g., `4i.yaml`) and `configs/models/cellot.yaml`. Example: `python ./scripts/train.py --outdir ./results/4i/drug-cisplatin/model-cellot --config ./configs/tasks/4i.yaml --config ./configs/models/cellot.yaml --config.data.target cisplatin`.
- No public pretrained checkpoints documented; ConfPert must train per (dataset, perturbation, cell-type) cell. Approx 3 hours CPU per perturbation per their README. Preprocessed data hosted externally (ETH research collection).
- Pushforward prediction interface is not exposed as a clean API; we will need a thin wrapper that loads the trained `g_theta` ICNN, autograds grad_x g_theta(x) over a batch of control cells, and returns the predicted perturbed sample matrix.

## Failure modes / caveats

- **Per-perturbation training:** one ICNN pair per perturbation k. No transfer to unseen perturbations and no shared representation across drugs. Combinatorial chemistry / unseen-compound generalization is out of scope.
- **ICNN training instability:** input-convex constraints (non-negative W^z) plus min-max dual training are notoriously brittle; convergence depends on initialization and the inner-loop schedule. Public repo offers no checkpoints, increasing reproduction variance.
- **Identifiability assumption:** "If this principle is violated, however, and perturbations strongly disrupt the population to an unidentifiable level, the performance of CellOT as well as other methods drops". Strong perturbations that destroy cell-type structure break the OT correspondence assumption.
- **No proliferation/death modeling:** "the current system is not able to recover effects (other than cell flux) that change the distribution of cells between time points, for example, proliferation and death". OT is mass-preserving; unbalanced extensions are not in this paper.
- **Deterministic per-cell mapping:** "CellOT models cell responses as deterministic trajectories" despite the "stochastic nature of cell-fate decisions". The pushforward is one-to-one; intrinsic single-cell stochasticity in the response is not represented (this matters for ConfPert: any per-cell predictive uncertainty must come from ConfPert's calibration layer, not from CellOT itself).
- **Scale:** "CellOT's generalization capacity has been evaluated on relatively small datasets" - eight lupus patients, seven glioblastoma patients. Coverage claims at this scale need careful conformal exchangeability arguments.
- **Visualization caveat:** main-text qualitative results often shown as 2D UMAP/PCA projections; full multivariate fidelity is summarized only by MMD/l2/r^2 scalars.

## Code URLs

- https://github.com/bunnech/cellot
- Paper: https://www.nature.com/articles/s41592-023-01969-x
- Open-access full text: https://pmc.ncbi.nlm.nih.gov/articles/PMC10630137/
- Author correction: https://www.nature.com/articles/s41592-023-02073-w

## Verbatim quotes worth keeping

1. "CellOT outperforms current methods at predicting single-cell drug responses, as profiled by scRNA-seq and a multiplexed protein-imaging technology."
2. "CellOT then predicts the transformation of a control cell x_i^c upon perturbation k via x_hat_i^k = grad g*_{theta_k}(x_i^c)."
3. "CellOT learns the optimal pair of dual potentials (g_{theta_k}*, f_{phi_k}*) by solving [equation] for each perturbation k in K."
4. "CellOT outperforms the baselines in both metrics across all treatments, typically by one order of magnitude."
5. "If this principle is violated, however, and perturbations strongly disrupt the population to an unidentifiable level, the performance of CellOT as well as other methods drops."
6. "the current system is not able to recover effects (other than cell flux) that change the distribution of cells between time points, for example, proliferation and death."

ConfPert context: CellOT is the strongest published per-perturbation distributional baseline because its training loss is already Wasserstein-aligned. The K1 demonstration in ConfPert is therefore: take CellOT's good-but-uncalibrated pushforward, add a conformal calibration layer on a held-out split of perturbed cells, and show that the resulting prediction sets / distributional balls hit nominal coverage on held-out patients (lupus, glioblastoma) without degrading MMD or per-feature R^2. This is a clean coverage-on-top-of-SOTA narrative.
