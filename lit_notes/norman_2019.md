# Norman et al. 2019 - K562 Dual-CRISPRa Perturb-seq

**Venue:** Science 2019, 365(6455):786-793, doi:10.1126/science.aax4438. GSE133344.

## Claim
Norman 2019 reports a dual-gRNA CRISPRa Perturb-seq experiment in K562 cells covering 287 perturbations across roughly 110,000 single cells, with a linear genetic interaction (GI) decomposition delta_AB = c1 delta_A + c2 delta_B + epsilon for the transcriptional effect of double overexpression. The linear model captures most variance on average (mean R^2 = 0.71). The single-cell resolution exposes within-condition phenotypic heterogeneity that is not predicted by mean-on-mean models, including a multimodal distribution of GI residuals across pairs and within-pair variation along the GI manifold for specific examples like DUSP9/MAPK1.

## Method
The CRISPRa library targets 112 hit genes whose activation enhances or retards K562 growth. The growth screen tests 28,680 unique sgRNA pairs. The single-cell readout focuses on 287 conditions (132 gene pairs plus singles and controls) with median 273 cells per condition. K562 cells stably express the SunTag CRISPRa system and are transduced with the dual-sgRNA CRISPRa GI library, then profiled by Perturb-seq. The transcriptional GI score for each pair is computed from a linear regression of the doubly-perturbed mean profile on the two singly-perturbed mean profiles. The GI manifold is visualized by projecting mean profiles for the 287 perturbations through UMAP into 2D, with stable clusters defined via HDBSCAN. For DUSP9/MAPK1 the authors compute a principal curve through the joint single-cell distribution of unperturbed, singly-perturbed, and doubly-perturbed cells to order cells by phenotype.

## Datasets / experiments
GEO accession GSE133344. Library: 112 hit genes, 28,680 sgRNA pairs in growth screen, 287 conditions in single-cell readout, 132 gene pairs in the Perturb-seq subset. Roughly 110,000 single cells; median 273 cells per condition. The authors note that as few as 50 cells per perturbation could be sufficient for some downstream uses, implying ~10^6 cells would be enough to extend to all protein-coding genes. Phenotype categories described include buffering, suppression, synthetic sick/lethal-like (SSL), and epistasis modes. Concrete biological examples include CBL/CNN1 driving erythroid differentiation (hemoglobin induction 6-39-fold, SLC25A37 13-fold) and DUSP9/MAPK1 producing a graded distribution of doubly-perturbed states.

## Metrics
Mean-on-mean OLS R^2 for the linear GI model (mean R^2 = 0.71 across pairs, with a lower mode of pairs deviating from linear additivity). UMAP 2D embedding of mean profiles. HDBSCAN cluster stability over UMAP. Principal-curve fit through the joint single-cell distribution for selected GIs, with median-filtered gene expression along the curve used to read out classes of transcripts regulated by either perturbation.

## What we steal for ConfPert
Norman is dataset #1 for ConfPert. The within-perturbation heterogeneity quote, "Stochastic heterogeneity can cause individual cells (dots) bearing a given genetic perturbation to explore the space on the GI manifold surrounding the average direction of travel (arrows)," is direct support for ConfPert's distributional framing: the paper itself acknowledges that the mean-on-mean GI score collapses real per-cell variation that is not random noise but structured exploration of a manifold. The 287-perturbation, 273-median-cells/condition structure gives us per-perturbation distributional populations to calibrate against, and the 132 gene pairs in the Perturb-seq subset give a pre-curated set of held-out doubles for combinatorial evaluation. The DUSP9/MAPK1 principal-curve analysis is a precedent for treating per-cell phenotype as a continuous distribution rather than a point estimate, and the lower mode of GI residuals where the linear model fails is exactly the regime where ConfPert's calibrated prediction sets should be widest.

## What we wrap
Dataset, not model. Norman GSE133344 is one of 5 datasets in the ConfPert evaluation suite. We use it for K1 (single-perturbation calibration on the 112 single-CRISPRa conditions) and for combinatorial evaluation (held-out double from the 132 gene pairs). The linear GI model serves as a baseline conditional mean predictor whose residual structure ConfPert must calibrate over.

## Failure modes / caveats
Norman's own analysis is mean-on-mean. The single-cell heterogeneity is described qualitatively for a few examples (DUSP9/MAPK1, CBL/CNN1) rather than quantified across all 132 pairs, so the "ground truth" for distributional calibration must be reconstructed by us from the raw counts. Single cell line (K562, chronic myeloid leukemia, erythroid-biased) limits population-level inference. CRISPRa is overexpression, not loss-of-function, so coefficients c1, c2 measure additive overexpression effects only; this constrains transfer to KO/KD datasets like Replogle and Adamson. Median 273 cells/condition is fine for population means but limits within-condition tail estimation, which matters for the conformal coverage we need at high alpha. The lower mode of pairs where R^2 is poor is small in the Norman analysis but is exactly the subset where distributional methods should outperform; our eval needs to stratify on this.

## Code URLs
GEO GSE133344 (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE133344). scperturb.org mirror (Peidli et al. h5ad standardization). Original Weissman lab repo for analysis pipeline is referenced in the paper supplementary materials.

## Verbatim quotes worth keeping
1. "Stochastic heterogeneity can cause individual cells (dots) bearing a given genetic perturbation to explore the space on the GI manifold surrounding the average direction of travel (arrows)."
2. "This linear model of transcriptional GIs explained more than 70% of the variance in gene expression on average (mean R^2 = 0.71)."
3. "Cells overexpressing both DUSP9 and MAPK1 showed a range of phenotypes spanning the transcriptional states observed in cells overexpressing either DUSP9 or MAPK1 alone."
4. "This variation did not appear to be the result of stable differences in the expression of MAPK1 and DUSP9, suggesting a possible role either for historical differences or stochastic gene expression."
5. "A relatively small number of GIs (lower mode in Fig. 4B) deviated from the expectation given by the linear model."
6. "The single-cell resolution afforded by Perturb-seq can reveal phenotypic heterogeneity for some GIs that we reasoned could yield further insight into mechanism."
7. "As our single-cell approach is sensitive to multiple outcomes or perturbations with incomplete penetrance, it is a natural strategy to pursue combinatorial searches for factors driving (trans)differentiation."
8. "To order cells in an unbiased way by 'phenotype,' we computed a principal curve measuring the path of maximum variation in the data set."

ConfPert context: Norman is the canonical CRISPRa dual-perturbation Perturb-seq dataset. We use it for K1 (single-pert calibration) AND combinatorial (held-out double) evaluation. The paper's own admission of stochastic heterogeneity around the mean GI direction is the strongest in-paper justification for moving from point-prediction GI scores to calibrated distributional prediction.
