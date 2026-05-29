# Replogle et al. 2022 - Genome-Scale Perturb-seq

**Venue:** Cell 2022. DOI 10.1016/j.cell.2022.05.013.

## Claim
Replogle 2022 is the first genome-scale Perturb-seq study, profiling >2.5M single cells across three CRISPRi screens. The screens cover (a) K562 genome-wide targeting all expressed genes (n=9,866) sampled day 8 post-transduction, (b) K562 common-essential targeting 2,057 genes sampled day 6, and (c) RPE1 common-essential sampled day 7. The paper establishes that pooled CRISPRi + 10x scRNA-seq + nonparametric distributional testing recovers a high-resolution genotype-phenotype map, including rare regulators, mitochondrial-vs-nuclear partitioning, and cell-line-shared chromosomal-instability programs.

## Method
Library design: dual-sgRNA-per-element CRISPRi library, two distinct sgRNAs per gene per element, delivered lentivirally to dCas9-KRAB cells. 10x Genomics 3' scRNA-seq with direct sgRNA capture. Curated set of non-targeting control (NTC) sgRNAs in every batch enables internal z-normalization of expression to correct batch effects. For per-gene differential expression, they apply the Anderson-Darling (AD) two-sample test on z-scored single-gene distributions of perturbed vs control cells. For global transcriptional change, they apply a permuted energy-distance test in PC space (cells x principal components) comparing perturbed vs control populations. Leverage scores quantify single-cell deviation from the control manifold.

## Datasets / experiments
- K562 genome-wide (day 8): 9,866 genes, ~11,258 sgRNA-pair elements, mean ~183 / median ~171 cells per perturbation, >2.5M total high-quality cells.
- K562 essential (day 6): 2,057 genes, ~2,285 sgRNA-pair elements, mean ~148 / median ~124 cells per perturbation.
- RPE1 essential (day 7): ~2,679 sgRNA-pair elements, mean ~101 / median ~79 cells per perturbation.
- Cross-cell-line comparison: chromosomal-instability (CIN) score Pearson r = 0.69 between K562 and RPE1; cophenetic correlation between hierarchical clusterings = 0.37; median pairwise correlation across shared perturbations r = 0.23.
- Mitochondrial-vs-nuclear classification: random-forest accuracy 0.64 on mt-genome perturbations vs 0.25 on nuclear (gene-level partitioning).

## Metrics
Anderson-Darling p-values (per-gene DE), permuted energy-distance test p-values (global), single-cell leverage scores, perturbation-perturbation Pearson correlation matrices, cophenetic correlation across cell lines, random-forest classification accuracy.

## What we steal for ConfPert
CRITICAL load-bearing self-citation. Replogle 2022 is the canonical genome-scale Perturb-seq paper that EXPLICITLY chose nonparametric distributional tests over mean-based statistics. Their stated rationale "we chose to use conservative, nonparametric statistical tests to detect transcriptional changes rather than making specific assumptions about the underlying distribution of gene expression levels" and the AD test selection because it is "sensitive to transcriptional changes in a subset of cells, enabling us to find differences even with incomplete penetrance" together form the data-side prior art that ConfPert's calibrated coverage on prediction-side distributions wraps. The AD/energy-distance choice is the analog on the empirical side of what split-conformal calibration is on the predictive side: distribution-free, finite-sample-valid, and robust to the heterogeneity in perturbation response that mean-shift evaluation hides.

## What we wrap
Dataset, not model. Replogle K562 + RPE1 are 2 of 5 datasets in the ConfPert evaluation suite. The cross-cell-line K562 -> RPE1 split (cophenetic 0.37, median pairwise r 0.23) is the canonical covariate-shift test for ConfPert weighted conformal prediction: target conditional shifts on the marginal cell-state distribution while preserving shared CIN-type biology (r=0.69), giving a real distribution-shift benchmark with known structure rather than a synthetic perturbation.

## Failure modes / caveats
- Single sgRNA-pair element per gene in each screen: no within-screen replicate allowing direct biological-replicate variance estimation at the element level.
- Limited cell counts per perturbation, especially RPE1 essential (median 79). Coverage at the lower tail constrains how tight any conformal interval can be without aggressive weight stabilization.
- Off-target / neighbor knockdown: ~7.5% of perturbations cause knockdown of a neighboring gene. ConfPert must either filter or explicitly model these as confounded.
- Day-of-harvest differs across screens (6/7/8 days). Time-axis confound between K562-essential and K562-genome-wide is not the same biological state.
- Gene-program heterogeneity within a single perturbation is exactly what motivates AD/energy-distance over t-tests, but it also means any pointwise predictor evaluated only on means will systematically underestimate true uncertainty.

## Code URLs
- Interactive + processed data: https://gwps.wi.mit.edu/
- Raw SRA: BioProject PRJNA831566
- scperturb.org mirror with harmonized AnnData

## Verbatim quotes worth keeping
1. "In general, we chose to use conservative, nonparametric statistical tests to detect transcriptional changes rather than making specific assumptions about the underlying distribution of gene expression levels."
2. "To detect individual differentially expressed genes (DEGs), we applied the Anderson-Darling (AD) test that is sensitive to transcriptional changes in a subset of cells, enabling us to find differences even with incomplete penetrance."
3. "First, we examined global transcriptional changes using a permuted energy distance test."
4. "We compared cells bearing each genetic perturbation with control cells at the level of principal components."
5. "These allowed for internal z-normalization of expression measurements to correct for batch effects."
6. "Leverage scores quantify how outlying each perturbed cell's transcriptome is relative to control cells."
7. "We retained >2.5 million high-quality cells with a median coverage of >100 cells per perturbation."

ConfPert framing note: quotes (1) and (2) together are THE canonical citation justifying ConfPert's distributional-evaluation framing on the data side. Quote (1) establishes the field-standard rejection of parametric mean-based testing for Perturb-seq; quote (2) names "incomplete penetrance" as the explicit biological reason. ConfPert extends this from observed-data testing to predictive-distribution coverage: if Replogle 2022 already concedes that perturbation responses are heterogeneous enough to require nonparametric tests on observed cells, then any predictive model evaluated only by mean-shift correlation is being held to a strictly weaker bar than the data-side gold standard, and calibrated distributional prediction with finite-sample coverage guarantees is the natural prediction-side analog.
