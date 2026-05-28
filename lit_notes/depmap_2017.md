# Tsherniak et al. 2017 - DepMap

**Venue:** Cell 2017. DOI 10.1016/j.cell.2017.06.010.

## Claim

DepMap is a systematic, genome-scale loss-of-function effort that defines gene essentiality and differential cancer dependencies across hundreds of human cancer cell lines. The 2017 paper analyzed 501 RNAi screens to identify 769 strong differential dependencies and to predict, from cell-line molecular features, which genes a given cancer model relies on. The work proposes the cancer dependency map as a foundation for prioritizing therapeutic targets at scale, and was extended after 2017 with CRISPR-Cas9 screens in Project Achilles to grow the resource.

## Method

Genome-scale RNAi using pooled shRNA libraries delivered by lentivirus. Each cell line was screened in quadruplicate over at least 16 population doublings or 40 days, with shRNA abundance measured by sequencing to score depletion or enrichment. To deconvolve the well-known shRNA off-target problem, the authors developed DEMETER, a computational framework that decomposes observed shRNA dropout into gene-level on-target effects and seed-sequence-driven off-target effects across the entire panel of cell lines simultaneously. DEMETER replaces the earlier ATARiS scoring used in Project Achilles. Cell lines were jointly characterized by CCLE molecular data (expression, mutation, copy number, methylation), and elastic-net / random-forest predictive models were trained to map molecular features to dependency scores. A 6-sigma threshold relative to the panel mean was used to call strong differential dependencies.

## Datasets / experiments

- 501 cancer cell lines screened with genome-scale shRNA library (~17,098 genes targeted, multiple shRNAs per gene).
- 769 strong differential dependencies identified at the 6-sigma threshold.
- Predictive models built for 426 dependencies, with 55% successfully predicted from molecular features.
- 66,646 candidate molecular features evaluated as potential biomarkers.
- 74% of cell lines had at least one 6-sigma dependency that the authors describe as a readily druggable target.
- After 2017 the resource scaled into Project Achilles CRISPR-Cas9 screens (genome-wide Cas9 + sgRNA libraries, scored with CERES and later Chronos), now over 2000 cancer models on the DepMap portal, integrated with CCLE multi-omics, PRISM repurposing drug viability, Project Score (Sanger), and additional proteomics layers (Olink, Sanger).

## Metrics

- DEMETER score: gene-level dependency, on-target component recovered from shRNA dropout (replaces ATARiS).
- DEMETER2: refined version released later in the DepMap pipeline.
- CERES: gene effect score for CRISPR screens that corrects for copy-number-driven cutting toxicity.
- Chronos: current CRISPR scoring algorithm on DepMap, models read counts dynamically across screen time and reduces copy-number bias further than CERES.
- 6-sigma cutoff for "strong differential dependency" relative to the cross-line mean.
- Predictive model accuracy reported per dependency; 82% of best-predicting biomarkers are RNA-expression features, 16% DNA mutations.

## What we steal for ConfPert

DepMap is the cell-line-context substrate for K3 in ConfPert. It supplies the molecular characterization (CCLE expression, mutation, copy number, methylation) and the dependency phenotype that anchor "what this cell line is" beyond the perturbation label. Combined with PRISM drug-viability screens on the same cell-line panel, DepMap gives a drug-cell-line response surface against which ConfPert's calibrated cell-population predictions can be tested for whether predicted resistant subpopulations align with cell lines (or pseudo-subpopulations within a line) whose dependency profile is consistent with the observed insensitivity. Concretely, ConfPert can use DepMap-derived cell-line embeddings as conditioning covariates and use PRISM viability as an external validation signal that is independent of the transcriptomic perturbation atlas used for training.

## What we wrap

Dataset, not a model. The DepMap portal at depmap.org distributes preprocessed gene effect matrices, CCLE omics, and PRISM drug response. ConfPert wraps DepMap as a conditioning-feature provider and as an evaluation oracle, not as a learned component. We freeze a release tag (e.g. 26Q1 or whichever quarterly release we lock to) and document it.

## Failure modes / caveats

- Bulk-cell-line gene essentiality is not the same object as transcriptomic single-cell perturbation response. Linking K3 (DepMap context) to ConfPert's calibrated single-cell distribution prediction is an indirect bridge and must be justified explicitly: we are using DepMap as covariate context and as an orthogonal viability signal, not as a label that ConfPert directly predicts.
- Off-target shRNA effects: DEMETER mitigates but does not eliminate seed-mediated artifacts. CRISPR screens have their own copy-number-amplification artifacts that CERES and Chronos partially correct.
- Cell-autonomous viability only: DepMap measures cell-intrinsic survival in 2D culture. It does not capture immune, stromal, or microenvironmental dependencies, so claims about in-vivo resistance subpopulations from DepMap context are weak.
- Cell-line representativeness: the panel is biased toward established, easy-to-grow cancer lines and underrepresents some lineages.
- The 2017 paper itself flags that "the comprehensive identification and prediction of dependencies will require a substantial increase in the number and diversity of cell lines analyzed."
- Most cancer dependencies are not somatically mutated or focally amplified, so simple genomics-driven target prioritization misses most of the signal.

## Code URLs

- https://depmap.org/portal/
- https://depmap.org/rnai (2017 RNAi data release)
- https://www.broadinstitute.org/achilles (Project Achilles)

## Verbatim quotes worth keeping

1. "We analyzed 501 genome-scale loss-of-function screens performed in diverse human cancer cell lines."
2. "We developed DEMETER, an analytical framework that segregates on- from off-target effects of RNAi."
3. "We identified a set of 769 strong differential gene dependencies for which the DEMETER scores of at least one cell line were six standard deviations (6 sigma) or greater from the mean."
4. "The vast majority of predictable differential dependencies (82%) were best predicted by RNA expression levels, whereas DNA mutation accounted for only 16%."
5. "Together, these observations provide a foundation for a cancer dependency map that facilitates the prioritization of therapeutic targets."
6. "We propose that a concerted, international effort should be launched to create a definitive cancer dependency map."
7. "Our observations indicate that the comprehensive identification and prediction of dependencies will require a substantial increase in the number and diversity of cell lines analyzed."

ConfPert context: K3 background data resource. Pair with PRISM for drug-resistance subpopulation evaluation.
