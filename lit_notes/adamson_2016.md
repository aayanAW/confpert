# Adamson et al. 2016 - K562 UPR Perturb-seq

**Venue:** Cell 2016. GSE90546.

## Claim
Adamson 2016 introduces Perturb-seq, a multiplexed single-cell CRISPR screening platform that pairs droplet-based scRNA-seq with a genetically encoded guide barcode (GBC) so that pooled CRISPRi perturbations can be read out at single-cell transcriptome resolution. The paper applies the platform to dissect the mammalian unfolded protein response (UPR) in K562 cells, decoupling the three canonical UPR branches (ATF6, IRE1 alpha, PERK) and surfacing branch-specific and bifurcated activation patterns under perturbation.

## Method
CRISPRi (dCas9-KRAB) in K562 with a custom Perturb-seq lentiviral vector that carries two cassettes: an RNA polymerase III sgRNA expression cassette and an RNA polymerase II driven GBC expression cassette terminated by a strong BGH polyadenylation signal. The 3 prime end of the polyA transcript carries the guide barcode, so droplet 3 prime scRNA-seq (10x Chromium) recovers the GBC alongside cellular mRNA. Each droplet then contributes a triple of indices: cell barcode (CBC) for cell identity, unique molecular identifier (UMI) for transcript counting, and GBC for perturbation identity. Verbatim: "we built a platform to genetically encode a third type of index on a synthetic polyadenylated transcript. This index, which we term a 'guide barcode' (GBC), can mark specific cell perturbations (e.g., the identity of a Cas9-targeting single guide RNA, sgRNA) and thus allows complex pools of cells to be interrogated in parallel on existing droplet-based platforms." Two genome-scale CRISPRi screens (CRISPRi-v1: 15,977 genes / 20,899 TSSs at 10 sgRNAs per TSS; CRISPRi-v2: 18,905 genes / 20,526 TSSs at 5 sgRNAs per TSS) were used to nominate hits whose repression activates an ER-stress reporter; ~100 hits were then re-screened by Perturb-seq.

## Datasets / experiments
Three Perturb-seq experiments are reported. (1) Focused UPR Perturb-seq: ~65,000 transcriptomes, 91 sgRNAs targeting 82 genes plus 2 negative controls, with even sgRNA representation (457 +/- 108 cells per sgRNA). (2) Three-branch UPR epistasis experiment: ~15,000 cells perturbing ATF6, ERN1 (IRE1 alpha), and EIF2AK3 (PERK) singly and in all combinations, with and without thapsigargin treatment, to dissect crosstalk among the three sensors. (3) Genome-scale CRISPRi screens (bulk pooled, FACS readout) feeding hit selection. The scperturb / common-subset mirror of GSE90546 that downstream benchmarks load is much smaller than the full ~65k-cell focused experiment, and is the subset typically reported as ~6k cells across 8 perturbations after QC; consult the actual h5ad shipped by scperturb when reporting numbers in ConfPert tables. UPR-relevant genes profiled include HSPA5 (BiP), ATF6, ERN1, EIF2AK3, DDIT3, ATF4, XBP1, the SEC61 translocon subunits, and aminoacyl tRNA synthetases.

## Metrics
Mean transcriptome per perturbation across cells, hierarchical clustering of perturbation centroids in expression space, branch activation scores constructed from canonical UPR target genes for ATF6 / IRE1 alpha / PERK arms, low-rank ICA (LRICA) for unsupervised decoupling of co-varying transcriptional programs, and differential expression between perturbed cells and non-targeting controls. Recovery of known UPR biology (HSPA5 knockdown activating all three branches; SEC61 knockdown selectively activating IRE1 alpha) is used as the validation metric, alongside epistasis structure across single and double sensor knockdowns.

## What we steal for ConfPert
Adamson 2016 is dataset 3 for ConfPert. Smaller scale than Norman 2019 or Replogle 2022 but canonical for the UPR axis: HSPA5, ATF6, ERN1, EIF2AK3, DDIT3, ATF4, XBP1, SEC61 components. Useful as an across-dataset coverage replication test: a calibration / coverage method tuned on Norman or Replogle should still hit nominal coverage on Adamson UPR perturbations despite the dataset shift in cell count, perturbation set, and stress-response biology. Particularly attractive because the UPR has well-characterized branch-level ground truth, so coverage can be conditioned on biologically meaningful axes (per-branch activation score, per-perturbation centroid).

## What we wrap
Dataset, not model. Use the GSE90546 deposit, mirrored at scperturb.org and on Zenodo as a clean h5ad with cell x gene matrix and per-cell perturbation labels derived from GBC calls. ConfPert wraps the AnnData object only; we do not re-implement the GBC calling pipeline or the LRICA decomposition.

## Failure modes / caveats
Small effective sample size after QC (the commonly redistributed scperturb subset is ~6k cells across ~8 perturbations) limits per-perturbation statistical power for tail-of-distribution coverage estimation; a calibration set carved out of this dataset will be small. CRISPRi knockdown is incomplete and variable across cells, so per-cell perturbation strength is heterogeneous in a way the GBC label does not capture. Authors flag intermolecular provirus recombination during lentiviral pooling as a barcode-scrambling risk: "intermolecular provirus recombination during transduction can scramble barcode identities in pooled lentivirus preparations." Authors also flag analytical noise as the dominant scaling barrier: "By far the biggest barrier we anticipate is on the analytical side. Perturb-seq generates massive amounts of intrinsically noisy data." K562 is a CML line, so any UPR conclusions are conditional on that genetic background and may not transfer to non-cancer epithelial cells. Bulk-level analyses can mask the bifurcated single-cell structure the paper highlights, so any bulk-pseudobulked baseline ConfPert builds on this dataset will throw away the headline phenomenon.

## Code URLs
GEO GSE90546 (raw data and processed matrices). Mirror copies and curated h5ad versions at scperturb.org and on Zenodo. PMC: PMC5315571. DOI: 10.1016/j.cell.2016.11.048.

## Verbatim quotes worth keeping
"Functional genomics efforts face tradeoffs between number of perturbations examined and complexity of phenotypes measured. We bridge this gap with Perturb-seq, which combines droplet-based single-cell RNA-seq with a strategy for barcoding CRISPR-mediated perturbations, allowing many perturbations to be profiled in pooled format."

"Single-cell analyses decoupled the three UPR branches, revealed bifurcated UPR branch activation among cells subject to the same perturbation, and uncovered differential activation of the branches across hits, including an isolated feedback loop between the translocon and IRE1 alpha."

"Close proximity of the GBC and the BGH pA within this cassette favors faithful transmission of GBC sequences into single-cell RNA-seq libraries, which typically capture only the 3 prime ends of transcripts."

"Repression of HSPA5, which encodes the major ER chaperone BiP, robustly activated all three branches. Repression of aminoacyl tRNA synthetases activated both IRE1 alpha and ATF4 transcriptional programs. Finally, repression of all three subunits of the translocon (SEC61A1/SEC61G/SEC61B) appeared to selectively activate only the IRE1 alpha branch."

"These results suggest a model in which IRE1 alpha actively monitors the number of translocons (and perhaps function) and increases them as needed."

"Bulk RNA-seq would in this case describe a state that no cell actually occupies."

ConfPert context: Adamson 2016 is a small but canonical Perturb-seq dataset. Useful for sanity-check across-dataset coverage replication. Pair with Norman 2019 and Replogle 2022 to cover three K562 CRISPRi Perturb-seq datasets at three different scales and three different perturbation focuses (combinatorial gain-of-function, genome-scale loss-of-function, and UPR-focused loss-of-function).
