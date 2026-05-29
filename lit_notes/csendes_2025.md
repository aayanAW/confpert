# Csendes et al. 2025 - Foundation Cell Model Benchmark

**Venue:** BMC Genomics 26(1):393 (2025), DOI 10.1186/s12864-025-11600-2, PMID 40269681, PMCID PMC12016270. Originally posted on bioRxiv 2024-10-01 as "Benchmarking a foundational cell model for post-perturbation RNAseq prediction" (singular "foundational model"), bioRxiv DOI 10.1101/2024.09.30.615843. Authors: Gerold Csendes, Gema Sanz, Kristof Z. Szalay, Bence Szalai (Turbine Ltd., Budapest). The published BMC Genomics version expanded the bioRxiv preprint from one model (scGPT) to two (scGPT and scFoundation) and added Sanz as author.

## Claim

Csendes et al. 2025 benchmark single-cell foundation cell models (scGPT, scFoundation) against simple baselines (Train Mean, Random Forest, Elastic Net, k-NN) on post-perturbation RNA-seq prediction across four Perturb-Seq datasets. The headline finding, stated verbatim in the abstract, is: "Surprisingly, we found that even the simplest baseline model - taking the mean of training examples - outperformed scGPT and scFoundation. Furthermore, basic machine learning models that incorporate biologically meaningful features outperformed scGPT by a large margin." A second, equally load-bearing claim: "we identified that the current Perturb-Seq benchmark datasets exhibit low perturbation-specific variance, making them suboptimal for evaluating such models."

## Method

Models benchmarked:
- scGPT (transformer foundation model with perturbation tokens, fine-tuned for prediction).
- scFoundation (transformer foundation model paired with GEARS for prediction).
- Train Mean (predict the mean of training perturbation effects, a trivial baseline).
- Random Forest Regressor.
- Elastic Net Regression.
- k-Nearest Neighbors Regression.
- The basic ML models were also run with biologically meaningful input features: Gene Ontology (GO) annotations, scGPT embeddings, and scFoundation embeddings used as features for Random Forest.

Datasets (all Perturb-seq, single-cell CRISPR screens):
- Adamson 2016: 68,603 cells, 87 perturbations, CRISPRi, ER stress focused (homogeneous).
- Norman 2019: 91,205 cells, 284 perturbations, CRISPRa.
- Replogle K562: 162,751 cells, 1,093 perturbations, genome-wide CRISPRi.
- Replogle RPE1: 162,733 cells, 1,544 perturbations, genome-wide CRISPRi.

Splits: Perturbation Exclusive (PEX) split. Test perturbations are unseen during training. The authors note that current benchmarks only assess PEX, not Cell Exclusive (CEX).

Evaluation metrics:
- Pearson correlation on raw gene expression (dismissed as not meaningful).
- Pearson Delta: correlation on the differential expression vector (post-perturbation minus control). Primary metric.
- Pearson Delta DE: Pearson Delta restricted to top 20 differentially expressed genes.
- Pearson Delta excluding the targeted (knocked out) gene itself.

All evaluations are at pseudo-bulk level. There are no per-cell distributional metrics, no MMD, no Wasserstein, no energy distance.

## Datasets / experiments

Pearson Delta (full DE vector) leaderboard, approximate:

| Dataset | Train Mean | scGPT | scFoundation | RF + GO |
|---|---|---|---|---|
| Adamson | 0.711 | 0.641 | 0.552 | 0.739 |
| Norman | 0.557 | 0.554 | 0.459 | 0.586 |
| Replogle K562 | 0.373 | 0.327 | 0.269 | 0.480 |
| Replogle RPE1 | 0.628 | 0.596 | 0.471 | 0.648 |

Train Mean outperforms both foundation models on every dataset. RF + GO further outperforms Train Mean on every dataset. scFoundation is worst.

Dataset variance analysis: median pairwise correlation between perturbation effects within each dataset.
- Adamson: 0.662 (highly homogeneous, ER-stress focused).
- Replogle K562: 0.117 (most heterogeneous, genome-wide).
- Inverse relationship: higher intra-dataset correlation predicts higher Pearson Delta from all models, indicating models exploit the dataset's shared mean signal rather than perturbation-specific structure.

## Metrics

Distributional and per-cell metrics: absent. The paper operates entirely in the pseudo-bulk Pearson space. The authors explicitly justify the focus on Pearson Delta over raw Pearson because raw Pearson is dominated by baseline expression magnitudes:

> "In the raw gene expression space, all models performed similarly (Pearson > 0.95). However, the Pearson correlation values between raw gene expression profiles are strongly influenced by the baseline expression magnitudes of different genes, so we did not consider these metrics meaningful."

This is the same critique Ahlmann-Eltze et al. 2025 and Cell-Eval 2025 raise in stronger form. Csendes 2025 acknowledges the metric problem but does not propose distributional or calibration-based alternatives.

## What we steal for ConfPert

Csendes 2025 is the benchmarking-skeptic companion to Ahlmann-Eltze 2025. Both arrive at the same conclusion from different angles: foundation single-cell models do not beat trivial baselines on post-perturbation prediction, and the standard benchmarks are unable to discriminate model quality. Cite both, plus Ramakrishnan 2025, as "the field is realizing the standard benchmarks are inadequate." ConfPert positions as the constructive response: a distributional + calibration framework that (1) replaces inadequate Pearson Delta with proper distributional metrics (energy distance, MMD, per-cell Wasserstein), and (2) provides finite-sample calibrated coverage guarantees that no current benchmark offers. The Csendes critique of low-variance datasets motivates ConfPert's emphasis on dataset selection and stratified evaluation.

Specific transferable observations:
- Train Mean as a mandatory baseline: include in every ConfPert leaderboard.
- Pearson Delta only meaningful relative to a non-trivial baseline; raw Pearson useless.
- PEX-only evaluation is a known gap. ConfPert should report CEX and PEX separately.
- Intra-dataset correlation as a dataset diagnostic. Report for each evaluation dataset.

## What we wrap

Not a model. A benchmark paper used as a positioning citation. No code to wrap. The GitHub repo (https://github.com/turbine-ai/PerturbSeqPredBenchmark) hosts evaluation scripts that may be useful as reference for re-running their pseudo-bulk Pearson Delta on ConfPert outputs to demonstrate parity with prior leaderboards.

## Failure modes / caveats

Documented by Csendes:
- Low perturbation-specific variance in Adamson and Replogle RPE1 makes models hard to discriminate.
- Single-cell noise may obscure perturbation signal; pseudo-bulk performs comparably to single-cell methods, undercutting the rationale for single-cell foundation models.
- PEX-only split is incomplete.
- Foundation model embeddings, when repurposed as RF features, underperform GO annotation features. Verbatim from the paper: "perturbation effects ... may be better represented by functional categories such as GO terms than by the gene regulatory networks primarily learned by foundation models."
- The benchmark covers only four datasets, all CRISPRi/CRISPRa knockouts, no chemical perturbations, no dose response, no time course.
- Metrics are all aggregate Pearson; no distributional, calibration, or uncertainty metrics.

Caveats not addressed by Csendes that ConfPert should:
- No miscalibration analysis. The paper shows point predictions are bad but does not ask whether predicted distributions cover ground-truth single-cell distributions.
- No quantification of when foundation models would be expected to help (e.g., low-data regimes, transfer settings).

## Code URLs

- GitHub: https://github.com/turbine-ai/PerturbSeqPredBenchmark
- bioRxiv preprint: https://www.biorxiv.org/content/10.1101/2024.09.30.615843v1
- BMC Genomics article: https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-025-11600-2
- PMC: https://pmc.ncbi.nlm.nih.gov/articles/PMC12016270/

## Verbatim quotes worth keeping

1. "Surprisingly, we found that even the simplest baseline model - taking the mean of training examples - outperformed scGPT and scFoundation."

2. "basic machine learning models that incorporate biologically meaningful features outperformed scGPT by a large margin."

3. "we identified that the current Perturb-Seq benchmark datasets exhibit low perturbation-specific variance, making them suboptimal for evaluating such models."

4. "the Pearson correlation values between raw gene expression profiles are strongly influenced by the baseline expression magnitudes of different genes, so we did not consider these metrics meaningful."

5. "Moving forward, more rigorous and meaningful benchmarks that include higher variance and incorporate diverse datasets are needed to properly assess the applicability of machine learning models in post-perturbation prediction tasks."

6. "Our findings also raise questions about the utility of single-cell RNA-seq data for post-perturbation RNA-seq prediction ... the baseline models we used operated on pseudo-bulked data and performed comparably or better than foundation models."
