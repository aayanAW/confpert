# ConfPert compute estimate

**Budget cap:** $5,000 - $8,000 Modal/Lambda H100 spot, per user spec.
**Estimate target:** $3,000 - $5,000 with margin to $8,000.
**Locked:** before Phase 1 baseline freeze.

## Per-job estimate (Modal H100 spot pricing, ~$3.20/h on-demand, ~$1.50/h spot)

| # | Job | GPU-h estimate | Cost @ $1.50/h | Notes |
|---|-----|---------------:|---------------:|-------|
| **K1 lightweight predictors (laptop CPU)** | | | | |
| 1 | Mean baseline on 4 datasets x 3 splits x 4 noise variants | 0 GPU | $0 | CPU-tractable, ~30 min total |
| 2 | Bilinear ridge per Ahlmann-Eltze 2025 on 4 datasets x 3 splits | 0 GPU | $0 | CPU-tractable, ~1 h total |
| 3 | scGen retrain on 4 datasets | 4 GPU-h | $6 | Small VAE, fits CPU but H100 is faster |
| 4 | sVAE+ / SAMS-VAE retrain on 3 datasets | 12 GPU-h | $18 | Posterior sampling adds compute |
| **K1 medium predictors (Modal H100)** | | | | |
| 5 | CPA on 4 datasets (cpa-tools, retrain per dataset) | 16 GPU-h | $24 | scvi-tools backbone |
| 6 | biolord on 4 datasets (retrain per dataset) | 24 GPU-h | $36 | Latent optimization is slow |
| 7 | GEARS in uncertainty mode on 4 datasets | 12 GPU-h | $18 | cell-gears v0.1.2, uncertainty head exposed |
| **K1 heavy predictor (Modal H100)** | | | | |
| 8 | STATE inference using SE-600M public checkpoint + ST per dataset | 30 GPU-h | $45 | No retrain; inference + calibration only |
| **K2 (calibration-vs-capacity / -vs-data analysis)** | | | | |
| 9 | Aggregation + Spearman + permutation tests | 0 GPU | $0 | CPU-tractable, ~1 h |
| **K3 PRISM discovery pipeline** | | | | |
| 10 | Tahoe-100M subset preprocessing | 8 GPU-h | $12 | Subset to a tractable cell-line / drug subset, not full 100M cells |
| 11 | STATE fine-tune on Tahoe-subset for K3 | 24 GPU-h | $36 | Fine-tune to PRISM-relevant compounds |
| 12 | K3 enrichment pipeline (BH on Hallmark + DepMap) | 0 GPU | $0 | CPU-tractable, ~2 h |
| **Re-runs and debugging buffer (per HetPert pattern)** | | | | |
| 13 | Failed Modal launches and retries (4x attempts pattern from HetPert GEARS) | ~40 GPU-h | $60 | Empirical from prior session: 4 attempts per heavy predictor |
| 14 | Cross-cell-line K562 -> RPE1 calibration runs | 20 GPU-h | $30 | One per medium predictor |
| **Cell-Eval audit + plug-in** | | | | |
| 15 | Cell-Eval clone + grep + plug-in registration | 0 GPU | $0 | Engineering, ~30 min |

**Subtotal Modal:** ~$285 at spot, ~$610 at on-demand.

## Spot vs on-demand contingency

H100 capacity is transiently saturated per the prior HetPert session (CauseFlow Replogle v2 had to swap H100 -> A100 due to capacity). Plan for two scenarios:

- **Best case (spot, no queue):** ~$300 total. Well under budget.
- **Median case (50% spot + 50% on-demand):** ~$450.
- **Worst case (full on-demand + 2x retries):** ~$1,200.

All three are inside the $3,000 - $5,000 estimate target, with the $5,000 - $8,000 cap leaving 4x to 16x margin for unexpected debugging and re-runs.

## What is NOT included in the budget

- Wet-lab validation: out of scope per user directive.
- scDFM, CellFlow, Squidiff, PerturbDiff, CellOT, Unlasting in K1 sweep: cited as related work, not benchmarked. Adding them would push the budget by approximately $300 each due to per-dataset retraining.
- Full Tahoe-100M training: 100M cells is too large for our budget. We use a subset (top 50 drugs x 10 cell lines x ~500K cells ~= ~10% of full).
- Adamson 2016 in main K1 sweep: appendix only, runs on laptop CPU at ~$0.

## Disk and storage

- Modal volume `causeflow-artifacts` already provisioned (used by HetPert).
- 4 h5ads on /Volumes/STORAGE flash drive (~3.2 GB) carry over from HetPert.
- Tahoe-100M subset: ~8 GB additional (download to flash, symlink to repo).
- Cell-Eval clone: ~50 MB.

## Time estimate

- Phase 1 framework + library: 1-2 weeks of part-time work (1 person).
- Phase 1 baseline locking: 2-3 days once framework lands.
- Phase 1 K1 sweep complete: 1 week of Modal + analysis.
- Phase 1 K2 + K3 analysis: 1 week.
- Total to a workshop-grade draft: ~5-6 weeks.

## Audit trail

This estimate locks before Phase 1 baseline-freeze. Any line item that exceeds estimate by 2x triggers a session pause and re-estimate, not silent overshoot. Final actual cost will be appended to this document at Phase 1 close.
