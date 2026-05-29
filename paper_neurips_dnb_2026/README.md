# Phase 2 Paper Directory — ICLR / NeurIPS D&B 2027 Submission

**Status:** Phase 2A in progress.
**Working title:** *Capacity is not Calibration: A Pre-Registered Distribution-Aware Benchmark of Single-Cell Perturbation Predictors*
**Target venue:** ICLR 2027 D&B (primary) → NeurIPS 2027 D&B (auto-fallback). Per `PHASE_2_PLAN.md` Section 1.

## Directory layout

```
paper_neurips_dnb_2026/
├── README.md                       # this file
├── preregistration_v2.md           # human-readable K2 v2 pre-registration (Phase 2A.2)
├── preregistration_v2.yaml         # machine-readable pre-reg schema (Phase 2A.3)
├── preregistration_v2.sha256       # SHA-256 fingerprints of .md + .yaml (Phase 2A.4)
├── osf_upload_manifest.md          # OSF-ready upload bundle instructions (Phase 2A.4)
├── tex/                            # .tex sources (populated Phase 2E)
├── figures/                        # .pdf figures (populated Phase 2D/2E)
├── bib/                            # refs.bib (populated Phase 2E)
└── datasheet.md                    # Gebru 2021 datasheet (Phase 2A bonus)
└── model_card.md                   # Mitchell 2019 model card (Phase 2A bonus)
```

## Template choice — DEFERRED

The official ICLR 2027 D&B `.sty` template will be pulled when the CFP opens (typically July 2026). Until then, `tex/` is empty. The `main.tex` should NOT be created yet — premature template choice was a v1 mistake.

When CFP opens:
1. Download `iclr2027_dnb.sty` (or NeurIPS 2027 D&B if pivoting)
2. Place in `tex/`
3. Create `tex/main.tex` from template
4. Follow Section 5.1 of `PHASE_2_PLAN.md` for section structure

## Section structure (planned, per `PHASE_2_PLAN.md` Section 5.1)

1. Introduction — lead with K2 cell-line covariate finding; drop "saturation" motivation
2. Related work — honest acknowledgement of 2025-26 multivariate CP cluster + parallel benchmarks
3. Benchmark design — datasets, predictors, splits, discrepancies, pre-reg protocol
4. K1 calibration sweep — main result table + Pareto figure
5. K2 cell-line covariate analysis — H1/H1b/H2/H2b/H3 + power analysis. **EMPIRICAL HEADLINE.**
6. K3 downstream PRISM utility — calibrated signatures recover more BH-FDR pathways
7. Pre-registration as transferable methodology — template, dispositions, ablation panel. **METHODOLOGICAL HEADLINE if K2 v2 lands null.**
8. Library, Cell-Eval plug-in, reproducibility
9. Limitations and discussion
10. Appendix: Datasheet + Model Card + per-predictor + per-dataset detail

## Cross-references

- Pre-registration v1 (Phase 1, locked): `../preregistration.md`
- Phase 2 plan: `../PHASE_2_PLAN.md`
- Phase 1 paper sources (frozen): `../paper_ai4science/`
- Library: `../src/confpert/`
- Results: `../baselines/results.json`
