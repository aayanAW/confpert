# ConfPert

**Calibrated distributional uncertainty for single-cell perturbation predictors.**

ConfPert is a model-agnostic conformal wrapper that provides finite-sample marginal
coverage on six distributional discrepancies for cell-population perturbation
responses. It wraps any sample-producing predictor (CPA, scGen, biolord, sVAE+,
GEARS, scGPT, STATE, scDFM, ...) and ships:

1. **K1: a calibration benchmark and library.** Eight wrappable predictors,
   four datasets, three split types, six discrepancies, three nominal coverage
   levels. Pip-installable. Cell-Eval plug-in.
2. **K2: a pre-registered empirical claim** about how foundation-model calibration
   scales with parameter count and training-set size. Two hypotheses, both
   pre-registered.
3. **K3: a discovery pipeline** for drug-resistance subpopulations on PRISM,
   validated against DepMap CRISPR essentiality and Hallmark MSigDB pathways
   under BH-FDR correction.

## Pre-registration and reproducibility

All hypotheses, success criteria, baseline-freeze conditions, and reporting
requirements are git-stamped at commit
`c7046e4b4c08106acd52c37e7ea8410e126d2ae4` (timestamp `2026-05-03 01:30:12 -0400`).
See [preregistration.md](preregistration.md), [_synthesis.md](lit_notes/_synthesis.md),
and [compute_estimate.md](compute_estimate.md).

## Status

- **Phase 0** (lit reading): complete (35 paper notes in [lit_notes/](lit_notes/))
- **Phase 1** (initial benchmark): complete. 9 predictors × 5 datasets, K1/K2/K3 results in [baselines/results.json](baselines/results.json). Workshop paper submitted to ICML 2026 AI4Research 2026-05-14. Phase 1 paper frozen at [paper_ai4science/](paper_ai4science/). See [PHASE_1_PLAN.md](PHASE_1_PLAN.md).
- **Phase 2** (expanded archival paper): Phase 2A locked 2026-05-25; **Phase 2B/2C/2D/2E/2F scaffolding complete 2026-05-25** (session 6). Target ICLR 2027 D&B (Sept-Oct 2026). Working title: *Capacity is not Calibration: A Pre-Registered Distribution-Aware Benchmark of Single-Cell Perturbation Predictors*. See [PHASE_2_PLAN.md](PHASE_2_PLAN.md), [paper_neurips_dnb_2026/](paper_neurips_dnb_2026/), [paper_neurips_dnb_2026/preregistration_v2.md](paper_neurips_dnb_2026/preregistration_v2.md), and [paper_neurips_dnb_2026/tex/main.tex](paper_neurips_dnb_2026/tex/main.tex).

### Phase 2 CLI (verifier, ablation panel, Phase 2D runner, style audit)
```bash
# Verify the K2 v2 pre-registration (machine-enforced lock check)
python -m confpert.cli prereg verify \
    --prereg paper_neurips_dnb_2026/preregistration_v2.yaml \
    --results baselines/results.json \
    --dry-run

# Re-emit SHA-256 lock block (after any pre-reg edit)
python -m confpert.cli prereg emit-hashes \
    --prereg paper_neurips_dnb_2026/preregistration_v2.yaml

# Run pre-registration ablation panel (Pipeline A vs B counterfactual)
python -m confpert.cli prereg-ablate \
    --prereg paper_neurips_dnb_2026/preregistration_v2.yaml \
    --results baselines/results.json

# One-shot Phase 2D: verifier + ablation + power + bootstrap CIs
python -m confpert.cli phase2d \
    --prereg paper_neurips_dnb_2026/preregistration_v2.yaml \
    --results baselines/results.json \
    --repo-dir . \
    --out paper_neurips_dnb_2026/phase2d_report.json

# Editorial audit on the paper draft (em-dash budget, rule-of-three, etc.)
python -m confpert.cli style-audit paper_neurips_dnb_2026/tex/main.tex --strict
```

### Phase 2 paper draft (compiles via tectonic)
```bash
cd paper_neurips_dnb_2026/tex && tectonic -X compile main.tex --reruns 4
# -> main.pdf (~70 KB skeleton, fills in once Phase 2D results land)

# Auto-generate Phase 2 figures from a phase2d_report.json
python scripts/make_phase2_figures.py \
    --phase2d-report paper_neurips_dnb_2026/phase2d_report_phase1.json \
    --out paper_neurips_dnb_2026/figures/
```

See [CHANGELOG.md](CHANGELOG.md) for version history.

## Install

```bash
git clone https://github.com/alwaniaayan-png/confpert.git
cd confpert
pip install -e .[dev]
pytest -q
```

## Layout

```
confpert/
├── pyproject.toml
├── README.md
├── PHASE_1_PLAN.md
├── compute_estimate.md
├── preregistration.md
├── data/                       # symlinks to /Volumes/STORAGE flash drive h5ads
├── src/confpert/
│   ├── __init__.py
│   ├── metrics.py              # 6 distributional discrepancies + aggregator
│   ├── conformal.py            # 4 conformal heads + RCPS risk control
│   ├── predictors.py           # lightweight predictors + noise variants
│   ├── noise_models.py         # 4 pre-registered noise variants for point predictors
│   └── data/                   # 4 dataset loaders (Norman, Replogle K562/RPE1, Adamson)
├── tests/
│   ├── test_metrics.py
│   ├── test_conformal.py
│   └── test_predictors.py
├── scripts/
│   ├── run_k1_baseline.py      # K1 sweep harness
│   ├── cell_eval_audit.py      # Cell-Eval audit before metric implementation
│   └── modal_launch.py         # Modal --detach entrypoints for heavy predictors
├── baselines/
│   ├── results.json            # locked baseline results, SHA-256 in preregistration
│   └── cell_eval_audit.md      # Cell-Eval HEAD commit + present/missing list
└── lit_notes/                  # 35 paper notes + _synthesis.md
```

## Citation

Pre-publication. Internal project of Aayan Alwani, May 2026.
