# Phase 2D Analysis Report

- Timestamp: `2026-05-25T23:38:56+00:00`
- YAML: `paper_neurips_dnb_2026/preregistration_v2.yaml`
- Results: `baselines/results.json`

## Verifier report

- Lock check passed: **True**
- K2 v2 disposition: **`double_pass_no_corpus`**

| Hypothesis | Status | Observed |
|---|---|---|
| H1_capacity | PASS | n_datasets_passing=4, datasets_passing=['frangieh', 'norman', 'replogle_rpe1'... |
| H1b_data | FAIL | n_datasets_passing=0, datasets_passing=[], per_dataset={'adamson': {'rho': na... |
| H2_cell_line_covariate | PASS | f_pvalue_cell_line_context=0.006351499458553939, eta_squared_cell_line_contex... |
| H2b_corpus_diversity | FAIL | n_datasets_passing=0, datasets_passing=[], per_dataset={'adamson': {'rho': 0.... |
| H3_architecture_family | PASS | kw_statistic=31.152894956205465, kw_pvalue=1.71883747013994e-07, cliffs_delta... |

## Pre-registration ablation panel

_Skipped._

## Power analysis summary

```json
{
  "H1_capacity_per_dataset_phase2": "PowerResult(test_name='spearman', n=14, alpha=0.05, target_effect=0.5, estimated_power=0.411, min_detectable_effect_at_80=0.7225247192382811, notes='Monte-Carlo n_simulations=2000, seed=42; UNDERPOWERED: estimated power 0.41 < 0.80 at target rho 0.5')",
  "H1b_data_cross_dataset_phase2": "PowerResult(test_name='spearman', n=10, alpha=0.05, target_effect=-0.5, estimated_power=0.254, min_detectable_effect_at_80=0.8237905883789063, notes='Monte-Carlo n_simulations=2000, seed=42; UNDERPOWERED: estimated power 0.25 < 0.80 at target rho -0.5')",
  "H2_anova": "PowerResult(test_name='two_way_anova', n=120, alpha=0.0167, target_effect=0.1, estimated_power=0.8120158202642977, min_detectable_effect_at_80=0.09776043487433836, notes='FTestAnovaPower with k_groups=3')",
  "H2b_spearman_per_dataset": "PowerResult(test_name='spearman', n=14, alpha=0.0167, target_effect=-0.5, estimated_power=0.258, min_detectable_effect_at_80=0.7857186889648438, notes='Monte-Carlo n_simulations=2000, seed=42; UNDERPOWERED: estimated power 0.26 < 0.80 at target rho -0.5')",
  "H3_kruskal_wallis": "PowerResult(test_name='kruskal_wallis', n=120, alpha=0.0167, target_effect=0.5, estimated_power=0.923, min_detectable_effect_at_80=0.4374191284179687, notes='Monte-Carlo n_simulations=2000, k_groups=3, seed=42')"
}
```

## Bootstrap CI table (Wilson 95% binomial)

| predictor | dataset | score | alpha | n_perts_test | p̂ | 95% CI |
|---|---|---|---|---|---|---|
| ahlmann_bilinear_ridge | adamson | bimodality_mismatch | 0.050 | 16 | 0.750 | [0.505, 0.898] |
| ahlmann_bilinear_ridge | adamson | bimodality_mismatch | 0.100 | 16 | 0.750 | [0.505, 0.898] |
| ahlmann_bilinear_ridge | adamson | bimodality_mismatch | 0.200 | 16 | 0.750 | [0.505, 0.898] |
| ahlmann_bilinear_ridge | adamson | energy | 0.050 | 16 | 0.938 | [0.717, 0.989] |
| ahlmann_bilinear_ridge | adamson | energy | 0.100 | 16 | 0.938 | [0.717, 0.989] |
| ahlmann_bilinear_ridge | adamson | energy | 0.200 | 16 | 0.938 | [0.717, 0.989] |
| ahlmann_bilinear_ridge | adamson | ks | 0.050 | 16 | 0.875 | [0.640, 0.965] |
| ahlmann_bilinear_ridge | adamson | ks | 0.100 | 16 | 0.875 | [0.640, 0.965] |
| ahlmann_bilinear_ridge | adamson | ks | 0.200 | 16 | 0.875 | [0.640, 0.965] |
| ahlmann_bilinear_ridge | adamson | mmd_rbf | 0.050 | 16 | 0.875 | [0.640, 0.965] |
| ahlmann_bilinear_ridge | adamson | mmd_rbf | 0.100 | 16 | 0.875 | [0.640, 0.965] |
| ahlmann_bilinear_ridge | adamson | mmd_rbf | 0.200 | 16 | 0.875 | [0.640, 0.965] |
| ahlmann_bilinear_ridge | adamson | variance_ratio_dev | 0.050 | 16 | 0.938 | [0.717, 0.989] |
| ahlmann_bilinear_ridge | adamson | variance_ratio_dev | 0.100 | 16 | 0.938 | [0.717, 0.989] |
| ahlmann_bilinear_ridge | adamson | variance_ratio_dev | 0.200 | 16 | 0.938 | [0.717, 0.989] |
| ahlmann_bilinear_ridge | adamson | w1 | 0.050 | 16 | 1.000 | [0.806, 1.000] |
| ahlmann_bilinear_ridge | adamson | w1 | 0.100 | 16 | 1.000 | [0.806, 1.000] |
| ahlmann_bilinear_ridge | adamson | w1 | 0.200 | 16 | 1.000 | [0.806, 1.000] |
| ahlmann_bilinear_ridge | datlinger | ? | 0.050 | 1 | nan | (invalid) |
| ahlmann_bilinear_ridge | datlinger | ? | 0.100 | 1 | nan | (invalid) |
| ahlmann_bilinear_ridge | datlinger | ? | 0.200 | 1 | nan | (invalid) |
| ahlmann_bilinear_ridge | frangieh | bimodality_mismatch | 0.050 | 250 | 0.983 | [0.959, 0.993] |
| ahlmann_bilinear_ridge | frangieh | bimodality_mismatch | 0.100 | 200 | 0.955 | [0.917, 0.976] |
| ahlmann_bilinear_ridge | frangieh | bimodality_mismatch | 0.200 | 200 | 0.875 | [0.822, 0.914] |
| ahlmann_bilinear_ridge | frangieh | energy | 0.050 | 250 | 0.953 | [0.920, 0.973] |
| ahlmann_bilinear_ridge | frangieh | energy | 0.100 | 200 | 0.920 | [0.874, 0.950] |
| ahlmann_bilinear_ridge | frangieh | energy | 0.200 | 200 | 0.840 | [0.783, 0.884] |
| ahlmann_bilinear_ridge | frangieh | ks | 0.050 | 250 | 1.000 | [0.985, 1.000] |
| ahlmann_bilinear_ridge | frangieh | ks | 0.100 | 200 | 0.895 | [0.845, 0.930] |
| ahlmann_bilinear_ridge | frangieh | ks | 0.200 | 200 | 0.780 | [0.718, 0.832] |

_…and 1486 more cells in the JSON sidecar._