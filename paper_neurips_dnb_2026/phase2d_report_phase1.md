# Phase 2D Analysis Report

- Timestamp: `2026-05-25T09:48:42+00:00`
- YAML: `paper_neurips_dnb_2026/preregistration_v2.yaml`
- Results: `baselines/results.json`

## Verifier report

- Lock check passed: **True**
- K2 v2 disposition: **`triple_null_promote_ablation_panel`**

| Hypothesis | Status | Observed |
|---|---|---|
| H1_capacity | FAIL | n_datasets_passing=2, datasets_passing=['replogle_rpe1', 'tahoe'], per_datase... |
| H1b_data | FAIL | n_datasets_passing=0, datasets_passing=[], per_dataset={'adamson': {'rho': na... |
| H2_cell_line_covariate | FAIL | f_pvalue_cell_line_context=0.11450855606238489, eta_squared_cell_line_context... |
| H2b_corpus_diversity | FAIL | n_datasets_passing=0, datasets_passing=[], per_dataset={'adamson': {'rho': -0... |
| H3_architecture_family | FAIL | kw_statistic=5.435623475323703, kw_pvalue=0.06601906350602582, cliffs_delta=0... |

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
| ahlmann_bilinear_ridge | norman | bimodality_mismatch | 0.050 | 100 | 0.980 | [0.930, 0.994] |
| ahlmann_bilinear_ridge | norman | bimodality_mismatch | 0.100 | 100 | 0.980 | [0.930, 0.994] |
| ahlmann_bilinear_ridge | norman | bimodality_mismatch | 0.200 | 100 | 0.840 | [0.756, 0.899] |
| ahlmann_bilinear_ridge | norman | energy | 0.050 | 100 | 0.880 | [0.802, 0.930] |
| ahlmann_bilinear_ridge | norman | energy | 0.100 | 100 | 0.880 | [0.802, 0.930] |
| ahlmann_bilinear_ridge | norman | energy | 0.200 | 100 | 0.790 | [0.700, 0.858] |
| ahlmann_bilinear_ridge | norman | ks | 0.050 | 100 | 0.880 | [0.802, 0.930] |
| ahlmann_bilinear_ridge | norman | ks | 0.100 | 100 | 0.880 | [0.802, 0.930] |
| ahlmann_bilinear_ridge | norman | ks | 0.200 | 100 | 0.820 | [0.733, 0.883] |
| ahlmann_bilinear_ridge | norman | mmd_rbf | 0.050 | 100 | 0.970 | [0.915, 0.990] |
| ahlmann_bilinear_ridge | norman | mmd_rbf | 0.100 | 100 | 0.970 | [0.915, 0.990] |
| ahlmann_bilinear_ridge | norman | mmd_rbf | 0.200 | 100 | 0.870 | [0.790, 0.922] |

_…and 817 more cells in the JSON sidecar._