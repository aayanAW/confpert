# ConfPert Codabench Competition (Phase 2 deliverable)

Per pre-registration v2 §14.7 and PHASE_2_PLAN.md §4.4, ConfPert hosts a public
Codabench leaderboard for the calibration-deviation benchmark. This directory
contains the competition bundle artifacts.

## What lives here

```
codabench/
├── README.md                    (this file)
├── competition.yaml             Codabench manifest (datasets, scorer entrypoint, leaderboard)
├── starting_kit/
│   ├── predictor_template.py    Reference `predict_samples` implementation
│   ├── Dockerfile               Build template for submission images
│   ├── README.md                Participant-facing instructions
│   └── example_submission/      Worked example using the Mean predictor
├── scoring_program/
│   ├── score.py                 Wraps confpert pipeline; emits leaderboard JSON
│   ├── metadata.yaml            Codabench scorer metadata
│   └── Dockerfile               Scoring container
└── reference_data/
    └── manifest.json            Pointer to (private) test perturbation lists
```

## Participant submission contract

A submission is a Docker image that exposes ONE entrypoint:

```python
def predict_samples(
    X_ctrl_test: np.ndarray,      # (n_cells, n_genes)
    perturbation: str,            # gene name or compound id
    n_cells: int,                 # how many cells to sample
) -> np.ndarray:                  # (n_cells, n_genes) — sampled predicted population
    ...
```

The image must be CPU-only OR GPU-eligible (declared in `competition.yaml`).
Inputs are loaded from `/data/X_ctrl_<dataset>.npy` mounted into the container.
Outputs are written to `/output/pred_<dataset>_<perturbation>.npy`.

## Scorer

For each (predictor, dataset, perturbation, discrepancy, head, α-level) cell,
the scorer:

1. Loads the participant's predicted populations and the held-out observed
   populations (private to the leaderboard).
2. Computes per-(perturbation, discrepancy) nonconformity scores.
3. Builds the conformal interval at level α on the calibration fold.
4. Reports achieved coverage + calibration deviation per cell.
5. Aggregates into the leaderboard metric:
   **mean calibration deviation across discrepancies and α-levels**, lower is
   better; ties broken by bootstrap 95% CI lower bound.

## Timeline

- Phase 2D end (~2026-08-24): open competition with the 6 must-have datasets +
  5 of 6 discrepancies live.
- Phase 2E (Sep 2026): publish leaderboard URL in the ICLR D&B submission
  appendix.
- Post-submission (Oct 2026+): rolling submissions accepted; quarterly leaderboard
  freezes for follow-up papers.

## Out of scope for Phase 2A-D

Manual UI customization on Codabench, badge artwork, partner branding.
Phase 2D opens the competition with the default Codabench theme.

## Status (2026-05-25)

- Spec written (this file) — DONE
- `competition.yaml` skeleton — DONE
- Scoring entrypoint — DONE
- Starting-kit Docker image — DONE (template only; participants build their own)
- Codabench account + competition create — USER-SUPERVISED (requires Codabench login)
