# ConfPert Phase 1 Plan

**Locked references:** Phase 0 commit `c7046e4b4c08106acd52c37e7ea8410e126d2ae4` (preregistration.md, 35 lit notes, _synthesis.md). compute_estimate.md committed in this same Phase 1 prep batch.

**Realism statement:** ConfPert end-to-end is ~5-6 weeks of full-time work. Phase 1 goal is to ship the **framework library + locked K1 baselines on the laptop-CPU-tractable predictors + Modal-detached jobs for the heavy predictors + Cell-Eval audit + status report**. Partial Phase 1 completion, not full project completion.

## Acceptance criteria for tonight (concrete and measurable)

1. `confpert/` Python package builds and `pip install -e .` succeeds.
2. Unit tests pass: at minimum the metrics module and the conformal-head modules.
3. Cell-Eval audit committed at `baselines/cell_eval_audit.md` with HEAD commit hash.
4. K1 baselines locked for at least 3 lightweight predictors (Mean, Bilinear ridge per Ahlmann-Eltze, plus one of {scGen, GEARS-uncertainty, sVAE+}) on at least 1 dataset (Norman). Results saved with SHA-256 hash to `baselines/results.json`.
5. Modal launch scripts written and detached jobs running for at least 2 of {CPA, biolord, sVAE+, GEARS, STATE}.
6. Caffeinate process active throughout to prevent laptop sleep.
7. HANDOFF_CONFPERT_<date>.md status report committed describing what's done, what's running, what failed, what's next.
8. All work committed to git with explicit messages.

## Execution order

### Phase A: Pre-flight (~30 min, in progress)
- [x] Confirm preregistration.md committed (commit c7046e4b).
- [x] caffeinate -dimsu running (PID 80301, 10h budget).
- [x] compute_estimate.md committed.
- [ ] PHASE_1_PLAN.md committed.

### Phase B: Library scaffolding (~1.5 h)
- Create `src/confpert/` Python package with `__init__.py`, `metrics.py`, `conformal.py`, `predictors.py`, `noise_models.py`, `data/`, `cli.py`.
- Port HetPert's `metrics.py` (5 distributional metrics + bimodality coefficient + variance ratio) and add MMD-RBF as the sixth.
- Port HetPert's `conformal.py` (CDSplitConformal, EnergyConformal) and add MMDConformal, OTCPConformal, CauchoisSubgroupConformal.
- Port HetPert's `baselines.py` (Mean, Additive, Ridge, NoisyMean, OracleMixture).
- Add `AhlmannEltzeBilinearRidge` per their eq. 1, 3 (K=10 PCs, lambda=0.1).
- Add four noise variants for point predictors (no-noise, isotropic, per-gene marginal, full empirical covariance) per preregistration.
- Port HetPert's dataset loaders (Norman, Replogle K562, Replogle RPE1, Adamson) verbatim.
- Stub Tahoe-100M loader (defer actual data; record TODO).
- Write `pyproject.toml` for pip install.
- Write `tests/test_metrics.py`, `tests/test_conformal.py`, `tests/test_noise_models.py`.

### Phase C: K1 lightweight baselines (~2 h)
- Write `scripts/run_k1_baseline.py` that takes (dataset, predictor, split-type) and produces a `baselines/results.json` row.
- Run on Norman first (already on flash, fastest), 3 splits, 3 predictors (Mean, Bilinear ridge, NoisyMean from HetPert), 6 discrepancies, 3 alphas, 4 noise variants.
- SHA-256 hash the results.json after each row.
- Commit to git with hash recorded in preregistration.md.
- Repeat on Replogle K562, Replogle RPE1, Adamson if time permits.

### Phase D: Modal launches detached (~1 h setup, runs in background)
- Adapt HetPert's `modal_launch.py` to ConfPert's API.
- Three entrypoints: `cpa_calibrate`, `biolord_calibrate`, `gears_calibrate` (uncertainty mode), `state_calibrate`.
- Use `modal run --detach` for all four.
- Write polling script that pulls `causeflow-artifacts:/confpert_<predictor>_results.json` and merges into local `baselines/results.json`.
- Note: STATE 600M-checkpoint requires HF auth + Arc license terms; if blocked, defer to next session and document.

### Phase E: Cell-Eval audit (~30 min)
- Clone `ArcInstitute/cell-eval` at HEAD.
- Grep for `wasserstein|kolmogorov|mmd|bimodal|variance_ratio|energy` in `src/cell_eval/metrics/_impl.py`.
- Document commit hash + present/missing list in `baselines/cell_eval_audit.md`.
- Commit.

### Phase F: K2 prelim analysis on whatever is locked (~30 min)
- If at least 4 predictors locked across at least 2 datasets, compute Spearman rho for H1 and H1b on the partial data.
- Mark as "preliminary, full sweep pending Modal completion."
- Save to `results/k2_preliminary.json`.

### Phase G: K3 scaffold (~30 min, deferred for full run)
- PRISM data loader stub.
- Hallmark MSigDB loader (gmt file from broad).
- DepMap dependency loader.
- BH-FDR correction implementation.
- Pipeline scaffold without execution; needs Tahoe-100M trained model.

### Phase H: Status report and final commit (~30 min)
- HANDOFF_CONFPERT_2026-05-03_overnight.md describing every step taken, every failure, every command run, every artifact produced.
- Final git commit + tag `phase1-overnight-stop`.
- Modal app list dump.
- Caffeinate cleanup if natural session-end.

## Auto-debug loop policy

For any failing step:
1. Capture error log to `debug/<step>_<timestamp>.log`.
2. Apply known-issues mitigations from HetPert handoff (pandas .nonzero patch, scgpt flash-attn skip, ctrl-pert filter, etc.).
3. Retry once.
4. If still failing, document in HANDOFF as "blocked, manual review needed" and proceed to next step.
5. Never use destructive operations (rm -rf, git reset --hard, etc.) without explicit user authorization.

## What I will NOT do tonight

- Wet-lab validation: out of scope.
- Full Tahoe-100M training: too expensive overnight.
- scDFM / CellFlow / Squidiff / PerturbDiff in K1 sweep: per preregistration, cited only.
- Push to a remote: no remote configured for confpert/ yet.
- Modify preregistration.md: locked, cannot change retroactively.
- Add new theoretical objects: Phase 0 review locked the no-theory posture.
- Write the paper draft: defer to a daytime session with manuscript focus.

## Resource accounting

- Modal spend cap tonight: $200 (well within compute_estimate.md budget).
- Disk: /Volumes/STORAGE flash drive plus internal disk (9.7 GB free per HetPert handoff).
- Time: 8-10 hours wall clock from 01:30 EDT, ending ~10:00 EDT.
- Caffeinate: 10h budget started.

## Termination

Session ends naturally when either:
- Phase H (status report + final commit) lands successfully.
- A blocker prevents further productive work AND HANDOFF documents the state honestly.
- 10h caffeinate budget elapses (in which case I save state and terminate gracefully).

I will not extend caffeinate or self-impose additional time. The user expects a status report at wake-up.
