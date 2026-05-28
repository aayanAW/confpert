# Changelog

All notable changes to ConfPert documented here. Format per [Keep a Changelog](https://keepachangelog.com/en/1.1.0/). Versioning per [SemVer 2.0](https://semver.org/spec/v2.0.0.html).

---

## [Unreleased] — Phase 2 in progress

### Added (2026-05-25 session 6: Phase 2B/C/D/E/F scaffolding)
- `src/confpert/data/frangieh.py` — real Phase 2C loader (scperturb Zenodo h5ad)
- `src/confpert/data/datlinger.py` — real Phase 2C loader (replaces last session's wrong count_matrix path with scperturb DatlingerBock2017.h5ad)
- `src/confpert/data/schmidt.py` — real Phase 2C loader, assembles GEO cellranger 4-file bundle to h5ad on the fly
- `src/confpert/data/mcfaline_figueroa.py` — real Phase 2C loader (GSE225775 RAW.tar, ~12 GB; user-supervised extraction; resolves the prereg-v2 §1.2 "GSE_TBD" placeholder)
- `src/confpert/data/walker.py` — documented blocker (Walker GSE189574 URL unresolvable; 3-option resolution path surfaced in `NotImplementedError` with substitution candidates from scperturb LaraAstiasoHuntly2023)
- `src/confpert/heavyweight_helpers.py` — `safe_train_variance` + `sample_from_train_distribution` + `make_train_test_split_by_perturbation`; explicit train-fold enforcement prevents the scGPT v1 test-variance leak shape from recurring in future predictor wrappers
- `src/confpert/phase2d_runner.py` — one-shot Phase 2D pipeline (verifier + ablation + power + bootstrap CIs), CLI `confpert phase2d`
- `src/confpert/style_audit.py` — LLM-style audit (em-dash budget, rule-of-three lists, `\textbf{PASS}`, self-praising adjectives, "we introduce" repetition, etc.), CLI `confpert style-audit`
- `scripts/make_phase2_figures.py` — Phase 2 figure-generation script (Fig 2 K1 heatmap, Fig 3 K2 covariate scatter, Fig 4 Pareto)
- `paper_neurips_dnb_2026/tex/main.tex` — Phase 2 paper skeleton (10 sections + 5 appendices, compiles to 70 KB PDF via tectonic)
- `paper_neurips_dnb_2026/bib/refs.bib` — extended Phase 1 refs.bib with 2025-2026 multivariate CP cluster + foundation-model benchmark + dataset citations
- `paper_neurips_dnb_2026/codabench/` — Codabench leaderboard bundle: competition.yaml + starting_kit + scoring_program
- `paper_neurips_dnb_2026/phase2d_report_phase1.{json,md}` — first phase2d run on Phase 1 results.json (847 cells, Wilson 95% binomial CIs)
- `paper_neurips_dnb_2026/figures/{fig2_k1_heatmap,fig3_k2_covariate,fig4_pareto}.{pdf,png}` — auto-generated Phase 1 versions
- 40 new unit tests across `tests/test_phase2_loaders.py` (12) + `tests/test_heavyweight_helpers.py` (12) + `tests/test_phase2d_runner.py` (3) + `tests/test_style_audit.py` (10) + extended `tests/test_phase2_dataset_stubs.py` (+3)
- Heavyweight predictor wrappers (Modal stubs): Geneformer, scFoundation, CellFM in `scripts/modal_launch.py` (each: IMAGE + prestage entrypoint + calibrate-fn NotImplementedError with Phase 2B.2 TODO + heavyweight_helpers references)

### Changed
- `src/confpert/data/_phase2_stubs.py` — Phase 2C URL resolution: Datlinger/Frangieh switched to scperturb Zenodo, Schmidt to GEO cellranger bundle, McFaline-Figueroa to GSE225775 RAW.tar, Walker marked TBD with documented blocker; `loader_status` field added to spec
- `scripts/modal_launch.py` — `_phase2_dataset_prestage_fn` rewritten with per-dataset manifest (zenodo_h5ad / geo_cellranger_bundle / geo_raw_tar / blocked) replacing last session's blind URL guessing; multi-file Schmidt + 12 GB McFaline tar paths supported
- `src/confpert/cli.py` — added `phase2d` and `style-audit` subcommands
- `src/confpert/__init__.py` — version 0.2.0rc1 → 0.2.0rc2
- `pyproject.toml` — version 0.2.0rc1 → 0.2.0rc2; pytest `integration` marker registered

### Added (2026-05-25 session 1-5: Phase 2A lock)
- `paper_neurips_dnb_2026/` directory for Phase 2 paper artefacts
- `paper_neurips_dnb_2026/preregistration_v2.md` — full K2 v2 H2/H2b/H3 hypothesis family + falsification-intervention conditional design + K3 v2 sharpenings + pre-reg ablation panel contingency
- `paper_neurips_dnb_2026/preregistration_v2.yaml` — machine-readable schema, authoritative for verifier
- `paper_neurips_dnb_2026/preregistration_v2.sha256` — fingerprint file
- `paper_neurips_dnb_2026/datasheet.md` — Gebru 2021 template, 7 sections
- `paper_neurips_dnb_2026/model_card.md` — Mitchell 2019 template, 8 sections
- `paper_neurips_dnb_2026/osf_upload_manifest.md` — instructions for OSF pre-reg mirror
- `paper_neurips_dnb_2026/prereg_v2_dryrun_20260525.json` — first lock-time dry-run verification (Lock check passed: True)
- `paper_neurips_dnb_2026/README.md` — directory documentation
- `src/confpert/prereg.py` — pre-registration verifier module (`emit-hashes` + `verify --dry-run` + `verify` full modes; statistical-test implementations are stubs to be completed in Phase 2D)
- `src/confpert/prereg_ablation.py` + `src/confpert/bootstrap.py` + `src/confpert/power.py` — analysis modules per pre-reg v2 §3-§4
- `src/confpert/predictors_v2_stubs.py` + `src/confpert/data/_phase2_stubs.py` — Phase 2 predictor + dataset metadata
- `PHASE_2_PLAN.md` — full multi-venue pivot plan with ultrathink refinements
- `lit_notes/scoop_watch.md` — weekly literature scan log for scoop defense
- `CHANGELOG.md` — this file

### Changed (session 1-5)
- `src/confpert/__init__.py` — version bumped 0.1.0 → 0.2.0rc1
- `src/confpert/cli.py` — argparse dispatcher with `prereg verify` and `prereg emit-hashes` subcommands

### Phase 1 paper status (unchanged)
- `paper_ai4science/` frozen at AI4Research workshop submission (2026-05-14). Subsequent improvements (intro reframe, K2 honest reinterpretation, abstract trim, STATE-Tahoe unblock + K2 update) applied in this session 2026-05-24, recompiled into PDFs ready for any backup-venue submission.

---

## [0.1.0] — 2026-05-03 (Phase 1 lock-in)

### Added
- Phase 1 ConfPert library with 4 modules (metrics, conformal, predictors, noise_models)
- 4 dataset loaders (Norman, Replogle K562, Replogle RPE1, Adamson)
- 8 wrapped predictors + 4 noise variants (per `preregistration.md`)
- 33/33 unit tests
- Modal launcher (`scripts/modal_launch.py`)
- Cell-Eval plug-in (`confpert.cell_eval_plugin`)
- K1 sweep + K2 hypothesis testing + K3 PRISM downstream
- Pre-registration (`preregistration.md`, git commit `c7046e4b`, 2026-05-03)
