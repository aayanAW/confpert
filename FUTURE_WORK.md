# ConfPert — Deferred Future Work

Items not in current scope; tracked for next iteration. Last updated 2026-05-25.

---

## 1. CellFM heavyweight predictor (F3 family completion)

**Status:** Deferred. Documented as N/A per pre-reg v2 §1.6 missing-cell disclosure.

**Blocker:** No working HF mirror found 2026-05-25.
- `biomed-AI/cellfm-800M` → 404
- `biomap-research/CellFM` → 404
- Official repo is MindSpore-native (not torch). MindSpore install on Modal A100 with CUDA 12.x is the documented integration risk per `predictors_v2_stubs.CELLFM.notes`.

**Resolution paths:**
- Watch HF for re-uploads; check periodically.
- Invest 1-2 days porting MindSpore weights to PyTorch via `mindspore.load_checkpoint` + manual key remap.
- Wait for a community PyTorch port (none found 2026-05).

**Effort if mirror appears:** Wrapper code is ready (see `scripts/modal_launch.py::_cellfm_calibrate_fn`); replace `NotImplementedError` with embedding-baseline pattern (mirrors Geneformer).

---

## 2. McFaline-Figueroa lightweight sweep

**Status:** Deferred. Documented as N/A per pre-reg v2 §1.6.

**Blocker:** `mcfaline_figueroa.h5ad` assembled on Modal volume (~32 GB) from the 4 GSM sci-Plex-GxE samples, but the loader rejects it: `Only 0 vehicle/DMSO cells found`.

**Root cause:** Drug labels (DMSO, treatments, doses) live in the per-GSM `hashTable.out.txt.gz` files, NOT in `sample_meta`. Current assembler at `scripts/modal_launch.py::_mcfaline_assemble_fn` parses `sample_meta` field 3 as the drug column, which yields strings like `sci3_CRISPRi_kinome` (pathway, not drug).

**Fix needed:** Add `_merge_sciplex_hashtable` step to the assembler:
1. Read `<GSM>_hashTable.out.txt.gz` (per-cell hash counts).
2. Demux cells via dominant hash → drug treatment from `<GSM>_hash_sample_sheet.txt.gz`.
3. Attach `adata.obs["drug"]`, `adata.obs["dose"]` columns.

**Effort:** ~1-2 hours; pattern is documented in sci-Plex-GxE github (cole-trapnell-lab/sci-Plex-GxE/blob/main/notebooks/01_hash_demux.ipynb).

---

## 3. OSF mirror upload + DOI insertion

**Status:** Pending.

**Blocker:** Requires user OSF account + manual click-through. Handoff step 12.

**Resolution path:**
1. Follow `paper_neurips_dnb_2026/osf_upload_manifest.md` upload list.
2. Paste returned DOI into `paper_neurips_dnb_2026/preregistration_v2.yaml` field `lock.osf_doi`.
3. Re-emit hashes: `python -m confpert.cli prereg emit-hashes --prereg paper_neurips_dnb_2026/preregistration_v2.yaml`.
4. Commit + re-verify dry-run.

---

## 4. Remote push

**Status:** Pending.

**Blocker:** No `origin` remote configured.

**Resolution path:** User adds remote URL, then `git push -u origin main`.

---

## 5. NeurIPS D&B 2026 sty swap

**Status:** Placeholder.

**Current:** `paper_neurips_dnb_2026/tex/main.tex` uses `\usepackage{article} + geometry` until ICLR 2027 D&B CFP opens (~July 2026).

**Resolution path:** When CFP opens, swap preamble to official `neurips_data_2026.sty` (or whatever the official ICLR 2027 D&B template is named). Re-verify `\textwidth`, citation style, font.

---

## 6. scgpt full perturbation fine-tune (Phase 2B.3)

**Status:** Embedding-baseline mode shipped (Csendes 2025 precedent). Full fine-tune deferred.

**Why not done:** Full `scgpt.Trainer.train_perturb` integration costs $50/dataset × 6 datasets = $300 of $3K cap. Risk of multi-day failure across 4 wrappers. Embedding-baseline meets pre-reg §1.6 contract.

**Effort to add:** 2-3 days end-to-end including GenePerturbationDataset schema integration, gene-tokenization, batch-id handling, control cell pairing.

---

## 7. Test snapshot freshness

**Status:** All 126 tests pass on current data. Phase 1 fixture (`tests/fixtures/phase1_results.json` from commit dbcc769) is frozen.

**Future maintenance:** If pre-reg schema changes (new hypothesis added or threshold tweaked), the fixture-based snapshot expected values in `test_prereg_full_verify.py` and `test_prereg_ablation.py` need re-computation. Generator script in `scripts/` would help (currently manual `python -c "..."` invocations).

---

## 8. Style audit improvements

**Status:** Strict mode passes 0 violations. Bold-shouting false-positive filter added (table-cell + comment context).

**Future improvements:**
- Add `equation` and `align` environment context skip for any future rule that fires inside math mode.
- Add per-section budget overrides (intro vs methods have different prose density).
- Optional: switch from regex to LaTeX AST parser (`pylatexenc`) for cleaner context detection.
