# ConfPert — Context for Ultraplan (2026-05-31)

> Paste this into the running ultraplan web session, or commit it so the remote
> session reads it. **Without this, ultraplan is planning against stale git state**
> (all fixes below are uncommitted working-tree changes; the public repo is at the
> old squash commit).

## 0. Critical orientation (read first)

- **Reviewed paper = `paper_ai4science/main.tex`** (ICML 2026 style, 12pp PDF = ~5pp
  body + ~7pp appendix, 484 lines). NOT `main_short.tex` / `main_paper.tex` /
  `supplementary.tex` (post-submission reformats) and NOT `paper/main.tex`
  (NeurIPS-style Phase-2 reframe). Proven by section/table/figure numbering matching
  the verbatim reviews (§8 Limitations, Table 3 = cross-cell-line, Fig 4(b) = K3 ratio)
  + the ICML-only "Title Suppressed" header.
- **git root is `confpert/`, NOT `~/FLM4S`.** (This is why ultraplan failed to launch
  the first time — cwd `~/FLM4S` is not a git repo.)
- **All session fixes are UNCOMMITTED.** Public remote `github.com/aayanAW/confpert`
  is at the pre-fix state. Commit before relying on any remote tool, or the plan will
  duplicate completed work.
- **Decision already in hand:** Accept (Poster) at ICML 2026 AI4Research. Reviews:
  yhtC=3 (clear-reject, integrity), workshop-fit reviewer=4 (reject, fit/novelty),
  Unm6=7 (accept, validity cons), usGP=7 (accept, breadth).
- **Current gate state (verified this session):** compile exit 0 · "Title Suppressed"
  in PDF = 0 · literal `(?)` = 0 · undefined cites = 0 · braces balanced · body
  code-path tokens = 0 · `style-audit --strict` 0 block / 0 over-budget · 12pp.

## 1. Issues found (by source) and status

### From reviewer yhtC (integrity/presentation) — ALL FIXED
| Issue | Status | Where |
|---|---|---|
| Fabricated Csendes title ("Train-mean as a control…") | ✅ FIXED (real title; 22-entry bib re-verified vs primary sources) | `refs.bib:11` |
| 3 broken `(?)` cites (§2.1 energy, §3.1 Tahoe, §8 distr-free) | ✅ FIXED (szekely2013 / tahoe2025 / vovk2012); 0 undefined | bib + body |
| Fig 3 labels stacked 3/4 panels | ✅ FIXED (regenerated, de-overlapped) | `fig_k2_h1.pdf` |
| Fig 4(b) overlapping x-ticks | ✅ FIXED (dropped log scale) | `fig_k3_robustness.pdf` |
| Body dense w/ code identifiers (§5/6/8) | ✅ FIXED (all 4 named tokens out of body; §Reproducibility L207 command strings → plain prose + appendix ref) | `main.tex:207` |
| "Title Suppressed Due to Excessive Size" header | ✅ FIXED (re-assert `\@icmltitlerunning` after `\icmltitle`) | `main.tex:55` |

### From rating-4 reviewer (fit/novelty/style)
| Issue | Status |
|---|---|
| Off-topic for testing workshop | ⚠️ ACK — unfixable for this venue; reviewer says target single-cell / benchmarks track → **Track B** |
| Modified style file → ICML desk-reject | ✅ LIKELY FIXED (= the Title-Suppressed symptom; sty is official 2025-10-29, unmodified). Verify no separate workshop `.sty` at camera-ready |
| "reads very much in an LLM style" | ⚠️ PARTIAL — style-audit 0/0 (necessary, not sufficient); needs `humanizer` + read-aloud human pass |
| Methodological novelty unclear / off-the-shelf | ❌ NOT fixed by edits — real work, **Track B only** (see §3) |
| Reviewer's own steer: pre-registration as transferable methodology | 🎯 ADOPT — light reframe Track A, first-class section Track B |

### From Unm6 (validity/clarity)
| Issue | Status |
|---|---|
| Cross-cell-line guarantee not flagged; Energy→0 should be foregrounded | ✅ FIXED (tab:xcl caption caveat "not coverage-guaranteed" + body §8 Energy→0 sentence) |
| STATE-600M missing from Tahoe | ⚠️ PARTIAL — stated plainly in body; re-run = author/compute (Track B) |
| No glossary | ✅ FIXED (Appendix A glossary table) |
| Main 5pp / appendix 12pp inversion | ⚠️ FLAGGED — blocked by 4pp body limit; fix at D&B venue (Track B) |

### From usGP (breadth)
| Issue | Status |
|---|---|
| Broaden perturbation types / baselines | ⚠️ future work → Track B = Phase 2's 14 predictors × 10 datasets |

### Found by the audit (NOT raised by reviewers) — FIXED
| Issue | Status |
|---|---|
| "weighted CP exposed in the library" claimed ×3 — **VERIFIABLY FALSE** (no weighted-CP code exists; `conformal.py` has 5 heads only) | ✅ FIXED → "future work / not implemented" at L85/L203/L256 |
| Fig 3 stale + self-contradictory ("0/4 pass" vs Table 2's "2/5") | ✅ FIXED (regenerated from `baselines/results.json`, matches Table 2 to 6dp, no fabricated points) |

### Audit self-corrections (do NOT re-introduce these as "issues")
- "stray brace `}`" in Discussion — **does not exist** (braces balanced 699/699); was a hypothesis, never real.
- Wei 2025 / TISSUE 2023 / VCC citations — were the auditor's scoop-defense suggestions, **not reviewer requests**. Optional Track-B related-work, not blockers.
- Earlier "main.tex is corrupted" — false; the file was clean.

## 2. Still open — Track A (finish THIS poster; ~half day, $0)
1. **Confirm `results/k2_preliminary.json` is canonical K2 source.** Fig 3 + Table 2
   depend on it; `results/k2_state_n8.json` is STALE (Tahoe n=7). — **USER INPUT.**
2. **LLM-style pass** — run `humanizer` on body prose + read-aloud; re-run style-audit;
   diff every number after.
3. **Contribution-paragraph reframe (L68)** — lead with pre-registration protocol +
   sVAE+ output-head design rule; conformal heads named as off-the-shelf; K562 finding
   = explicitly post-hoc secondary. Author's words.
4. **Attribute VCC quote (§1, ~L68)** — add Roohani 2025 *Cell* cite for "almost all
   models performed worse than baseline on MAE" (currently unattributed). 1 bib + 1 cite.
5. **Verify workshop camera-ready template** (10-min; no separate `.sty`).
6. **Re-run gate** (compile / Title-Suppressed / (?) / undefined / style-audit) + eyeball both figs.
7. **(Optional) blind→named** if camera-ready is non-blind: `\usepackage{icml2026}` →
   `[accepted]` + author block.
8. **(Optional polish)** move 2 borderline body tokens (`confpert.cell_eval_plugin`,
   `ArcInstitute/cell-eval`) to appendix.

## 3. Still open — Track B (archival paper; folds into existing `PHASE_2_PLAN.md`)
This is the only track that touches the novelty critique. Target ICLR/NeurIPS **D&B**
or single-cell methods — NOT a testing workshop.
1. **Close the coverage gap (the #1 novelty lever).** Implement
   `WeightedPerturbationConformal` (Tibshirani 2019 likelihood-ratio weighting) and/or
   Mondrian-by-cell-line in `src/confpert/conformal.py` + a K562↔RPE1 domain classifier;
   show cross-cell-line coverage recovers. Reuses cached scores → near-zero Modal.
   **Run `/code-review ultra` on this before trusting any coverage number.**
2. **Pre-registration as a first-class methodological section** (both rating-4 + Unm6
   point here — most defensible novel core at D&B).
3. **Breadth (usGP):** Phase 2's 14×10, more perturbation types.
4. **STATE-600M on Tahoe (Unm6):** re-attempt the Modal-heartbeat-failed download via
   CPU pre-stage; completes the most H1-informative cell.
5. **Held-out / wet-lab-anchored K3** to replace the predictor-derived substrate.
6. **Main/appendix rebalance (Unm6):** at D&B (no page limit) pull the 4 conformal head
   levels + their relation to CD-split/CQR/OT-CP/Mondrian into main text.
7. **Optional comparators** in related work: Wei 2025 (*Nat Methods*, 27×29), TISSUE
   (Sun 2023, calibration→FDR).

## 4. Honest novelty assessment (do not oversell)
ConfPert *measures* the coverage gap; it does not *close* it. Conformal heads + 6
discrepancies are off-the-shelf (paper itself says so). Novelty ~6–6.5/10 as-is; ceiling
~6.5–7/10 even with Track B (weighted CP = Tibshirani 2019, not new theory). It will NOT
become a testing-methods paper — the rating-4 reviewer is correct. Defensible cores:
(a) pre-registration protocol, (b) first coverage framework for this domain, (c) sVAE+
sparse-additive exact-nominal result. The post-hoc K562-vs-non-K562 finding must stay a
flagged secondary (it is not pre-registered).

## 5. Key files for ultraplan
| Path | Role |
|---|---|
| `paper_ai4science/main.tex` | THE reviewed paper (edit target) |
| `paper_ai4science/refs.bib` | 22 entries, all verified |
| `paper_ai4science/REVIEWER_RESPONSE_AUDIT.md` | full per-item audit + Tier-0 verdict + reviewer skeleton |
| `paper_ai4science/REVISION_PLAN_2026-05-31.md` | Track A / Track B plan keyed to verbatim reviews |
| `src/confpert/conformal.py` | 5 heads, NO weighted-CP (do not claim it has one) |
| `scripts/make_paper_figures.py` | Fig 3 generator |
| `scripts/make_fig_k3_robustness.py` | Fig 4 generator |
| `results/k2_preliminary.json` | live K2 (confirm canonical) |
| `results/k2_state_n8.json` | STALE — do not use |
| `baselines/results.json` | K1/K2 score rows; Fig 3 reproduces Table 2 to 6dp |
| `PHASE_2_PLAN.md` | Track B target (14×10, ICLR D&B) |

## 6. Verification gate (run after any edit)
```bash
cd confpert/paper_ai4science
tectonic -X compile main.tex --reruns 4                  # exit 0
pdftotext main.pdf - | grep -c "Title Suppressed"         # 0
pdftotext main.pdf - | grep -c "(?)"                      # 0
grep -ci undefined /tmp/cf.log                            # 0 (from compile log)
awk 'NR<213 && /scripts\/|\.py|modal\\_launch/' main.tex   # empty (body code-paths)
awk '{o+=gsub(/{/,x);c+=gsub(/}/,y)} END{print o-c}' main.tex   # 0
cd .. && python -m confpert.cli style-audit paper_ai4science/main.tex --strict  # 0/0
```
Plus human: `humanizer` + read-aloud (LLM-style), open each fig PDF (legibility).

## 7. Decisions needed from user (carry into the plan)
1. Run Track A now? (safe, $0)
2. Is `k2_preliminary.json` canonical? (blocks Fig 3 / Table 2)
3. Track B = fold into existing Phase 2, or separate?
4. Camera-ready blind or named?
5. Draft the reframed contribution paragraph, or author writes it?
6. Commit the uncommitted fixes (so remote tools see them)?
