# ConfPert — Revision Plan (from verbatim AI4Research reviews)

**Date:** 2026-05-31 · **Decision:** Accept (Poster), ICML 2026 AI4Research ·
**Reviews:** yhtC (3, clear-reject, integrity), workshop-fit reviewer (4, reject —
venue/novelty/LLM-style), Unm6 (7, accept, validity cons), usGP (7, accept, breadth).
**Reviewed file:** `paper_ai4science/main.tex` (verified). Prior pass:
`REVIEWER_RESPONSE_AUDIT.md`. State today: compiles exit 0, Title-Suppressed 0,
literal `(?)` 0, braces balanced, style-audit 0/0, 12pp. Body still has clean
library/API identifiers + 2 git hashes (`confpert`, `PerturbationConformal`,
`CDSplitConformal`, `CauchoisSubgroupConformal`, `RCPSCalibrator`, `SCORES`,
`confpert.cell_eval_plugin`, `ArcInstitute/cell-eval`, `c7046e4b`, `2db6b9c6`) —
defensible in a tool paper; yhtC's 4 NAMED offenders are all in the appendix.

---

## 0. The strategic call that organizes everything

The paper is **already accepted as a poster**, and the rating-4 reviewer said plainly:
*"for a single-cell methods venue or a benchmarks track, this would be a stronger fit;
... I cannot see how the paper can be developed into a testing-focused methodological
paper."* Two consequences:

1. **Testing-workshop framing is a confirmed dead end.** Do not try to make this a
   hypothesis-testing methods paper — the reviewer is right, and chasing it wastes effort.
2. **Two cleanly separable goals — do not conflate:**

| | **Track A — camera-ready THIS poster** | **Track B — the archival paper** |
|---|---|---|
| Goal | ship a clean published AI4Research poster | survive the novelty critique |
| Venue | AI4Research (already in) | ICLR/NeurIPS **D&B** or single-cell methods (= existing Phase 2) |
| Cost | low, ~half a day, $0 | weeks, Modal $ |
| New science | none | yes (shift-valid CP, broader benchmark) |
| Addresses | yhtC (integrity), Unm6 (clarity/validity), LLM-style | workshop-fit (novelty), usGP (breadth) |

The reviews are free, high-quality guidance for **Phase 2** (already specified in
`PHASE_2_PLAN.md`: 14 predictors × 10 datasets, ICLR D&B). Track B = fold reviews into
Phase 2; don't start a third thing.

**Honest novelty ceiling (research-grade standards):** even with Track B done, this is a
*benchmark + pre-registration protocol + first coverage framework for the domain* —
ceiling ~6.5–7/10 at a D&B venue. Weighted CP uses Tibshirani 2019 (off-the-shelf), so
it adds no new conformal theory and will NOT make this a testing-methods paper. State
this plainly; do not oversell.

---

## 1. Audit workflow (the verification gate)

A reviewer point is "closed" only when it passes its objective check. Run after every edit.

```bash
cd confpert/paper_ai4science
tectonic -X compile main.tex --reruns 4                          # exit 0
pdftotext main.pdf - | grep -c "Title Suppressed"                 # 0  (yhtC header)
pdftotext main.pdf - | grep -c "(?)"                              # 0  (yhtC broken cites)
grep -c "undefined" /tmp/cf.log                                   # 0  (cites resolve)
grep -nE 'tahoe\\_subset\\_build|score\\_variance\\_ratio\\_dev|categorical\\_attributes\\_keys|modal\\_launch' main.tex | awk -F: '$1<213'  # empty (yhtC's 4 NAMED tokens out of body; library/API names + git hashes may remain — defensible)
awk '{o+=gsub(/{/,x);c+=gsub(/}/,y)} END{print o-c}' main.tex      # 0  (brace balance)
cd .. && python -m confpert.cli style-audit paper_ai4science/main.tex --strict   # 0/0
```
Two human-judgment checks the script cannot do:
- **LLM-style** — `humanizer` skill + read-aloud (style-audit 0/0 is necessary, not
  sufficient; the rating-4 reviewer judged this by eye).
- **Figure legibility** — open each regenerated PDF, confirm at print size.

---

## 2. Per-reviewer status (verbatim-keyed, verified today)

### yhtC (clear-reject — integrity & presentation): ALL closed.
| Item (verbatim) | Status | Evidence |
|---|---|---|
| Fabricated Csendes title | ✅ FIXED | real title in refs.bib; 22 entries re-verified |
| 3 broken `(?)` (§2.1 energy, §3.1 Tahoe, §8 distr-free) | ✅ FIXED | 0 `(?)`, 0 undefined |
| Fig 3 labels stacked 3/4 panels | ✅ FIXED | regenerated, de-overlapped, visually confirmed |
| Fig 4(b) overlapping x-ticks | ✅ FIXED | single tick set, visually confirmed |
| Body dense w/ code identifiers (§5/6/8): `tahoe_subset_build`, `score_variance_ratio_dev`, `modal_launch.py::…`, `categorical_attributes_keys=[pert_col]` | ✅ FIXED | all 4 NAMED tokens out of body. The last body offender (§Reproducibility L207 — `run_k1_cross_cell_line.py`, `modal_launch.py::tahoe_subset_build`, `::state_calibrate`) was rewritten to plain prose pointing to Appendix~\ref{app:compute} (full recipe already duplicated there). Body now has zero command strings / script paths; only clean library/API names + 2 git hashes remain (defensible). Optional polish: move `confpert.cell_eval_plugin` / `ArcInstitute/cell-eval` to appendix too |
| "Title Suppressed" header | ✅ FIXED | re-assert `\@icmltitlerunning`; PDF count 0 |

### workshop-fit reviewer (rating 4 — venue / novelty / LLM-style / style file):
| Item (verbatim) | Status | Track |
|---|---|---|
| "not quite on-topic for the workshop" | ⚠️ ACK, unfixable for this venue | B: target D&B / single-cell |
| "style file … modified layout → desk-reject at ICML proper" | ✅ LIKELY FIXED (= the Title-Suppressed symptom; sty is the official 2025-10-29 ICML file, unmodified) | A: 10-min template check (§3.5) |
| "reads very much in an LLM style" | ⚠️ PARTIAL | A: `humanizer` + read-aloud |
| "methodological novelty unclear … off-the-shelf" | ❌ NOT fixed by edits | B: real work (§4); reframe = positioning only |
| reviewer's own steer: "pre-registration as transferable methodology for honest ML benchmarking under multiple testing" | 🎯 ADOPT | A: light reframe; B: first-class section |

### Unm6 (accept — validity & clarity):
| Item (verbatim) | Status |
|---|---|
| Cross-cell-line guarantee not flagged; Energy→0 should be foregrounded | ✅ FIXED (caption caveat + body §8) |
| STATE-600M missing from Tahoe | ⚠️ PARTIAL — stated in body; **re-run = Track B** |
| No glossary | ✅ FIXED (appendix glossary) |
| Main 5pp / appendix 12pp inversion | ⚠️ FLAGGED — blocked by 4pp limit; **fix in Track B** (D&B has room) |
| (pro) sVAE+ exact-nominal theoretically grounded | keep as lead empirical result |

### usGP (accept — breadth):
| Item | Status |
|---|---|
| "broader range of perturbation types and baselines" | ⚠️ future work → **Track B** = Phase 2's 14×10 |

### Correction to prior audit — Wei/TISSUE/VCC citations
Wei 2025 / TISSUE 2023 were **my** scoop-defense suggestions; **no verbatim reviewer
asked for them.** Reclassify as *optional* Track-B related-work strengthening, **not**
camera-ready blockers. The one cheap honesty item: the VCC quote in §1 ("almost all
models performed worse than baseline on MAE") is currently unattributed — add the
Roohani 2025 *Cell* cite in Track A.

---

## 3. Track A — finish the camera-ready poster (do first; ~half a day, $0)

Mostly done. Remaining, in order:

1. **Confirm `results/k2_preliminary.json` is the canonical K2 source.** Fig 3 + Table 2
   depend on it; `k2_state_n8.json` is stale. — **needs user.**
2. **LLM-style pass.** Run `humanizer` on body prose, then read-aloud edit; re-run
   style-audit `--strict` (keep 0/0), diff for any changed number. Targets the rating-4
   reviewer's most actionable presentation complaint. File: `main.tex`.
3. **Contribution-paragraph reframe (L68).** Adopt the reviewer's own suggestion: lead
   with the pre-registration protocol (two reviewers value it), sVAE+ as strongest
   empirical result, K562 as explicitly post-hoc secondary; name the conformal heads as
   deliberately off-the-shelf. Author's words preferred; I can draft.
4. **Attribute the VCC quote (§1)** — add Roohani 2025 *Cell* citation (1 bib entry + 1
   `\citep`).
5. **Verify workshop template (10 min).** Confirm AI4Research camera-ready uses the
   standard ICML style with no separate workshop `.sty`/banner. Low probability of an
   issue (header comment already says ICML 2026); just close it.
6. **Re-run the §1 gate** — all green; visually confirm both figures.
7. **(Optional) author block** — if camera-ready is non-blind, `\usepackage{icml2026}` →
   `[accepted]` + real author block.
8. **(Optional) finish yhtC code-id cleanup** — the §Reproducibility body command
   strings (the last named offenders) are already moved to the appendix this session;
   remaining is pure polish: the 2 borderline body tokens (`confpert.cell_eval_plugin`,
   `ArcInstitute/cell-eval`) could also go to the appendix. Not a blocker.

NOT in Track A: Wei/TISSUE cites, weighted-CP code, STATE re-run, main/appendix rebalance.

---

## 4. Track B — archival paper (fold reviews into Phase 2; where novelty is actually addressed)

Maps 1:1 onto `PHASE_2_PLAN.md`.

1. **Close the coverage gap, don't just measure it (the #1 novelty lever).** Implement
   `WeightedPerturbationConformal` (Tibshirani 2019 likelihood-ratio weighting) and/or
   Mondrian-by-cell-line in `src/confpert/conformal.py` (+ a K562↔RPE1 domain
   classifier), and **show cross-cell-line coverage recovers.** Converts the honest "it
   collapses, fix is future work" into "we restore finite-sample validity under shift."
   Reuses cached scores → near-zero Modal. **Run `/code-review ultra` on this before
   trusting any coverage number** (coverage math is subtly easy to get wrong).
2. **Pre-registration as a first-class contribution.** Both the rating-4 reviewer and
   Unm6 point here. Write the transferable-protocol section (git-stamp + categorical
   disposition + multiple-testing correction as a reusable recipe). Most defensible novel
   core at a D&B venue.
3. **Breadth (usGP).** Phase 2's 14 predictors × 10 datasets; more perturbation types.
4. **STATE-600M on Tahoe (Unm6).** Re-attempt the documented Modal-heartbeat failure via
   the CPU pre-stage workaround; if it lands, completes the most H1-informative cell.
5. **Held-out / wet-lab-anchored K3** to replace the predictor-derived substrate (turns
   "more pathways pass FDR" into "more *true* pathways").
6. **Main/appendix rebalance (Unm6).** At D&B (unlimited length) pull the four conformal
   head levels + their relation to CD-split/CQR/OT-CP/Mondrian into the main text.
7. **Venue (rating-4 reviewer's explicit steer).** ICLR D&B 2027 / single-cell methods —
   NOT a testing workshop. Add optional comparators (Wei 2025, TISSUE 2023) to related
   work here.

---

## 5. Decisions needed from user

1. **Run Track A now?** (Safe, $0, ~half a day. I can do 2–6 autonomously; #1 needs you.)
2. **Is `k2_preliminary.json` the canonical K2 result?** (Track A blocker for figure/table.)
3. **Track B = fold into existing Phase 2, or treat separately?** (Recommend: fold.)
4. **Camera-ready blind or named?** (Switches the style option + author block.)
5. **Draft the reframed contribution paragraph (A3) for your edit, or leave fully to you?**

## 6. Risks / honest caveats

- **The novelty critique is real; only Track B touches it.** Track A makes the poster
  clean but does not change the paper's standing as "a benchmark, not a method." Don't
  expect Track A to move a future reviewer's novelty score.
- **`humanizer` can drift numbers/claims if run carelessly** — restrict to prose, re-run
  the full gate after, diff every number.
- **Weighted-CP (B1)** is highest-value + highest-care — mandatory `/code-review ultra` +
  self-check agents before trusting coverage numbers.
- **`/ultraplan` failed** because cwd `~/FLM4S` isn't a git repo (git root is `confpert/`).
  Re-run remote planning from inside `confpert/` for the cloud version.
- **Prior-audit drift caught:** Wei/TISSUE/VCC were over-stated as needed; reclassified
  here as optional. No verbatim reviewer requested them.
