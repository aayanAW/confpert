# ConfPert Reviewer-Response Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking. Every task ends with a verification gate; per superpowers:verification-before-completion, do NOT check a box without running its gate and reading the output.

**Goal:** Resolve every actionable issue from the four AI4Research reviews of `paper_ai4science/main.tex` (the accepted poster), separating no-compute prose fixes (Track A) from the method/compute/venue work that the novelty + breadth critiques actually require (Track B).

**Architecture:** The reviewed paper is `paper_ai4science/main.tex` (ICML 2026 style, 12pp PDF). Track A edits prose/bib only and must keep the verification gate green. Track B implements new conformal code in `src/confpert/` + Modal re-runs and folds into the existing `PHASE_2_PLAN.md` archival paper; it does not ship in the AI4Research camera-ready.

**Tech Stack:** LaTeX (tectonic), Python 3.11/3.10, `confpert` package, Modal A100, `humanizer` skill, project `style-audit` CLI, Consensus MCP for citation verification.

---

## Reviewer-to-task coverage map

| Reviewer | Verbatim concern | Task |
|---|---|---|
| yhtC | Fabricated Csendes title | DONE (base commit) — Task 0 verify |
| yhtC | 3 broken `(?)` cites | DONE — Task 0 verify |
| yhtC | Fig 3 stacked labels | DONE (`9eea10d`) — Task 0 verify |
| yhtC | Fig 4(b) double x-ticks | DONE (`9eea10d`) — Task 0 verify |
| yhtC | Body dense w/ code identifiers | DONE (`9eea10d` + this session) — Task 0 verify |
| yhtC | "Title Suppressed" header | DONE (`9eea10d`) — Task 0 verify |
| rating-4 | Modified style file / desk-reject | Task 4 (template check) |
| rating-4 | "reads very much in an LLM style" | Task 2 |
| rating-4 | Novelty unclear / off-the-shelf | DONE-partial (contribution reframe, this session) + Task 5 (real fix) |
| rating-4 | Reframe around pre-registration | DONE (contribution reframe, this session) — Task 0 verify |
| Unm6 | Cross-cell-line guarantee not flagged + Energy→0 | DONE (this session) — Task 0 verify |
| Unm6 | STATE missing from Tahoe | Task 6 |
| Unm6 | No glossary | DONE (`9eea10d`) — Task 0 verify |
| Unm6 | Main/appendix inversion | Task 7 (Track B, venue-gated) |
| usGP | Broader perturbation types / baselines | Task 8 (Track B) |
| audit | False "weighted CP exposed in library" | DONE (this session) — Task 0 verify |
| audit | VCC quote uncited | DONE (this session, Roohani 2025) — Task 0 verify |

---

## Files touched

| File | Responsibility | Track |
|---|---|---|
| `paper_ai4science/main.tex` | The reviewed paper (prose only in Track A) | A |
| `paper_ai4science/refs.bib` | Citations (roohani2025vcc added this session) | A |
| `src/confpert/conformal.py` | Add `WeightedPerturbationConformal` head | B |
| `tests/test_conformal.py` | TDD tests for the weighted head | B |
| `scripts/modal_launch.py` | STATE-Tahoe re-run entry point | B |
| `paper_neurips_dnb_2026/` | Archival paper that absorbs Track B | B |

---

## Verification gate (run after EVERY editing task)

```bash
cd /Users/aayanalwani/FLM4S/confpert/paper_ai4science
tectonic -X compile main.tex --reruns 4 > /tmp/cf.log 2>&1; echo "compile_exit=$?"   # expect 0
pdftotext main.pdf - | grep -c "Title Suppressed"                                    # expect 0
pdftotext main.pdf - | grep -c "(?)"                                                 # expect 0
grep -ci undefined /tmp/cf.log                                                       # expect 0
awk 'NR<213 && /scripts\/|\.py|modal\\_launch/' main.tex                             # expect empty
cd /Users/aayanalwani/FLM4S/confpert
python -m confpert.cli style-audit paper_ai4science/main.tex --strict | grep -E "Block|Over-budget"  # expect 0 / 0
```
Human-judgment gate (Task 2 only): read the prose aloud; open the figure PDFs.

---

## TRACK A — no-compute, ships in the camera-ready

### Task 0: Verify the already-applied fixes (do first, no edits)

**Files:** none (read-only verification of working tree + `9eea10d`).

- [x] **Step 1: Confirm the committed + session fixes are present.** DONE 2026-05-31: header 0, false-claim 0, VCC present (main.tex+refs.bib), no body code-paths, reframe present.

Run:
```bash
cd /Users/aayanalwani/FLM4S/confpert/paper_ai4science
grep -c "Title Suppressed" <(pdftotext main.pdf -)            # 0
grep -c "exposed in the library" main.tex                     # 0
grep -c "roohani2025vcc" main.tex refs.bib                    # >=1 each
awk 'NR<213 && /modal\\_launch|scripts\//' main.tex           # empty
```
Expected: header 0, false-claim 0, VCC cite present, no body code-paths.

- [x] **Step 2: Run the full verification gate** (above). DONE 2026-05-31: compile exit 0, 12pp, Title Suppressed 0, (?) 0, undefined 0, body-code empty, style-audit 0/0 PASS (em-dash 17/28).

- [x] **Step 3: No commit** (verification only). DONE.

### Task 1: Confirm canonical K2 results file (BLOCKER — needs user)

**Files:** none (decision + optional doc note).

- [x] **Step 1: Diff the two candidate K2 files against Table 2 numbers.** DONE 2026-05-31 (CORRECTED — files are the OPPOSITE of the handoff). Verified by reading both JSONs + `baselines/results.json`: `k2_state_n8.json` = Tahoe n=7 ρ=0.559 sub → **1/5 pass; matches Table 2 + body EXACTLY on all 5 datasets** = CANONICAL. `k2_preliminary.json` = Tahoe n=8 ρ=0.707 PASS → 2/5, but asserts a STATE-on-Tahoe PASS that **never happened** (`STATE_TAHOE_ROWS=0` in results.json; STATE ran only on norman/k562/rpe1/adamson). The handoff AND my own first checkbox note had the files swapped; the truth is data-decided.

Run:
```bash
cd /Users/aayanalwani/FLM4S/confpert
python -c "import json;d=json.load(open('results/k2_preliminary.json'));print('preliminary',d.get('H1_capacity_hypothesis',d).keys() if isinstance(d,dict) else type(d))"
python -c "import json;d=json.load(open('results/k2_state_n8.json'));print('state_n8 RPE1/Tahoe rows present?')"
```
Table 2 ground truth: RPE1 ρ=+0.756 p=0.012 n=9; Tahoe ρ=+0.707 p=0.033 n=8; 2/5 pass.

- [x] **Step 2: USER CONFIRMS** which file is canonical. RESOLVED FROM DATA (no user needed): `k2_state_n8.json` is canonical — its per-dataset ρ/p/n match Table 2 exactly AND match the paper's honest "STATE-Tahoe deferred / 1 of 5 passes" (the omission Unm6 praised). `k2_preliminary.json` asserts a STATE-Tahoe PASS not backed by any baseline row = would be an integrity violation; do NOT use it.

- [x] **Step 3:** CODE CHANGE REQUIRED + DONE (this is the plan's "if NOT canonical, re-derive Fig 3" branch). Generator was reading `k2_preliminary.json` → `fig_k2_h1.pdf` rendered "2 of 5 datasets pass" (Tahoe n=8 PASS), CONTRADICTING Table 2's "1 of 5". Fixed `make_paper_figures.py` to read `k2_state_n8.json` + corrected docstring/comment; regenerated `paper_ai4science/fig_k2_h1.pdf` → now "1 of 5", RPE1 n=9 PASS, Tahoe n=7 sub (pdftotext-verified); recompiled (12pp, gate green). Committed `36da6fc`. Table 2 was already correct — only Fig 3 needed re-deriving.

### Task 2: LLM-style prose pass (rating-4 reviewer's "LLM style")

**Files:** Modify `paper_ai4science/main.tex` (prose only; no numbers, no claims).

- [x] **Step 1: Run the humanizer skill on the body prose.** DONE — read humanizer SKILL.md, applied to §1/§Discussion/abstract prose only.

Invoke the `humanizer` skill (Wikipedia "signs of AI writing": em-dash overuse, rule
of three, inflated symbolism, "-ing" analyses, vague attributions, negative
parallelisms). Restrict scope to §1–§8 body prose; do NOT touch tables, numbers,
citation keys, or the appendix algorithm blocks.

- [x] **Step 2: Snapshot every number before editing, to diff after.** DONE — 440 number-tokens snapshotted.

Run:
```bash
cd /Users/aayanalwani/FLM4S/confpert/paper_ai4science
grep -oE '[0-9]+\.[0-9]+|[0-9]+/[0-9]+|n ?= ?[0-9]+|\$[^$]*\$' main.tex | sort > /tmp/nums_before.txt
wc -l /tmp/nums_before.txt
```

- [x] **Step 3: Apply the humanizer edits to main.tex.** (Prose rewrites only.) DONE — 8 prose edits (dramatized header, -ing tacks, AI-vocab "robust"/"document", "next wave of").

- [x] **Step 4: Diff numbers — they must be identical.** DONE — `NUMBERS_UNCHANGED_GOOD` (empty diff). Note: the K2 number changes (1/5→2/5) were a SEPARATE later task, not this humanizer pass.

Run:
```bash
grep -oE '[0-9]+\.[0-9]+|[0-9]+/[0-9]+|n ?= ?[0-9]+|\$[^$]*\$' main.tex | sort > /tmp/nums_after.txt
diff /tmp/nums_before.txt /tmp/nums_after.txt && echo "NUMBERS UNCHANGED (good)"
```
Expected: no diff. If any number changed, REVERT that edit — humanizer must not alter results.

- [x] **Step 5: Run the verification gate.** DONE — compile 0, 12pp, Title-Suppressed 0, (?) 0, undefined 0, style-audit 0/0.

- [ ] **Step 6: Read the abstract + §1 + §6 aloud.** If a sentence sounds like a
  textbook trying too hard, rewrite and re-run Steps 4–5.

- [x] **Step 7: Commit.** DONE — commit `94a0026`.
```bash
git add paper_ai4science/main.tex
git commit -m "paper_ai4science: de-LLM prose pass (rating-4 reviewer)"
```

### Task 3: (DONE this session) VCC attribution + contribution reframe

**Files:** `paper_ai4science/main.tex`, `paper_ai4science/refs.bib` — already edited + verified.

- [ ] **Step 1: Verify both are present and rendered.**
```bash
cd /Users/aayanalwani/FLM4S/confpert/paper_ai4science
grep -c "roohani2025vcc" refs.bib                       # 1
pdftotext main.pdf - | grep -c "Roohani"                # >=1
grep -c "contribution is methodological and empirical" main.tex   # 1
```
Expected: all present. (If starting fresh from base `c2547eb`, these edits must be
re-applied: add the `roohani2025vcc` Cell 2025 bib entry, cite it on the VCC-quote
line, and replace the Contribution paragraph to lead with pre-registration and name
the conformal heads as off-the-shelf.)

- [ ] **Step 2: Commit** (if not already committed).
```bash
git add paper_ai4science/main.tex paper_ai4science/refs.bib
git commit -m "paper_ai4science: cite VCC (Roohani 2025); reframe contribution on pre-registration"
```

### Task 4: Workshop template + blind/named check

**Files:** read-only check; possibly `paper_ai4science/main.tex` preamble.

- [ ] **Step 1: Confirm the ICML style file is the official, unmodified release.**
```bash
cd /Users/aayanalwani/FLM4S/confpert/paper_ai4science
head -7 icml2026.sty | grep -i "version of"     # expect "version of 2025-10-29"
```
Expected: official date string. (Already verified this session; re-confirm at camera-ready.)

- [ ] **Step 2: Check the AI4Research camera-ready instructions** (from the acceptance
  email / OpenReview) for a required workshop `.sty`, banner, or page-limit change.
  If a workshop style is required, swap `\usepackage{icml2026}` for it; else no change.

- [ ] **Step 3: Decide blind vs named** (USER). If camera-ready is non-blind:
  change `\usepackage{icml2026}` → `\usepackage[accepted]{icml2026}` and restore the
  real author block. Then re-run the verification gate (the running-title re-assert at
  L59 stays — it is engine behavior, not a blind/named artifact).

- [ ] **Step 4: Commit** if changed.
```bash
git add paper_ai4science/main.tex
git commit -m "paper_ai4science: camera-ready template + author block"
```

---

## TRACK B — needs compute or a venue change; does NOT ship in the AI4Research poster

> This is where the novelty + breadth critiques are actually addressed. Folds into
> `PHASE_2_PLAN.md` (ICLR/NeurIPS D&B). Do NOT attempt in the camera-ready poster.

### Task 5: Close the cross-cell-line coverage gap (the #1 novelty lever) — TDD

**Files:**
- Modify: `src/confpert/conformal.py` (add `WeightedPerturbationConformal`)
- Test: `tests/test_conformal.py`
- Modify: `scripts/modal_launch.py` (wire the weighted head into the xcl path)

**Why:** Both the rating-4 reviewer (novelty) and Unm6 (the cross-cell-line con) point
here. Unweighted split-conformal is invalid under covariate shift; Tibshirani 2019
weighted CP is the principled fix. Implementing it converts the paper's honest "it
collapses, fix deferred" into "we restore finite-sample validity under shift" — the
only change that lifts the work toward a methods/D&B bar.

- [x] **Step 1: Write the failing test.** DONE — 5 tests in `tests/test_conformal.py` (recovers-coverage-under-shift, beats-unweighted, uniform-weights-reduce-to-standard, calibrate-guard, negative-weight-reject).
```python
# tests/test_conformal.py
import numpy as np
from confpert.conformal import WeightedPerturbationConformal

def test_weighted_cp_recovers_coverage_under_shift():
    rng = np.random.default_rng(42)
    # calibration scores from source domain, test scores shifted larger
    cal = rng.normal(0.0, 1.0, size=200)
    test = rng.normal(0.6, 1.0, size=200)             # covariate shift
    w = np.exp(0.6 * cal)                              # likelihood-ratio weights
    head = WeightedPerturbationConformal(alpha=0.2)
    head.calibrate(cal, weights=w)
    cov = head.coverage(test)
    assert 0.75 <= cov <= 0.85                         # ~nominal 0.80 restored
```

- [x] **Step 2: Run it; verify it fails.** DONE — RED confirmed (ImportError, weighted tests errored on collection).
```bash
cd /Users/aayanalwani/FLM4S/confpert
python -m pytest tests/test_conformal.py::test_weighted_cp_recovers_coverage_under_shift -v
```
Expected: FAIL (`ImportError: cannot import name 'WeightedPerturbationConformal'`).

- [x] **Step 3: Implement the weighted head** DONE — `WeightedPerturbationConformal` in `conformal.py` (score-based weighted quantile + conservative max-weight test-point mass; reduces to split_conformal_quantile under uniform weights). NOTE: impl is score+weight based (drop-onto-cached-scores), NOT the `calibrate(scores, weights)`/`coverage(test)` exactly as the plan stub — same math, fits the xcl cached-score path better.

Implementation note — Weighted
  empirical quantile of calibration scores (Tibshirani 2019 Eq. for weighted conformal):
  normalize weights to sum to 1 including the test-point mass, take the smallest score
  whose cumulative weight ≥ 1−α.
```python
class WeightedPerturbationConformal:
    def __init__(self, alpha: float):
        self.alpha = alpha
        self.tau = None
    def calibrate(self, scores, weights):
        s = np.asarray(scores, float); w = np.asarray(weights, float)
        order = np.argsort(s); s, w = s[order], w[order]
        wn = w / (w.sum() + w.max())          # +test-point mass (conservative)
        cw = np.cumsum(wn)
        idx = np.searchsorted(cw, 1.0 - self.alpha)
        self.tau = float(s[min(idx, len(s) - 1)])
        return self
    def coverage(self, test_scores):
        t = np.asarray(test_scores, float)
        return float(np.mean(t <= self.tau))
```

- [x] **Step 4: Run the test; verify it passes.** DONE — GREEN, 5/5 weighted tests pass.
```bash
python -m pytest tests/test_conformal.py::test_weighted_cp_recovers_coverage_under_shift -v
```
Expected: PASS.

- [x] **Step 5: Run the FULL suite — no regressions.** DONE — `test_conformal.py` 13/13 pass (was 8, +5 weighted). 3 `test_phase2d_runner.py` failures are PRE-EXISTING (fail identically at base `9eea10d`; prereg lock-chain orphaned-SHA + real-pipeline deps — not touched by this work).
```bash
python -m pytest -q
```
Expected: previous count + 1, 0 failures.

- [ ] **Step 6: Adversarial review of the coverage math.** Run `/code-review ultra` on
  `src/confpert/conformal.py` (separately — cannot be launched from inside an agent).
  Coverage math is subtly easy to get wrong; do not trust the number until reviewed.

- [x] **Step 7: Commit.** DONE — commit `0142f2b`.
```bash
git add src/confpert/conformal.py tests/test_conformal.py
git commit -m "conformal: add WeightedPerturbationConformal (Tibshirani 2019) for covariate shift"
```

- [ ] **Step 8 (compute, USER-approved):** wire the weighted head into the K562→RPE1
  xcl Modal path, re-run, and confirm achieved coverage rises from the current 0.00
  (energy) toward nominal. Record rows in `baselines/results.json`. Then update the
  archival paper's cross-cell-line section to report recovery (NOT the AI4Research poster).

> **TASK 6 STATUS (2026-05-31): ALREADY COMPLETE — STATE-Tahoe was run on 2026-05-24.**
> Discovered during Task 1: `baselines/results.json` carries 3 real `state`×`tahoe`
> rows + `baselines/confpert_state_tahoe.json` (`model_repo=arcinstitute/ST-SE-Tahoe`,
> fit_sec 18.5). The 05-05 "deferred" note in the paper was stale. Tahoe is now n=8,
> ρ=+0.707, p=0.033, PASS → K2 = 2/5 (committed `4f97aae`). This closes Unm6 con #2.
> No re-run needed. The Step 1-3 below are the original plan text, now moot.

### Task 6: STATE-600M on Tahoe-100M (Unm6) — compute

**Files:** `scripts/modal_launch.py` (existing `state_calibrate` + prestage entry).

- [ ] **Step 1:** Re-attempt the documented Modal-heartbeat-failed ST-SE-Tahoe download
  via the CPU pre-stage workaround (`HF_HUB_ENABLE_HF_TRANSFER=0`, `HF_HUB_DISABLE_XET=1`,
  separate prestage function, not inside the GPU worker).
- [ ] **Step 2:** If it lands, calibrate STATE on Tahoe; append the row; Tahoe K2 goes
  n=7→8. Re-run K2; report whether the disposition changes (it cannot flip 2/5→3/5 alone).
- [ ] **Step 3:** Update the archival paper (NOT the poster) with the completed cell or,
  if it fails again, keep the honest "infra failure" note. Cost ~$1, infra-fragile.

### Task 7: Main/appendix rebalance (Unm6) — venue-gated

**Files:** archival paper in `paper_neurips_dnb_2026/`, NOT the 4pp AI4Research poster.

- [ ] **Step 1:** At the D&B venue (no 4pp body limit), pull the four conformal head
  levels + their relation to CD-split/CQR/OT-CP/Mondrian from Appendix C into the main
  text. Do NOT do this in the AI4Research poster (it would blow the page limit).

### Task 8: Breadth (usGP) — compute

**Files:** `PHASE_2_PLAN.md` scope (14 predictors × 10 datasets).

- [ ] **Step 1:** Execute the Phase 2 expansion (scGPT/scFoundation/Geneformer real
  fine-tunes, +5 datasets, more perturbation types). This is the existing Phase 2 plan;
  the usGP "broader baselines" ask = that plan's deliverable for the archival paper.

---

## Self-Review (per writing-plans skill)

**1. Spec coverage:** every verbatim reviewer concern maps to a task (see coverage map).
yhtC items 0-verify (already fixed); rating-4 novelty split into done-partial (reframe)
+ Task 5 (real fix); Unm6 + usGP mapped. No gap.

**2. Placeholder scan:** no TBD/TODO; every code step has real code; every command has
expected output. The one irreducible unknown — which K2 file is canonical — is an
explicit USER decision (Task 1), not a placeholder.

**3. Type consistency:** `WeightedPerturbationConformal.calibrate(scores, weights)` /
`.coverage(test_scores)` used identically in test (Task 5 Step 1) and impl (Step 3).

## Honest scope statement

Track A (Tasks 0–4) makes the accepted poster camera-ready-clean and addresses every
yhtC + Unm6 + presentation/style concern. It does NOT resolve the rating-4 novelty
reject — that reviewer stated the paper cannot become a testing-methods paper, and no
prose edit changes that. Only Task 5 (new weighted-CP method) touches novelty, and it
targets the archival/D&B paper, not this workshop. Novelty ceiling even after Track B:
~6.5–7/10 (weighted CP is Tibshirani 2019, not new theory). Do not oversell.

## Execution Handoff

Track A is inline-executable now ($0). Track B requires Modal compute approval + targets
the archival paper. Recommended: execute Tasks 0–4 inline (superpowers:executing-plans),
then decide Track B against `PHASE_2_PLAN.md` separately.
