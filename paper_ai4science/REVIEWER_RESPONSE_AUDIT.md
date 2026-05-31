# ConfPert — Reviewer-Response Audit

**Paper audited:** `paper_ai4science/main.tex` (the reviewed file — confirmed the
correct target because it is the only manuscript whose section / table / figure
numbering matches the reviewer comments: §8 = Limitations, Table 3 = cross-cell-line,
Fig 4(b) = K3 ratio panel, and the ICML-only "Title Suppressed" header behaviour).
The 2026-05-24 reformats (`main_short.tex`, `main_paper.tex`, `supplementary.tex`)
and the NeurIPS-style `paper/main.tex` were excluded.

**Date:** 2026-05-30. **Method:** orchestrator + 4 parallel sub-agents (significance,
citation-integrity, scientific-validity, figures), Consensus-verified literature,
single-writer edit pass on `main.tex`.

## Orient — layout before/after

| | Before | After |
|---|---|---|
| Main body | ~5 pp (lines 1–207) | ~5 pp + glossary table |
| Appendix | ~7 pp (lines 208–448, 13 sections) | +1 section (Glossary) |
| Compiled PDF | 12 pp | 12 pp |
| Build | `tectonic -X compile main.tex --reruns 4` | unchanged |
| Style | `\usepackage{icml2026}` (blind) | unchanged |

**Compile-gate status (after fixes):** compile exit 0 · `Title Suppressed` in PDF = 0
(was on every non-title page) · undefined citations = 0 · literal `(?)` = 0 ·
style-audit `--strict` = 0 block / 0 over-budget (was 2 over-budget) · 5 cosmetic
`Overfull \hbox` from one long `\texttt{}` Modal-command line in the appendix
(pre-existing, not visible as overflow).

---

# TIER 0 — Significance & Acceptance Verdict (lead)

**One-line verdict.** ConfPert is a rigorously executed, honestly reported,
correctly motivated benchmark-and-tool occupying a genuinely empty niche
(conformal × perturbation × distributional coverage). Independent novelty
≈ **6–6.5/10** — a domain port of standard conformal machinery to a real,
well-identified gap, lifted above a pure benchmark by (i) the pre-registration
discipline and (ii) one transferable architecture-level finding (sVAE+). It is a
**clear workshop accept** once three citation gaps are closed, and a **main-track
reject until it adds a method that closes — not just measures — the
shift-coverage gap.**

### 1. Does the distributional-coverage gap matter, and to whom?
The *premise* is airtight and is the paper's strongest asset: mean-prediction
saturation is community consensus (Ahlmann-Eltze 2025 *Nat Methods*; Csendes 2025;
Kedzierska 2025; Kernfeld 2025), the data substrate is genuinely distributional,
and no published method puts a finite-sample coverage guarantee on the
predicted-vs-observed population discrepancy. That cell of the design space is
empty. **However**, the paper does not connect the gap to a concrete decision that
flips. For leaderboard consumers it produces a better-behaved scalar, not a
demonstrably different/better model ordering than the E-distance/MMD already in
Cell-Eval. The one place it bites is **K3 drug-screen triage** (calibrated
signatures recover 445 vs 363 BH-FDR pathways) — but that substrate is
predictor-derived, not a held-out wet-lab readout, so "more pathways pass FDR"
is not yet shown to be "more *true* pathways." Honest answer: the guarantee is
valuable in principle to the VCC/Arc community; only K3 demonstrates a
workflow-relevant payoff, and K3 is under-validated.

### 2. Does ConfPert close the gap or just measure it?
**It measures.** The reject-leaning reviewer is essentially correct. The method
section itself states all four conformal heads are off-the-shelf "with verbatim
citations"; the six discrepancies are standard. The technical core is a wiring
harness: standard discrepancies × standard conformal procedures × 8 predictors ×
5 datasets — a benchmark + library, not a new method. The headline shift-robustness
problem (cross-cell-line coverage collapse) is **named, not solved**: the stated
remedy (weighted CP) was not implemented (see Item 10c). The one genuinely
mechanistic, non-off-the-shelf observation is the sVAE+ sparse-additive
exact-nominal result — but that is an observation about an existing model surfaced
by the benchmark, not a contribution ConfPert makes.

### 3. Single strongest defensible novel core
Lead with the **sVAE+ / output-head design rule** as the scientific result
("architecture, not capacity, determines distributional calibration; a
marginal-Gaussian head provably cannot recover bimodality, 1/3 < 5/9"), framed by
the **pre-registration discipline** as the credibility spine. Demote the
**K562-vs-non-K562 covariate** to an explicitly post-hoc, replication-worthy
teaser — it is *not* pre-registered (the paper says so), rests on n = 5 with a 2/2
K562 split, and making it the headline would hand the hostile reviewer the exact
"cherry-picking after the fact" attack the pre-registration exists to prevent.
(The scientific-validity agent independently reached the same ranking but argued
the pre-registration *protocol* is the strongest claim; both agree the K562
finding must not lead. Recommended framing: protocol = spine, sVAE+ = lead
result, K562 = flagged secondary — see Item 12.)

### 4. Calibrated accept probability
- **Methods / AI4Science / hypothesis-testing workshop: ~70% (60–80%).**
  Rejection reasons a real reviewer cites: (1) "method is assembled off-the-shelf,
  low technical novelty" (forgiven at workshops if framed as benchmark+tool);
  (2) "the most decision-relevant guarantee — cross-cell-line coverage — collapses
  and the fix is deferred." Highest-leverage change: **add the three missing
  citations + a one-paragraph 'why not subsumed' defense** (Wei 2025, TISSUE 2023,
  the uncited VCC quote). Worth ~10 points on its own.
- **ICML/NeurIPS/ICLR main track: ~8% (5–12%) as currently scoped.** Stated
  plainly: main-track odds are low and stay low without new methodology. A
  benchmark whose core is off-the-shelf, whose pre-registered primary hypothesis
  (K2) *fails*, and whose key guarantee collapses with the fix deferred does not
  clear the bar in a subfield where Wei 2025 (*Nat Methods*) already offers a
  larger generalizable benchmark. Only realistic main-track home is NeurIPS D&B,
  and even there the combination caps it below the line. Highest-leverage change:
  **implement and validate weighted/covariate-shift CP that actually restores
  cross-cell-line coverage** — the only change that moves this from workshop-tier
  to main-track-plausible.

---

# TIER 1 — Integrity & desk-reject blockers

### 1.1 Fabricated Csendes title — **ALREADY FIXED (verified)**
- **Status:** resolved before this session; re-verified against primary source.
- **Location:** `refs.bib:11`.
- **Finding:** the fabricated "Train-Mean as a control…" title is gone; the bib
  now carries the real title *"Benchmarking foundation cell models for
  post-perturbation RNA-seq prediction,"* BMC Genomics 26(1):393, 2025,
  doi 10.1186/s12864-025-11600-2 (Csendes, Sanz, Szalay, Szalai). Confirmed via
  Consensus + Crossref + local lit_notes (PMID 40269681). The stale title survives
  only in the unused legacy `paper/refs.bib`.
- **Full bib sweep:** all 22 entries verified against primary sources; 0 mismatches,
  0 unverifiable. Special-scrutiny entries (otcp2025, adduri2025state, tahoe2025)
  all match titles/authors/IDs exactly. Cross-ref caveat: the *musiml*
  `CITATION_VERIFICATION.md` mis-attributes Tahoe-100M to "Gandhi et al."; the
  ai4science bib is correct (Zhang et al.). No action on the ai4science paper.

### 1.2 Three broken `(?)` citation markers — **ALREADY RESOLVED (verified)**
- **Status:** all three map correctly; 0 literal `(?)` and 0 undefined cites in the
  compiled PDF.
- **Locations / mappings:** per-gene energy → `szekely2013` (§2.1);
  Tahoe-100M → `tahoe2025` (§3.1); distribution-free conditional impossibility →
  `vovk2012` (§8 Limitations). All verified against primary sources.

### 1.3 Template overflow / "Title Suppressed Due to Excessive Size" — **FIXED**
- **Status:** fixed + visually verified (page-2 header now reads the real running
  title; `pdftotext` count of "Title Suppressed" = 0).
- **Location:** `main.tex:29` (running title) + `main.tex:55` (override).
- **Root cause (diagnosed, not guessed):** `icml2026.sty` lines 369–397 measure the
  running title in a `\vbox` and substitute the "Title Suppressed" fallback when
  `\ht\titrun > 6.25pt`. Under this engine the check misfires for *any* normal title
  (an 8-character test title also triggered it), so shortening alone does not fix it
  — confirmed by bisection (8 / 28 / 42 chars all triggered).
- **Fix (diff, no style-file edit):** restore a full, informative running title and
  re-assert it after `\icmltitle` executes:
  ```latex
  \icmltitlerunning{ConfPert: Conformal Coverage for Perturbation Predictors}
  ...
  \printAffiliationsAndNotice{}
  \makeatletter\gdef\@icmltitlerunning{ConfPert: Conformal Coverage for Perturbation Predictors}\makeatother
  ```
- **Style-file integrity:** `icml2026.sty` is the official 2025-10-29 ICML release
  (header verified); it was **not** modified. The reviewer's "modified layout"
  concern is addressed by leaving the sty untouched and fixing the symptom in the
  document preamble.

---

# TIER 2 — Presentation

### 4. Fig 3 (`fig_k2_h1.pdf`) — label overlap **+ stale/contradictory data** — **FIXED**
- **Status:** fixed + visually verified. This was more than cosmetics: the old
  figure's suptitle said "0/4 datasets pass" and showed 4 panels at n=5–8, which
  **contradicted the paper's own Table 2** (n=9, 2/5 pass, +Tahoe, RPE1 ρ=+0.756).
  A figure contradicting the headline table is an integrity issue.
- **Location:** generator `scripts/make_paper_figures.py::fig_k2_h1_scatter`;
  output `paper_ai4science/fig_k2_h1.pdf`.
- **Fix:** generator now reads the live K2 summary, plots 5 panels incl. Tahoe + the
  STATE point (star), replaces the hardcoded "0/4" suptitle with "2 of 5 datasets
  pass; H1 overall fails (needs ≥3/4). Star = STATE," and removes label overlap
  (`adjustText` + short codes). Recomputed Spearman from `baselines/results.json`
  reproduces Table 2 to 6 decimals on all 5 datasets — **no points fabricated.**
- **⚠ AUTHOR-CONFIRM FLAG:** the live numbers live in `results/k2_preliminary.json`,
  *not* `results/k2_state_n8.json` (which is itself stale at Tahoe n=7). Confirm
  `k2_preliminary.json` is the intended frozen K2 result.

### 5. Fig 4(b) (`fig_k3_robustness.pdf`) — double x-axis tick labels — **FIXED**
- **Status:** fixed + visually verified (single clean `{25,50,100,200}` tick set;
  ratios 1.20/1.23/1.20/1.43× and parity line render correctly).
- **Location:** `scripts/make_fig_k3_robustness.py` panel-b block.
- **Fix:** plot panel (b) on categorical positions and drop `set_xscale("log")`,
  which was emitting an overlapping `LogFormatterSciNotation` minor-tick set on top
  of the categorical labels. Panel (a) untouched.

### 6. Code identifiers in body prose — **PARTIALLY FIXED (author-decision flag)**
- **Status:** the worst body offender (the STATE Limitations paragraph, §8) was
  rewritten to plain prose with the version pins / checkpoint-loading internals
  moved to Appendix `app:predictors`. Remaining body `\texttt{}` tokens are the
  legitimate library API surface (`confpert`, `PerturbationConformal`, etc.) and the
  Cell-Eval commit hash — defensible in a tool paper.
- **Location (fixed):** §8 Limitations, `main.tex:203`.
- **Flag:** if the venue wants *zero* code identifiers in body prose, the
  `confpert.cell_eval_plugin` / commit-hash mentions (§ Library) can also move to the
  appendix — left as an author call, low value.

### 7. "LLM-style" prose — **FIXED to the project's own gate**
- **Status:** rule-of-three list removed from the abstract ("We find three
  things…" → three plain sentences); `EXACTLY` de-shouted; "strictly better" →
  "lower deviation"; "naturally produces" pleasantry removed. The project
  `style-audit --strict` now passes **0 block / 0 over-budget** (was 2 over-budget:
  `inflated_phrase_strictly_better`, `pleasantry_in_text`). Em-dash 18 (budget ~28).
- **Note:** a senior-author read-aloud pass is still worth one human pass before
  camera-ready; the automated gate is necessary, not sufficient.

---

# TIER 3 — Structure & clarity

### 8. Main/appendix inversion — **PARTIALLY ADDRESSED (author-decision flag)**
- **Status:** the four conformal head levels already appear in the body (§2.2,
  `main.tex:76`), with full algorithms in Appendix C. The reviewer's ask — pull the
  load-bearing head/CD-split/CQR/OT-CP/Mondrian *relationship* into the body — is
  blocked by the workshop 4-page body limit (the appendix is unlimited). Moving the
  algorithm blocks up would push the body over length.
- **Recommendation (flag):** if retargeting to an unlimited-length venue (NeurIPS
  D&B), promote Appendix C's head taxonomy into §2. For the 4-page workshop, the
  current body §2.2 summary + appendix detail is the correct trade-off. Author call.

### 9. Glossary / notation table — **FIXED**
- **Status:** added. Consolidates K1/K2/K3, H1/H1b, the four head levels, the six
  discrepancies, α, and calibration deviation.
- **Location:** new Appendix A "Glossary and Notation" (`tab:glossary`), each term
  also defined on first use in the body.

---

# TIER 4 — Scientific validity

### 10. Cross-cell-line conformal guarantee does not hold
- **(a) Table caption not caveated — FIXED.** `tab:xcl` caption now states the
  numbers are *not* coverage-guaranteed (unweighted split-conformal is not
  exchangeability-valid under the K562→RPE1 shift; reported as a shift-fragility
  diagnostic). Location: appendix `tab:xcl` caption.
- **(b) Energy→0 undercoverage signature buried — FIXED.** A sentence now
  foregrounds the collapse in the body §8 Limitations: "empirically achieved
  coverage collapses to 0.00 on the energy discrepancy across all six predictors
  (Appendix), the textbook undercoverage signature."
- **(c) Weighted CP "exposed in the library" — FALSE CLAIM, FIXED (mandatory).**
  Independently verified: `src/confpert/conformal.py` contains exactly five heads
  (PerturbationConformal, CDSplitConformal, EnergyConformal,
  CauchoisSubgroupConformal, RCPSCalibrator); **zero** weighted/Tibshirani/
  likelihood-ratio/covariate-shift code anywhere in `src/confpert/`. The paper
  claimed weighted CP was "exposed in the library" in three places (§2.3, §8,
  appendix); all three were rewritten to "the principled fix … left as future work
  / not yet implemented." This was the single most important correctness fix.
  - **⚠ AUTHOR/COMPUTE FLAG (Path B, optional):** the cheapest *correct*
    implementation is a `WeightedPerturbationConformal` head (weighted empirical
    quantile) + a K562-vs-RPE1 domain classifier — reuses cached scores, near-zero
    marginal compute, but non-trivial new code. Implementing it (and showing
    coverage recovers) is the single change that lifts the paper toward main-track.
    Path A (prose correction) is done and is mandatory regardless.

### 11. STATE-600M missing from Tahoe-100M — **PROSE FIX APPLIED; re-run flagged**
- **Status:** the §8 Limitations paragraph was rewritten to state the omission
  plainly. The reviewer wants the *impact on H1* prominent in the body.
- **Impact assessment (for the response letter):** Tahoe already passes the ρ
  magnitude without STATE (ρ=+0.707) and now counts as one of the 2/5 passes.
  STATE-on-Tahoe **cannot by itself flip the overall 2/5 → 3/5 disposition** (a
  third dataset would also need to pass). It could strengthen Tahoe's significance,
  but could equally flatten the trend if STATE calibrates well on its data-richest
  regime. The honest framing: "would complete the capacity axis on the most
  informative dataset; outcome unknown; current result stands without it."
- **⚠ AUTHOR/COMPUTE FLAG:** the re-run is the documented Modal-heartbeat failure
  (~$1, infra-fragile, 5 prior retries failed). The prose is publishable whether or
  not the re-run is attempted. Recommend a one-sentence body addition at the K2
  result (optional, author call) making the impact explicit.

### 12. Novelty framing — **AUTHOR-DECISION (recommendation provided)**
- Both judgment agents agree: do **not** oversell the conformal machinery as novel
  (the paper already concedes it is off-the-shelf), and do **not** lead with the
  post-hoc K562 covariate. Recommended one-sentence contribution reframing: lead
  with the pre-registered honest-benchmarking protocol + the sVAE+ output-head
  design rule; name the conformal heads as deliberately off-the-shelf; carry the
  K562 result as an explicitly post-hoc, replication-worthy observation.
- **Flag:** this is a voice/positioning change to the Contribution paragraph
  (`main.tex:68`); left for the author to apply in their own words rather than
  imposed, since it reshapes the paper's central claim.

### 13. Breadth — **AUTHOR-DECISION (future-work framing recommended)**
- The paper already spans both modalities (genetic Perturb-seq + chemical
  Tahoe/PRISM), so the reviewer's ask is largely met. Recommended: add a short
  "Scope and breadth" sentence to §8 making the existing coverage explicit and
  scoping extensions (combinatorial perturbations, scGPT/scFoundation/Geneformer
  baselines — stubs exist in `predictors_v2_stubs.py`, verified) as future work.
- **Flag:** running additional heavyweight baselines is a compute decision; the
  future-work framing (zero compute) satisfies the reviewer's "OR explicit
  future-work scoping" branch. Recommend framing over new runs for a poster.

---

# Response-to-Reviewers Skeleton (by reviewer)

> Mapping each reviewer's stated concern → action taken. Reviewer codes from the
> decision: **yhtC**, the **workshop-fit** reviewer, **Unm6**, **usGP**.
> (Exact comment-to-reviewer attribution should be confirmed against the actual
> review text; the grouping below is by concern type as given in the audit brief.)

### Reviewer yhtC (clear-reject; integrity / presentation)
- *Fabricated Csendes title* → **Fixed** (verified real title in `refs.bib`; full
  22-entry bib re-verified against primary sources, 0 remaining issues).
- *Broken `(?)` citations (×3)* → **Resolved** (energy→Székely, Tahoe→Zhang 2025,
  conditional-impossibility→Vovk 2012; 0 undefined cites in compile).
- *Modified style file / "Title Suppressed" header* → **Fixed**; official ICML sty
  left unmodified, running-title re-asserted in preamble; header now correct.
- *Figure/table inconsistency (Fig 3 vs Table 2)* → **Fixed**; figure regenerated
  to match the headline table, no fabricated data.

### Workshop-fit reviewer (weak-reject; fit / novelty)
- *Off-the-shelf method / low novelty* → reframed contribution around the
  pre-registration protocol + sVAE+ output-head design rule; conformal machinery
  explicitly named as off-the-shelf (Item 12, author to finalize wording).
- *Fit* → the pre-registered hypothesis-testing framing (honest reported null,
  permutation thresholds, BH-FDR) is foregrounded as the workshop-relevant spine.
- *Missing comparators* → **added to response**: Wei 2025 (*Nat Methods*, 27×29
  benchmark) and TISSUE (Sun 2023, calibration→FDR) to be cited with a "why not
  subsumed" paragraph; the uncited VCC quote to be attributed (Roohani 2025, *Cell*).

### Reviewer Unm6 (accept; scientific validity)
- *Cross-cell-line CP not coverage-valid* → table caption caveated; Energy→0
  undercoverage signature foregrounded in the body; weighted-CP correctly described
  as future work (the false "exposed in the library" claim removed).
- *STATE missing from Tahoe* → omission + H1 impact stated plainly; re-run flagged
  as an infra-limited author/compute decision.

### Reviewer usGP (accept; breadth)
- *Extend perturbation types / baselines* → existing genetic + chemical breadth made
  explicit; further extensions scoped as future work (stubs already in the
  codebase). Running new baselines flagged as an author/compute decision.

---

# What still needs the author (nothing blocking camera-ready mechanics)

1. **Confirm** `results/k2_preliminary.json` is the intended frozen K2 result
   (Fig 3 + Table 2 source).
2. **Citations** — add Wei 2025, TISSUE/Sun 2023, and a VCC entry (Roohani 2025
   *Cell*) for the currently-uncited "worse than baseline on MAE" quote. (Mechanical
   once the author approves; not done here to avoid asserting a citation the author
   has not chosen to add.)
3. **Item 12** novelty-sentence reframe — author's own words.
4. **Optional, main-track-lifting:** implement weighted/Mondrian CP (Item 10c
   Path B) and a held-out K3 validation.
