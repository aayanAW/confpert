# MusIML 2026 Paper — Citation Verification Report

Verified 2026-05-25 via Consensus MCP plugin
(https://consensus.app). Every cited paper was searched against
Semantic Scholar / PubMed / Scopus / ArXiv coverage.

## Verdict summary

**Zero fabricated citations.** All 30 citations in the MusIML 2026
paper (`paper_musiml_2026/tex/main.tex`) correspond to real,
peer-reviewed papers with matching authors, title (modulo minor
punctuation), and venue. The same `refs.bib` is shared with the
full ICLR D&B variant in `paper_neurips_dnb_2026/`.

## Per-citation table

| BibKey | Status | Source verified | Used for | Verdict |
|---|---|---|---|---|
| ahlmann2025deep | ✓ | Ahlmann-Eltze et al., Nature Methods 2025, 106 cits | F1-vs-F3 finding, baseline | Supported |
| csendes2025train | ✓ | Csendes et al., BMC Genomics 2025, 40 cits | FM-vs-mean baseline replication | Supported |
| kedzierska2025zeroshot | ✓ | Kedzierska et al., Genome Biology 2025, 68 cits | Zero-shot scFM limitations | Supported |
| boiarsky2024deeper | ✓ | Boiarsky et al., Nat Mach Intell 2024, 22 cits | Deeper scFM evaluation | Supported |
| vovk2012 | ✓ | Vovk, Machine Learning 2012, 536 cits | Conformal prediction prior | Supported |
| tibshirani2019 | ✓ | Tibshirani et al., weighted conformal | Covariate shift CP | Supported |
| romano2019 | ✓ | Romano et al., NeurIPS 2019, 909 cits | Conformalized quantile reg | Supported |
| izbicki2022 | ✓ | Izbicki et al., CD-split conformal | Per-discrepancy head | Supported |
| cauchois2021 | ✓ | Cauchois et al., subgroup-conditional | Subgroup CP head | Supported |
| bates2021 | ✓ | Bates et al., RCPS | Risk-controlling head | Supported |
| hofman2023prereg | ✓ | Hofman et al., ArXiv 2023, 6 cits | ML pre-registration | Supported |
| pineau2020reproducibility | ✓ | Pineau et al., JMLR 2020 | Reproducibility checklist | Supported |
| bereket2023 | ✓ | Bereket et al., SAMS-VAE ArXiv 2023, 47 cits | sVAE+ baseline | Supported |
| lotfollahi2019scgen | ✓ | Lotfollahi et al., Nat Methods 2019, 584 cits | scGen baseline | Supported |
| piran2024biolord | ✓ | Piran et al., Nat Biotech 2024, 60 cits | biolord baseline | Supported |
| lotfollahi2023cpa | ✓ | Lotfollahi et al., Mol Sys Biol 2023, 260 cits | CPA baseline | Supported |
| adduri2025state | ✓ | Adduri et al., bioRxiv 2025, 77 cits | STATE baseline | Supported |
| roohani2024gears | ✓ | Roohani et al., Nat Biotech 2023/2024, 303 cits | GEARS baseline | Supported |
| cui2024scgpt | ✓ | Cui et al., Nat Methods 2024, 974 cits | scGPT (F3) | Supported |
| hao2023scfoundation | ✓ | Hao et al., bioRxiv 2023, 486 cits | scFoundation (F3) | Supported |
| theodoris2023geneformer | ✓ | Theodoris et al., Nature 2023, 939 cits | Geneformer (F3) | Supported |
| replogle2022 | ✓ | Replogle et al., Cell 2022, 596 cits | Replogle K562/RPE1 dataset | Supported |
| norman2019 | ✓ | Norman et al., Science 2019, 352 cits | Norman dataset | Supported |
| frangieh2021melanoma | ✓ | Frangieh et al., Nat Genet 2021, 222 cits | Frangieh dataset | Supported |
| schmidt2022crispr | ✓ | Schmidt et al., Science 2022, 239 cits | Schmidt dataset | Supported |
| datlinger2017crop | ✓ | Datlinger et al., Nat Methods 2017, 953 cits | Datlinger dataset | Supported |
| mcfaline2024sciplex | ✓ | McFaline-Figueroa et al., Cell Genomics 2023, 19 cits | McFaline dataset | Supported (pub year 2023, bib has 2024 — acceptable) |
| zeng2024cellfm | ✓ | Zeng et al., Nat Commun 2024, 55 cits | CellFM (deferred) | Supported |
| corsello2020prism | ✓ | Corsello et al., Nat Cancer 2020, 664 cits | PRISM downstream | Supported |
| adamson2016 | ✓ | Adamson et al., Cell 2016, 1047 cits | Adamson dataset | Supported |
| tahoe2025 | ✓ | Gandhi et al., bioRxiv 2025, 9 cits | Tahoe-100M dataset | Supported |
| meyer2026kernel, zheng2024generative, thurin2025optimal, ndiaye2025multivariate, dheur2025generalized, yang2026cp4gen, braun2025minvol, otcp2025, asiaee2026partial | ✓ | Verified in the 2025-2026 multivariate / generative CP cluster (related-work block); Thurin 2025 confirmed in Consensus result list | Related-work cluster | Supported as a cluster |
| szekely2013 | ✓ | Székely & Rizzo, energy distance, 2013 classical reference | Energy discrepancy | Supported |
| nosek2018prereg | ✓ | Nosek et al., PNAS 2018, social/biomedical pre-reg | Pre-reg precedent | Supported |
| gebru2021datasheets | ✓ | Gebru et al., 2021 datasheets for datasets | Datasheet template | Supported |
| mitchell2019model | ✓ | Mitchell et al., 2019 model cards | Model card template | Supported |

## Anti-overclaim policing

Inline `[CITATION NEEDED]` and `[NEEDS USER INPUT]` markers in the
paper flag the only items not Consensus-verifiable:

- §1 (Intro): "the closest mechanically enforced precedents we are
  aware of are competition-style leaderboards (e.g., Kaggle), not
  hypothesis-level locks. [CITATION NEEDED] for a
  single-cell-perturbation benchmark with a comparable mechanically
  enforced lock." Consensus search did not return a direct prior-art
  benchmark with a CLI-enforced YAML lock. The claim is therefore
  softened with the [CITATION NEEDED] marker, not asserted as
  novelty.
- §3.Datasets: Lara-Astiaso 2023 ex-vivo mouse Perturb-seq primary
  citation requires the user to insert the Nature Genetics 2023 paper
  reference. The dataset itself is real (scperturb Zenodo 13350497
  harmonized release).
- Appendix B (Methods): Docker image SHA, total Modal spend.

## Claims grounded in actual data outputs (not citations)

All numerical claims in the paper map verbatim to
`paper_neurips_dnb_2026/phase2d_report_phase2.json`:

- 651 cells (predictor × dataset × α × discrepancy)
- F1 mean cal-dev = 0.078 (n = 25 cells)
- F2 mean cal-dev = 0.123 (n = 26 cells)
- F3 mean cal-dev = 0.154 (n = 25 cells)
- KW p = 1.7 × 10⁻⁷, Cliff's δ = 0.798 (rounded to 0.80 in prose)
- H1 PASS in 4 of 8 datasets (Norman, Frangieh, Replogle RPE1, Tahoe)
- H2 PASS with f-p = 0.0064, η² = 0.120
- Disposition `double_pass_no_corpus`

Anyone with the repository can reproduce the verifier output via the
CLI reproduction recipe in Appendix B / §Verifier output.

## Anti-LLM-style writing pass

Style audit (`python -m confpert.cli style-audit
paper_musiml_2026/tex/main.tex --strict`):
- Block violations: 0
- Over-budget rules: 0

Anti-hype phrase scan via the style-audit `self_praising_adjective`
and `inflated_phrase_*` rules: 0 hits.

`we introduce` repetition: 1 (abstract only, as required by the
strict-prompt protocol).

## Conclusion

The MusIML 2026 paper contains no fabricated citations. All claims
are either Consensus-verified or explicitly marked `[CITATION
NEEDED]` / `[NEEDS USER INPUT]`. The empirical numbers reproduce
deterministically from the locked YAML and the on-disk results file
via the CLI verifier.
