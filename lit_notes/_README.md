# ConfPert literature notes

Phase 0 deliverable: one note file per paper in the project bibliography, plus `_synthesis.md`.

Note file format (per paper):
- **Claim** (one sentence)
- **Method** (compact technical summary, equations where load-bearing)
- **Datasets** (concrete numbers)
- **Metrics** (exact metric definitions)
- **What we steal** (specific ideas/equations to import into ConfPert)
- **What we wrap** (specific predictors / scoring functions to wrap as one of our nine wrappable predictors)
- **Failure modes** (caveats, limitations, where it breaks)
- **Code URLs** (canonical repo, checkpoint hosting)

Reading order matches the Phase 0 spec:

## Tier 1: Conformal for distributions (core method)
1. `chernozhukov_2021_dcp.md` - Chernozhukov, Wuthrich, Yinchu PNAS 2021
2. `angelopoulos_bates_2021_gentle.md` - Angelopoulos, Bates 2021
3. `bates_2021_rcps.md` - Bates et al. JACM 2021 (risk-controlling prediction sets)
4. `romano_2019_cqr.md` - Romano, Patterson, Candes NeurIPS 2019
5. `izbicki_2022_cd_hpd.md` - Izbicki, Shimizu, Stern JMLR 2022
6. `mvcp_gaussian_2025.md` - arXiv:2507.20941
7. `otcp_2025.md` - arXiv:2502.03609
8. `wcp_shift_2025.md` - arXiv:2501.13430
9. `tibshirani_2019_covariate_shift.md` - Tibshirani et al. NeurIPS 2019

## Tier 2: Wrappable predictors (reproduce these)
10. `cpa_2023.md` - Lotfollahi et al. Mol Syst Biol 2023
11. `scgen_2019.md` - Lotfollahi et al. Nat Methods 2019
12. `cellot_2023.md` - Bunne et al. Nat Methods 2023
13. `gears_2024.md` - Roohani et al. Nat Biotechnol 2024
14. `scgpt_2024.md` - Cui et al. Nat Methods 2024
15. `biolord_2024.md` - Piran et al. Nat Biotechnol 2024
16. `svaeplus_2023.md` - Bereket and Karaletsos NeurIPS 2023
17. `state_2025.md` - Adduri et al. Arc Institute 2025
18. `scdfm_2026.md` - Yu et al. ICLR 2026 (closest competitor, full read)
19. `cellflow_2025.md` - Klein, Fleck, Theis et al. bioRxiv 2025
20. `squidiff_2025.md` - He et al. Nat Methods 2025
21. `perturbdiff_2026.md` - Yuan et al. arXiv:2602.19685
22. `unlasting_2025.md` - Chi et al. arXiv:2506.21107

## Tier 3: Benchmarks and skeptics (positioning)
23. `perturbench_2024.md` - Wu et al. NeurIPS Datasets 2024
24. `cell_eval_2025.md` - Cell-Eval (in STATE paper)
25. `ahlmann_2025.md` - Ahlmann-Eltze, Huber, Anders 2024/2025
26. `csendes_2025.md` - Csendes et al. 2025
27. `ramakrishnan_2025.md` - Ramakrishnan et al. 2025
28. `vcc_2025.md` - Virtual Cell Challenge Cell commentary 2025

## Tier 4: Perturb-seq data and downstream
29. `norman_2019.md` - Norman et al. Science 2019
30. `replogle_2022.md` - Replogle et al. Cell 2022
31. `adamson_2016.md` - Adamson et al. Cell 2016
32. `tahoe100m_2025.md` - Tahoe-100M Vevo/Arc/Biohub 2025
33. `prism_2020.md` - Corsello et al. Nat Cancer 2020
34. `depmap_2017.md` - Tsherniak et al. Cell 2017

Then `_synthesis.md` covering: gap proof, predictor selection, conformal-method choice, narrative arc, top three risks.

Stop condition: no project code until `_synthesis.md` is approved.
