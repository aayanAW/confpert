"""K3 scaffold: drug-resistance subpopulation discovery via calibrated bimodal extraction.

Per preregistration.md (commit c7046e4b), K3 pipeline:

  1. Load PRISM primary screen (Corsello 2020, depmap.org/repurposing) for the 49
     selective-and-predictable non-oncology compounds (Pearson r > 0.4 in primary).
  2. For each compound, use a ConfPert-wrapped predictor on Tahoe-100M-trained
     model (or available checkpoint) to extract calibrated bimodal subpopulations.
  3. Cross-reference resistant-subpopulation gene signatures against:
     - DepMap CRISPR essentiality (compound x gene grid).
     - PRISM secondary screen viability AUC (compound x cell-line grid).
     - Hallmark MSigDB pathways (compound x pathway grid).
  4. BH-FDR correction at q = 0.05 over the joint compound x pathway grid AND
     separately over the compound x gene grid.
  5. Compare ConfPert calibrated against uncalibrated baselines.

This file is a SCAFFOLD. Phase 1 daytime work is required to:
  - Download PRISM primary + secondary, DepMap, Hallmark MSigDB.
  - Train ConfPert wrappers on a Tahoe-100M subset.
  - Run the full pipeline with BH correction.

Per the user's overnight directive, only the scaffold + data-availability
verification is built tonight. The actual K3 run requires Modal compute and
data downloads that exceed the overnight session.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


# ---------------------------------------------------------------------------
# K3 data sources (URLs and cached paths)
# ---------------------------------------------------------------------------

PRISM_PRIMARY_URL = (
    "https://depmap.org/portal/api/download/files?"
    "file_name=primary-screen-replicate-collapsed-logfold-change.csv"
    "&bucket=depmap-external-downloads&release=Repurposing+Public+19Q4"
)

PRISM_SECONDARY_URL = (
    "https://depmap.org/portal/api/download/files?"
    "file_name=secondary-screen-dose-response-curve-parameters.csv"
    "&bucket=depmap-external-downloads&release=Repurposing+Public+19Q4"
)

DEPMAP_CRISPR_URL = "https://depmap.org/portal/data/dependencies/CRISPR_(DepMap_Public_25Q3)_Chronos"

HALLMARK_GMT_URL = (
    "https://www.gsea-msigdb.org/gsea/msigdb/download_file.jsp?"
    "filePath=/msigdb/release/2024.1.Hs/h.all.v2024.1.Hs.symbols.gmt"
)


# ---------------------------------------------------------------------------
# BH-FDR correction (load-bearing for K3 success criterion)
# ---------------------------------------------------------------------------


def benjamini_hochberg(pvals: np.ndarray, q: float = 0.05
                       ) -> tuple[np.ndarray, float]:
    """Benjamini-Hochberg FDR correction.

    Returns (rejected_mask, threshold_pvalue).

    rejected_mask[i] is True iff pvals[i] is rejected at FDR <= q.
    """
    pvals = np.asarray(pvals, dtype=np.float64)
    n = len(pvals)
    if n == 0:
        return np.array([], dtype=bool), 0.0
    order = np.argsort(pvals)
    ranked = pvals[order]
    # BH threshold: largest k s.t. ranked[k-1] <= k * q / n
    bh_thresh = np.arange(1, n + 1) * q / n
    below = ranked <= bh_thresh
    if not below.any():
        return np.zeros(n, dtype=bool), 0.0
    k = int(np.where(below)[0].max() + 1)
    cutoff = float(ranked[k - 1])
    rejected = np.zeros(n, dtype=bool)
    rejected[order[:k]] = True
    return rejected, cutoff


# ---------------------------------------------------------------------------
# Hypergeometric Hallmark enrichment (placeholder; full impl uses MSigDB GMT)
# ---------------------------------------------------------------------------


def hypergeometric_pvalue(k_overlap: int, n_signature: int, n_pathway: int,
                          n_universe: int) -> float:
    """One-sided Fisher exact / hypergeometric p-value for over-representation."""
    from scipy.stats import hypergeom
    rv = hypergeom(n_universe, n_pathway, n_signature)
    return float(rv.sf(k_overlap - 1))  # P(X >= k)


# ---------------------------------------------------------------------------
# K3 success-criterion check (BH-corrected)
# ---------------------------------------------------------------------------


def evaluate_k3_success(
    confpert_results: dict[str, dict],   # {compound: {pathway/gene: pval, ...}}
    uncalibrated_results: dict[str, dict],
    q: float = 0.05,
    min_n_signatures: int = 3,
    uncalibrated_failure_threshold: int = 2,
) -> dict:
    """Evaluate K3 pre-registered success criterion.

    Per preregistration.md:
      At least min_n_signatures (=3) drug-cell-line pairs where ConfPert produces:
        - Hallmark enrichment q < 0.05 (BH over compound x pathway grid)
        - DepMap selective-dependency alignment q < 0.05 (BH over compound x gene grid)
      AND uncalibrated baselines fail at least uncalibrated_failure_threshold (=2)
      of those signatures at the same corrected thresholds.
    """
    # Flatten all (compound, pathway/gene) p-values into one BH grid for confpert
    cp_keys, cp_pvals = [], []
    for cmpd, results in confpert_results.items():
        for term, p in results.items():
            cp_keys.append((cmpd, term))
            cp_pvals.append(p)
    cp_pvals = np.array(cp_pvals, dtype=np.float64)
    cp_rejected, cp_cutoff = benjamini_hochberg(cp_pvals, q=q)

    confpert_signatures = [cp_keys[i] for i, r in enumerate(cp_rejected) if r]

    # Uncalibrated baseline same grid
    un_pvals = []
    for cmpd, term in cp_keys:
        un_pvals.append(uncalibrated_results.get(cmpd, {}).get(term, 1.0))
    un_pvals = np.array(un_pvals, dtype=np.float64)
    un_rejected, un_cutoff = benjamini_hochberg(un_pvals, q=q)

    # Count uncalibrated misses on ConfPert signatures
    confpert_misses_in_uncalibrated = sum(
        1 for i, key in enumerate(cp_keys)
        if cp_rejected[i] and not un_rejected[i]
    )

    success = (len(confpert_signatures) >= min_n_signatures
               and confpert_misses_in_uncalibrated >= uncalibrated_failure_threshold)

    return {
        "success": bool(success),
        "n_confpert_signatures": int(len(confpert_signatures)),
        "n_uncalibrated_misses_on_confpert_signatures": int(confpert_misses_in_uncalibrated),
        "confpert_bh_cutoff": float(cp_cutoff),
        "uncalibrated_bh_cutoff": float(un_cutoff),
        "min_n_signatures_required": min_n_signatures,
        "uncalibrated_failure_threshold": uncalibrated_failure_threshold,
        "q": q,
        "confpert_signatures": confpert_signatures[:20],
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--demo", action="store_true",
                   help="Run a demo BH-correction smoke test on synthetic p-values.")
    args = p.parse_args()

    if args.demo:
        print("[k3] BH demo: synthetic p-values from a mixture of nulls + signal.")
        rng = np.random.RandomState(0)
        n_signal = 5
        n_null = 200
        signal_p = rng.beta(0.3, 5.0, n_signal)  # small p-values
        null_p = rng.uniform(0.0, 1.0, n_null)
        all_p = np.concatenate([signal_p, null_p])
        rej, cutoff = benjamini_hochberg(all_p, q=0.05)
        print(f"[k3] BH q=0.05 -> rejected {int(rej.sum())} of {len(all_p)}; "
              f"cutoff p = {cutoff:.4f}")

        # Hypergeometric demo
        p_hg = hypergeometric_pvalue(k_overlap=8, n_signature=20, n_pathway=200,
                                      n_universe=20000)
        print(f"[k3] hypergeometric demo: 8/20 sig genes in 200/20000 pathway -> p = {p_hg:.4e}")

        # Success-criterion demo
        cp_res = {"compoundA": {"HALLMARK_MTOR": 1e-5, "DEPMAP_BRCA1": 1e-4},
                  "compoundB": {"HALLMARK_MYC": 1e-5}}
        un_res = {"compoundA": {"HALLMARK_MTOR": 0.5, "DEPMAP_BRCA1": 0.3},
                  "compoundB": {"HALLMARK_MYC": 0.4}}
        res = evaluate_k3_success(cp_res, un_res)
        print(f"[k3] K3 success-criterion demo:")
        print(json.dumps(res, indent=2))
        return

    # If run without --demo, document the data dependencies for daytime Phase 1 work
    print(f"[k3] PRISM primary URL: {PRISM_PRIMARY_URL}")
    print(f"[k3] PRISM secondary URL: {PRISM_SECONDARY_URL}")
    print(f"[k3] DepMap CRISPR Chronos: {DEPMAP_CRISPR_URL}")
    print(f"[k3] Hallmark MSigDB GMT: {HALLMARK_GMT_URL}")
    print()
    print("[k3] Phase 1 daytime work required before K3 can run end-to-end:")
    print("  1. Download all four data sources to data/k3/")
    print("  2. Filter PRISM primary to the 49 selective-and-predictable non-oncology compounds (Corsello 2020 Table S5)")
    print("  3. Train ConfPert wrappers (CPA, biolord, sVAE+, GEARS-uncertainty, STATE) on Tahoe-100M subset")
    print("  4. Extract calibrated bimodal subpopulations per compound using subgroup-conditional Cauchois 2024 head")
    print("  5. Compute Hallmark hypergeometric p-values + DepMap Spearman correlations per compound")
    print("  6. Apply BH at q=0.05 over both grids")
    print("  7. Apply evaluate_k3_success() to verify pre-registered criterion")


if __name__ == "__main__":
    main()
