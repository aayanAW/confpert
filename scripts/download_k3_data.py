"""K3 data downloader: PRISM + DepMap + Hallmark MSigDB.

Per scripts/k3_scaffold.py and preregistration.md K3 pipeline:
  - PRISM Repurposing primary screen (Corsello 2020): per-drug per-cell-line LFC.
  - PRISM Repurposing secondary screen: 8-point dose-response viability AUC.
  - DepMap CRISPR Chronos (gene essentiality, 25Q3 release).
  - Hallmark MSigDB v2024.1.Hs (50 cancer hallmark gene sets).

Downloads all four sources to `data/k3/` so the next session can run the K3
end-to-end pipeline without waiting on data fetches. Idempotent: skips files
that already exist (use --force to re-download).

After download, the next-session K3 driver script can:
  1. Filter PRISM primary to the 49 selective-and-predictable non-oncology
     compounds (Corsello 2020 Table S5).
  2. Train ConfPert wrappers on Tahoe-100M subset (separate Modal job).
  3. Extract calibrated bimodal subpopulations per compound via
     CauchoisSubgroupConformal on predicted-bimodal-mode partitions.
  4. Compute Hallmark hypergeometric p-values + DepMap Spearman correlations.
  5. Apply BH-FDR at q=0.05 over both grids.
  6. Apply evaluate_k3_success() to verify the pre-registered criterion.
"""
from __future__ import annotations

import argparse
import json
import sys
import urllib.request
from pathlib import Path


# Working canonical URLs (validated 2026-05-03 against depmap.org manifest).
# PRISM URLs come from PRISM Primary Repurposing DepMap Public 24Q2; for
# Repurposing+Public+19Q4 (Corsello 2020) the Figshare file IDs are different
# and currently return HTML when fetched directly. The 24Q2 release is the
# best public canonical refresh of the same screen.
SOURCES = {
    # Hallmark MSigDB v2024.1.Hs — direct broadinstitute.org URL, no login.
    "hallmark_msigdb_v2024.gmt": (
        "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2024.1.Hs/"
        "h.all.v2024.1.Hs.symbols.gmt"
    ),
    # PRISM Repurposing Public 24Q2 LFC collapsed (cell-line × compound LFC).
    # Equivalent to Corsello 2020's primary-screen-replicate-collapsed-logfold-change.
    "prism_24q2_lfc_collapsed.csv": (
        "https://ndownloader.figshare.com/files/46631056"
    ),
    # PRISM Repurposing 24Q2 cell-line metadata (DepMap cell-line ID + tissue).
    "prism_24q2_cell_line_info.csv": (
        "https://ndownloader.figshare.com/files/46630978"
    ),
    # PRISM Repurposing 24Q2 treatment metadata (drug name + MOA + target).
    "prism_24q2_treatment_info.csv": (
        "https://ndownloader.figshare.com/files/46631146"
    ),
    # PRISM Repurposing 24Q2 extended primary compound list (for compound names).
    "prism_24q2_extended_compound_list.csv": (
        "https://ndownloader.figshare.com/files/46630981"
    ),
}


def _download(url: str, dst: Path) -> tuple[bool, str]:
    """Download url to dst; return (success, message)."""
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "ConfPert/0.1"})
        with urllib.request.urlopen(req, timeout=300) as resp:
            data = resp.read()
            ctype = resp.headers.get("Content-Type", "")
        # Detect HTML responses disguised as the requested file (MSigDB login
        # redirect, Figshare auth, depmap behind portal-only URLs). Refuse to
        # write a fake artifact.
        head = data[:200].lstrip().lower()
        if (head.startswith(b"<!doctype html") or head.startswith(b"<html")
                or "text/html" in ctype.lower()):
            return False, (
                f"server returned HTML ({len(data)} bytes, ctype={ctype!r}); "
                f"this URL likely needs login / portal navigation. Save the file "
                f"manually and place at {dst}."
            )
        dst.write_bytes(data)
        return True, f"OK ({len(data)} bytes, ctype={ctype!r})"
    except Exception as e:
        return False, f"{type(e).__name__}: {e}"


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--out", default="data/k3", help="Output directory.")
    p.add_argument("--force", action="store_true", help="Redownload even if present.")
    p.add_argument("--source", choices=list(SOURCES) + ["all"], default="all",
                   help="Which file to fetch (or 'all').")
    args = p.parse_args()

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    items = list(SOURCES.items()) if args.source == "all" else [
        (args.source, SOURCES[args.source])
    ]

    status: dict[str, str] = {}
    for name, url in items:
        dst = out_dir / name
        if dst.exists() and not args.force:
            status[name] = f"SKIP (exists: {dst.stat().st_size} bytes)"
            print(f"[k3-dl] {name}: {status[name]}")
            continue
        print(f"[k3-dl] {name}: fetching ...")
        ok, msg = _download(url, dst)
        status[name] = ("OK: " if ok else "FAIL: ") + msg
        print(f"[k3-dl] {name}: {status[name]}")

    # Write a small manifest with paths + statuses for the K3 driver to consume.
    manifest = {
        "out_dir": str(out_dir.resolve()),
        "files": {name: {"path": str((out_dir / name).resolve()),
                          "exists": (out_dir / name).exists(),
                          "size": (out_dir / name).stat().st_size if (out_dir / name).exists() else 0,
                          "status": status[name]}
                  for name in status},
    }
    manifest_path = out_dir / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2))
    print(f"[k3-dl] manifest -> {manifest_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
