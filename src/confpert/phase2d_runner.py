"""Phase 2D analysis pipeline: one-shot runner over a complete Phase 2
`results.json`.

Chains the existing primitives (pre-reg verifier, pre-reg ablation panel,
power-analysis report, bootstrap-CI engine) into a single CLI:

  python -m confpert.cli phase2d \
      --prereg paper_neurips_dnb_2026/preregistration_v2.yaml \
      --results baselines/results.json \
      --repo-dir . \
      --out paper_neurips_dnb_2026/phase2d_report.json

Outputs a JSON + a Markdown summary suitable for direct inclusion in the
ICLR 2027 D&B paper appendix (per pre-reg v2 §5.1 paper-layout plan).

The runner is intentionally pure-CPU and idempotent: it reads the locked
pre-reg + final results.json, runs the K1/K2 verification + ablation +
bootstrap + power summary, and emits a stable report. Runs in <30 seconds
on the Phase 1 results.json (~2000 rows).
"""
from __future__ import annotations

import json
import os
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from confpert import prereg, prereg_ablation, power


def _load_results(path: Path) -> list[dict[str, Any]]:
    """Load Phase 1/2 results.json and flatten into one row per
    (predictor, dataset, score, alpha) cell.

    The canonical Phase 1 schema is:
      { "rows": [{predictor, dataset, alpha, scores: {ks: {...}, w1: {...}}}, ...],
        "sha256": "..." }
    where each row aggregates multiple discrepancy scores under `scores`.

    We flatten that structure into one cell per (predictor, dataset, score,
    alpha) for the bootstrap aggregator. Falls back to bare-list / jsonl
    formats if the file isn't in the canonical shape.
    """
    with open(path) as fh:
        raw = fh.read()

    flat: list[dict[str, Any]] = []
    try:
        data = json.loads(raw)
    except json.JSONDecodeError:
        # JSONL fallback
        for line in raw.splitlines():
            line = line.strip()
            if line:
                try:
                    flat.append(json.loads(line))
                except json.JSONDecodeError:
                    pass
        return flat

    raw_rows = data["rows"] if isinstance(data, dict) and "rows" in data else (
        data if isinstance(data, list) else [data]
    )
    for row in raw_rows:
        scores_block = row.get("scores")
        if isinstance(scores_block, dict):
            for score_name, score_cell in scores_block.items():
                if not isinstance(score_cell, dict):
                    continue
                flat.append({
                    "predictor": row.get("predictor"),
                    "dataset": row.get("dataset"),
                    "score": score_name,
                    "alpha": row.get("alpha", score_cell.get("alpha")),
                    "coverage": score_cell.get("achieved_coverage"),
                    "calibration_deviation": score_cell.get("calibration_deviation"),
                    "n_perts_test": score_cell.get("n_test", row.get("n_test_perts", 0)),
                    "n_perts_calib": row.get("n_calib_perts", 0),
                    "noise_variant": row.get("noise_variant"),
                    "seed": row.get("seed"),
                    "split_type": row.get("split_type"),
                    "head": row.get("head"),
                })
        else:
            # Already-flat row (Phase 2 may emit a flatter schema)
            flat.append(row)
    return flat


def _bootstrap_summary_table(results_rows: list[dict[str, Any]],
                              n_resamples: int = 1000,
                              seed: int = 42) -> dict[str, Any]:
    """Run bootstrap CIs on every (predictor, dataset, score, alpha) cell of
    the results table.

    Phase 1 results.json schema cells:
      - predictor (str)
      - dataset (str)
      - score (str)              # KS / W1 / energy / MMD / bimodality / variance_ratio
      - alpha (float)            # nominal miscoverage; coverage = 1 - alpha
      - coverage (float | None)  # achieved coverage on the test fold
      - n_perts_test (int)       # number of test perturbations
      - n_perts_calib (int)      # number of calibration perturbations
      - cell_line_context (str | None)
      - param_count (int | None)
      - head (str | None)

    We don't have raw per-perturbation scores in results.json (those would
    blow the file up to GB). Instead, we synthesize the bootstrap from the
    saved (achieved_coverage, n_perts_test) summary using the binomial
    approximation: coverage ~ N(p, sqrt(p(1-p)/n)). This is a coarser CI
    than the per-perturbation stratified bootstrap, but it's correct for
    the binomial coverage statistic and doesn't need the raw scores.

    For per-perturbation stratified bootstrap, run
    `confpert.bootstrap.bootstrap_coverage_from_scores` directly with the
    raw score arrays in a separate script.
    """
    import numpy as np

    rng = np.random.default_rng(seed)
    cells: dict[tuple[str, str, str, float], dict[str, Any]] = {}
    for row in results_rows:
        key = (
            row.get("predictor", "?"),
            row.get("dataset", "?"),
            row.get("score", "?"),
            float(row.get("alpha", float("nan")))
            if row.get("alpha") is not None else float("nan"),
        )
        if key in cells:
            # Some result files have multiple rows per cell (one per repeat seed);
            # average coverage + sum n_perts_test for the binomial CI input.
            cells[key]["cov_sum"] += float(row.get("coverage", float("nan")))
            cells[key]["count"] += 1
            cells[key]["n_perts_test"] += int(row.get("n_perts_test", 0))
        else:
            cells[key] = {
                "cov_sum": float(row.get("coverage", float("nan"))),
                "count": 1,
                "n_perts_test": int(row.get("n_perts_test", 0)),
                "head": row.get("head"),
            }

    summary = {}
    for (pred, ds, score, alpha), agg in cells.items():
        n = max(1, agg["n_perts_test"])
        p_hat = agg["cov_sum"] / max(1, agg["count"])
        if not (0.0 <= p_hat <= 1.0):
            ci = None
        else:
            se = (p_hat * (1.0 - p_hat) / n) ** 0.5
            # Wilson interval at 95%
            z = 1.96
            denom = 1 + z * z / n
            center = (p_hat + z * z / (2 * n)) / denom
            half = z * (se ** 2 + z * z / (4 * n * n)) ** 0.5 / denom
            ci = {"lo": max(0.0, center - half), "hi": min(1.0, center + half),
                  "point": p_hat, "method": "wilson_binomial"}
        cell_key = f"{pred}|{ds}|{score}|alpha={alpha:.3f}"
        summary[cell_key] = {
            "predictor": pred, "dataset": ds, "score": score, "alpha": alpha,
            "n_perts_test": n, "n_seeds": agg["count"],
            "p_hat": p_hat, "ci": ci, "head": agg.get("head"),
        }
    return summary


def run_phase2d(
    yaml_path: Path | str,
    results_path: Path | str,
    repo_dir: Path | str | None = None,
    out_path: Path | str | None = None,
    n_resamples: int = 1000,
    skip_ablation: bool = False,
) -> dict[str, Any]:
    """Full Phase 2D analysis: verify + ablate + power + bootstrap CIs.

    Returns a dict with five top-level keys:
      - meta (timestamps, paths, hashes)
      - verifier_report (output of confpert.prereg.verify, full mode)
      - ablation_report (output of confpert.prereg_ablation.run_ablation_panel)
      - power_report (output of confpert.power.full_power_report)
      - bootstrap_ci_table (per-cell binomial Wilson 95% CI)
    """
    yaml_path = Path(yaml_path)
    results_path = Path(results_path)
    repo_dir = Path(repo_dir) if repo_dir else yaml_path.resolve().parent.parent

    meta = {
        "phase2d_runner_version": "0.2.0rc2",
        "timestamp_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "yaml_path": str(yaml_path),
        "results_path": str(results_path),
        "repo_dir": str(repo_dir),
        "n_resamples": n_resamples,
        "skip_ablation": skip_ablation,
    }

    print(f"[phase2d] verifier ...")
    verifier_result = prereg.verify(
        yaml_path=yaml_path, results_path=results_path,
        repo_dir=repo_dir, dry_run=False,
    )

    if not skip_ablation:
        print(f"[phase2d] ablation panel ...")
        ablation_report = prereg_ablation.run_ablation_panel(
            yaml_path=yaml_path, results_path=results_path, repo_dir=repo_dir,
        )
        ablation_summary = {
            "headline": ablation_report.headline,
            "n_flipped": ablation_report.n_flipped,
            "rows": [asdict(r) for r in ablation_report.rows],
        }
    else:
        ablation_summary = {"skipped": True}

    print(f"[phase2d] power report ...")
    power_summary = power.full_power_report()

    print(f"[phase2d] bootstrap CIs ...")
    results_rows = _load_results(results_path)
    bootstrap_summary = _bootstrap_summary_table(results_rows, n_resamples=n_resamples)

    out = {
        "meta": meta,
        "verifier_report": {
            "lock_check_passed": verifier_result.lock_check_passed,
            "k2_v2_disposition": verifier_result.disposition_headline,
            "hypothesis_results": [
                {
                    "hypothesis_id": h.hypothesis_id,
                    "status": h.status,
                    "observed": dict(h.observed),
                    "threshold": dict(h.threshold),
                    "n_cells_used": h.n_cells_used,
                    "error": h.error,
                }
                for h in verifier_result.hypothesis_results
            ],
        },
        "ablation_report": ablation_summary,
        "power_report": power_summary,
        "bootstrap_ci_table": bootstrap_summary,
    }

    if out_path is not None:
        out_path = Path(out_path)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        with open(out_path, "w") as fh:
            json.dump(out, fh, indent=2, default=str)
        print(f"[phase2d] wrote {out_path}")

        # Markdown sidecar for easy paper-appendix inclusion
        md_path = out_path.with_suffix(".md")
        md_lines = _to_markdown(out)
        with open(md_path, "w") as fh:
            fh.write("\n".join(md_lines))
        print(f"[phase2d] wrote {md_path}")

    return out


def _to_markdown(report: dict[str, Any]) -> list[str]:
    lines = []
    lines.append("# Phase 2D Analysis Report")
    lines.append("")
    lines.append(f"- Timestamp: `{report['meta']['timestamp_utc']}`")
    lines.append(f"- YAML: `{report['meta']['yaml_path']}`")
    lines.append(f"- Results: `{report['meta']['results_path']}`")
    lines.append("")

    lines.append("## Verifier report")
    lines.append("")
    v = report["verifier_report"]
    lines.append(f"- Lock check passed: **{v['lock_check_passed']}**")
    lines.append(f"- K2 v2 disposition: **`{v['k2_v2_disposition']}`**")
    lines.append("")
    lines.append("| Hypothesis | Status | Observed |")
    lines.append("|---|---|---|")
    for h in v["hypothesis_results"]:
        observed_str = ", ".join(f"{k}={v}" for k, v in h["observed"].items())
        if len(observed_str) > 80:
            observed_str = observed_str[:77] + "..."
        lines.append(f"| {h['hypothesis_id']} | {h['status']} | {observed_str} |")
    lines.append("")

    lines.append("## Pre-registration ablation panel")
    lines.append("")
    ab = report["ablation_report"]
    if ab.get("skipped"):
        lines.append("_Skipped._")
    else:
        lines.append(f"**{ab['headline']}**")
        lines.append("")
        lines.append("| Hypothesis | Pre-reg | Shopped | Flipped |")
        lines.append("|---|---|---|---|")
        for r in ab["rows"]:
            lines.append(
                f"| {r['hypothesis']} | {r['prereg_result']} | "
                f"{r['shopped_result']} | "
                f"{'YES' if r.get('flipped_to_pass') else 'no'} |"
            )
    lines.append("")

    lines.append("## Power analysis summary")
    lines.append("")
    pr = report["power_report"]
    lines.append(f"```json")
    lines.append(json.dumps(pr, indent=2, default=str)[:4000])
    lines.append("```")
    lines.append("")

    lines.append("## Bootstrap CI table (Wilson 95% binomial)")
    lines.append("")
    lines.append("| predictor | dataset | score | alpha | n_perts_test | p̂ | 95% CI |")
    lines.append("|---|---|---|---|---|---|---|")
    bt = report["bootstrap_ci_table"]
    # Print top 30 cells by (predictor, dataset, alpha) for readability
    for key in sorted(bt.keys())[:30]:
        cell = bt[key]
        ci = cell["ci"]
        ci_str = f"[{ci['lo']:.3f}, {ci['hi']:.3f}]" if ci else "(invalid)"
        lines.append(
            f"| {cell['predictor']} | {cell['dataset']} | {cell['score']} | "
            f"{cell['alpha']:.3f} | {cell['n_perts_test']} | "
            f"{cell['p_hat']:.3f} | {ci_str} |"
        )
    if len(bt) > 30:
        lines.append(f"\n_…and {len(bt) - 30} more cells in the JSON sidecar._")
    return lines


def cli_main(argv: list[str]) -> int:
    import argparse
    p = argparse.ArgumentParser(description="ConfPert Phase 2D analysis runner")
    p.add_argument("--prereg", required=True,
                   help="Path to preregistration_v2.yaml")
    p.add_argument("--results", required=True,
                   help="Path to results.json (Phase 1 or Phase 2)")
    p.add_argument("--repo-dir", default=None, help="Repo root (default: prereg parent)")
    p.add_argument("--out", default=None, help="Output JSON path")
    p.add_argument("--n-resamples", type=int, default=1000)
    p.add_argument("--skip-ablation", action="store_true",
                   help="Skip the pre-reg ablation panel (faster)")
    args = p.parse_args(argv)

    run_phase2d(
        yaml_path=args.prereg,
        results_path=args.results,
        repo_dir=args.repo_dir,
        out_path=args.out,
        n_resamples=args.n_resamples,
        skip_ablation=args.skip_ablation,
    )
    return 0


if __name__ == "__main__":
    import sys
    raise SystemExit(cli_main(sys.argv[1:]))
