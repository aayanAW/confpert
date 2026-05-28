"""Pre-registration ablation panel: Pipeline A vs Pipeline B counterfactual.

Per `preregistration_v2.md` §4. This module implements the empirical
demonstration that pre-registration prevents threshold-shopping bias.

Pipeline A: pre-registered thresholds (the YAML's locked values).
Pipeline B: analyst threshold-shops post-hoc to maximize pass count.
           Each hypothesis gets the smallest single-step relaxation of one
           threshold parameter that would flip FAIL → PASS, drawn from a
           pre-declared "plausible alternative thresholds" set.

The threshold-relaxation set is defined here (not in pre-reg YAML) because
Pipeline B is the *counterfactual analyst*, not the pre-committed analysis.
Analyst-defined relaxations are: looser p (0.05 → 0.10), looser ρ (0.5 → 0.4),
fewer required datasets-passing (3-of-4 → 2-of-4 or 2-of-5).

Output:
  side-by-side table {hypothesis, prereg_threshold, prereg_result,
                      shopped_threshold, shopped_result, delta}
  + headline summary string ("without pre-registration N hypotheses would
  have flipped to PASS").

If K2 v2 lands triple-null (per pre-reg v2 §2.6 outcome table), this panel
becomes the Tier-1 paper contribution.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from confpert import prereg


# Pipeline B's analyst-counterfactual threshold relaxations
# (NOT pre-registered; this is THE thing the panel demonstrates is the bias)
SHOPPED_RELAXATIONS = {
    # For Spearman-based tests (H1, H1b, H2b)
    "spearman_per_dataset": [
        # Each entry is one threshold relaxation to try, in order from gentle to aggressive
        {"name": "p<0.10 instead of p<0.05", "p_value": 0.10},
        {"name": "rho>0.4 instead of rho>0.5", "rho_threshold": 0.4},
        {"name": "n_passing>=2 of 4 instead of >=3 of 4", "n_passing": 2},
        {"name": "p<0.10 AND rho>0.4", "p_value": 0.10, "rho_threshold": 0.4},
        {"name": "p<0.10 AND n_passing>=2", "p_value": 0.10, "n_passing": 2},
        {"name": "rho>0.4 AND n_passing>=2", "rho_threshold": 0.4, "n_passing": 2},
        {"name": "all three relaxed (p<0.10, rho>0.4, n>=2)",
         "p_value": 0.10, "rho_threshold": 0.4, "n_passing": 2},
    ],
    "two_way_anova_type2": [
        {"name": "p<0.05 instead of p<0.0167", "p_value": 0.05},
        {"name": "eta_sq>0.05 instead of eta_sq>0.10", "eta_threshold": 0.05},
        {"name": "both relaxed (p<0.05 AND eta_sq>0.05)",
         "p_value": 0.05, "eta_threshold": 0.05},
    ],
    "kruskal_wallis": [
        {"name": "p<0.05 instead of p<0.0167", "p_value": 0.05},
        {"name": "delta>0.30 instead of delta>0.50", "delta_threshold": 0.30},
        {"name": "delta>0.20 instead of delta>0.50", "delta_threshold": 0.20},
        {"name": "p<0.05 AND delta>0.30", "p_value": 0.05, "delta_threshold": 0.30},
        {"name": "p<0.05 AND delta>0.20", "p_value": 0.05, "delta_threshold": 0.20},
    ],
}


@dataclass
class AblationRow:
    hypothesis: str
    test_type: str
    prereg_threshold: str
    prereg_result: str
    shopped_threshold: str
    shopped_result: str
    flipped_to_pass: bool
    observed_stats: dict[str, Any] = field(default_factory=dict)


@dataclass
class AblationReport:
    rows: list[AblationRow]
    n_flipped: int
    headline: str


def _apply_spearman_shopped_criteria(
    per_dataset_observed: dict[str, dict],
    rho_threshold: float,
    rho_operator: str,
    p_threshold: float,
    n_passing_threshold: int,
) -> tuple[bool, int]:
    """Re-evaluate Spearman per-dataset results against shopped thresholds.
    Returns (overall_pass, n_datasets_passing)."""
    n_passing = 0
    for ds, info in per_dataset_observed.items():
        if "rho" not in info:
            continue
        rho = info["rho"]
        p = info["p_value"]
        if rho_operator == ">":
            passes_rho = rho > rho_threshold
        else:
            passes_rho = rho < rho_threshold
        passes_p = p < p_threshold
        if passes_rho and passes_p:
            n_passing += 1
    return n_passing >= n_passing_threshold, n_passing


def _shop_spearman_hypothesis(
    hyp_id: str, hyp_spec: dict, prereg_result: prereg.VerificationResult
) -> AblationRow:
    """For a FAIL spearman hypothesis, find smallest relaxation that flips it."""
    # Extract pre-reg thresholds
    rho_op = ">"
    rho_thresh = 0.5
    p_thresh = 0.05
    n_passing_prereg = 3
    for crit in hyp_spec.get("success_criteria", []):
        if crit.get("test") == "rho":
            rho_op = crit.get("operator", ">")
            rho_thresh = crit.get("threshold", 0.5)
        elif crit.get("test") == "p_value":
            p_thresh = crit.get("threshold", 0.05)
        elif crit.get("test") == "n_datasets_passing":
            n_passing_prereg = crit.get("threshold", 3)

    prereg_label = (
        f"rho {rho_op} {rho_thresh}, p < {p_thresh}, "
        f">= {n_passing_prereg} datasets"
    )
    prereg_pass_label = "PASS" if prereg_result.status == "PASS" else "FAIL"

    per_ds = prereg_result.observed.get("per_dataset", {})
    if prereg_result.status == "PASS":
        return AblationRow(
            hypothesis=hyp_id, test_type="spearman_per_dataset",
            prereg_threshold=prereg_label, prereg_result=prereg_pass_label,
            shopped_threshold="(N/A — pre-reg already passes)",
            shopped_result="(N/A)", flipped_to_pass=False,
            observed_stats={"per_dataset": per_ds},
        )

    # Try relaxations in order
    for relax in SHOPPED_RELAXATIONS["spearman_per_dataset"]:
        rho_t = relax.get("rho_threshold", rho_thresh)
        p_t = relax.get("p_value", p_thresh)
        n_t = relax.get("n_passing", n_passing_prereg)
        passes, n_pass = _apply_spearman_shopped_criteria(
            per_ds, rho_t, rho_op, p_t, n_t
        )
        if passes:
            return AblationRow(
                hypothesis=hyp_id, test_type="spearman_per_dataset",
                prereg_threshold=prereg_label, prereg_result=prereg_pass_label,
                shopped_threshold=relax["name"],
                shopped_result=f"PASS ({n_pass} datasets pass shopped threshold)",
                flipped_to_pass=True,
                observed_stats={
                    "shopped_rho_threshold": rho_t,
                    "shopped_p_threshold": p_t,
                    "shopped_n_passing_threshold": n_t,
                    "n_datasets_passing_shopped": n_pass,
                    "per_dataset": per_ds,
                },
            )

    return AblationRow(
        hypothesis=hyp_id, test_type="spearman_per_dataset",
        prereg_threshold=prereg_label, prereg_result=prereg_pass_label,
        shopped_threshold="(no relaxation flips to PASS)",
        shopped_result="FAIL", flipped_to_pass=False,
        observed_stats={"per_dataset": per_ds},
    )


def _shop_anova_hypothesis(
    hyp_id: str, hyp_spec: dict, prereg_result: prereg.VerificationResult
) -> AblationRow:
    p_prereg = 0.0167
    eta_prereg = 0.10
    for crit in hyp_spec.get("success_criteria", []):
        if crit.get("test") == "f_test_pvalue":
            p_prereg = crit.get("threshold", 0.0167)
        elif crit.get("test") == "eta_squared":
            eta_prereg = crit.get("threshold", 0.10)
    prereg_label = f"F p < {p_prereg}, eta^2 > {eta_prereg}"
    prereg_pass_label = "PASS" if prereg_result.status == "PASS" else "FAIL"

    if prereg_result.status == "PASS":
        return AblationRow(
            hypothesis=hyp_id, test_type="two_way_anova_type2",
            prereg_threshold=prereg_label, prereg_result=prereg_pass_label,
            shopped_threshold="(N/A — pre-reg already passes)",
            shopped_result="(N/A)", flipped_to_pass=False,
            observed_stats=dict(prereg_result.observed),
        )

    f_pvalue = prereg_result.observed.get("f_pvalue_cell_line_context", 1.0)
    eta_sq = prereg_result.observed.get("eta_squared_cell_line_context", 0.0)
    for relax in SHOPPED_RELAXATIONS["two_way_anova_type2"]:
        p_t = relax.get("p_value", p_prereg)
        eta_t = relax.get("eta_threshold", eta_prereg)
        if f_pvalue < p_t and eta_sq > eta_t:
            return AblationRow(
                hypothesis=hyp_id, test_type="two_way_anova_type2",
                prereg_threshold=prereg_label, prereg_result=prereg_pass_label,
                shopped_threshold=relax["name"],
                shopped_result=f"PASS (observed F p={f_pvalue:.4f}, eta^2={eta_sq:.3f})",
                flipped_to_pass=True,
                observed_stats={
                    "shopped_p_threshold": p_t, "shopped_eta_threshold": eta_t,
                    **prereg_result.observed,
                },
            )

    return AblationRow(
        hypothesis=hyp_id, test_type="two_way_anova_type2",
        prereg_threshold=prereg_label, prereg_result=prereg_pass_label,
        shopped_threshold="(no relaxation flips to PASS)",
        shopped_result="FAIL", flipped_to_pass=False,
        observed_stats=dict(prereg_result.observed),
    )


def _shop_kruskal_hypothesis(
    hyp_id: str, hyp_spec: dict, prereg_result: prereg.VerificationResult
) -> AblationRow:
    p_prereg = 0.0167
    delta_prereg = 0.50
    for crit in hyp_spec.get("success_criteria", []):
        if crit.get("test") == "kw_pvalue":
            p_prereg = crit.get("threshold", 0.0167)
        elif crit.get("test") == "cliffs_delta_best_vs_worst":
            delta_prereg = crit.get("threshold", 0.50)
    prereg_label = f"KW p < {p_prereg}, |delta| > {delta_prereg}"
    prereg_pass_label = "PASS" if prereg_result.status == "PASS" else "FAIL"

    if prereg_result.status == "PASS":
        return AblationRow(
            hypothesis=hyp_id, test_type="kruskal_wallis",
            prereg_threshold=prereg_label, prereg_result=prereg_pass_label,
            shopped_threshold="(N/A — pre-reg already passes)",
            shopped_result="(N/A)", flipped_to_pass=False,
            observed_stats=dict(prereg_result.observed),
        )

    kw_p = prereg_result.observed.get("kw_pvalue", 1.0)
    delta = abs(prereg_result.observed.get("cliffs_delta", 0.0))
    for relax in SHOPPED_RELAXATIONS["kruskal_wallis"]:
        p_t = relax.get("p_value", p_prereg)
        delta_t = relax.get("delta_threshold", delta_prereg)
        if kw_p < p_t and delta > delta_t:
            return AblationRow(
                hypothesis=hyp_id, test_type="kruskal_wallis",
                prereg_threshold=prereg_label, prereg_result=prereg_pass_label,
                shopped_threshold=relax["name"],
                shopped_result=f"PASS (observed KW p={kw_p:.4f}, |delta|={delta:.3f})",
                flipped_to_pass=True,
                observed_stats={
                    "shopped_p_threshold": p_t, "shopped_delta_threshold": delta_t,
                    **prereg_result.observed,
                },
            )

    return AblationRow(
        hypothesis=hyp_id, test_type="kruskal_wallis",
        prereg_threshold=prereg_label, prereg_result=prereg_pass_label,
        shopped_threshold="(no relaxation flips to PASS)",
        shopped_result="FAIL", flipped_to_pass=False,
        observed_stats=dict(prereg_result.observed),
    )


def run_ablation_panel(
    yaml_path,
    results_path,
    repo_dir,
) -> AblationReport:
    """Top-level ablation entry point.

    1. Run prereg.verify (full mode) to get Pipeline A results.
    2. For each FAIL hypothesis, find the smallest SHOPPED_RELAXATIONS entry
       that flips it to PASS. If none does, report unflippable.
    3. Aggregate into AblationReport.
    """
    pa_result = prereg.verify(
        yaml_path=yaml_path, results_path=results_path,
        repo_dir=repo_dir, dry_run=False,
    )

    # Load YAML to get per-hypothesis specs for threshold extraction
    spec = prereg._load_yaml(yaml_path)
    hyp_specs = spec.get("k2_v2", {}).get("hypotheses", {})

    rows = []
    for hr in pa_result.hypothesis_results:
        hyp_id = hr.hypothesis_id
        hyp_spec = hyp_specs.get(hyp_id, {})
        test_kind = hyp_spec.get("test", "")
        if test_kind == "spearman_per_dataset":
            rows.append(_shop_spearman_hypothesis(hyp_id, hyp_spec, hr))
        elif test_kind == "two_way_anova_type2":
            rows.append(_shop_anova_hypothesis(hyp_id, hyp_spec, hr))
        elif test_kind == "kruskal_wallis":
            rows.append(_shop_kruskal_hypothesis(hyp_id, hyp_spec, hr))
        else:
            rows.append(AblationRow(
                hypothesis=hyp_id, test_type=test_kind or "unknown",
                prereg_threshold="(unknown)", prereg_result=hr.status,
                shopped_threshold="(not implemented for this test type)",
                shopped_result="(N/A)", flipped_to_pass=False,
            ))

    n_flipped = sum(1 for r in rows if r.flipped_to_pass)
    if n_flipped > 0:
        headline = (
            f"Pre-registration prevented {n_flipped} of {len(rows)} hypotheses "
            f"from being reported as PASS via post-hoc threshold relaxation. "
            f"Without pre-registration, the headline would have included "
            f"{n_flipped} additional PASS claims."
        )
    else:
        if all(r.prereg_result == "PASS" for r in rows):
            headline = "Pre-reg passes on all hypotheses; ablation Δ=0."
        else:
            headline = (
                "No relaxation in the SHOPPED_RELAXATIONS set flips any FAIL "
                "hypothesis to PASS. Pre-reg result is robust against the "
                "considered counterfactual analyst thresholds. Stronger "
                "demonstration of pre-reg discipline."
            )
    return AblationReport(rows=rows, n_flipped=n_flipped, headline=headline)


def _format_text_report(report: AblationReport) -> str:
    lines = [
        "=" * 78,
        "Pre-Registration Ablation Panel (Pipeline A vs Pipeline B)",
        "=" * 78,
        "",
        f"{'Hypothesis':<28} {'Pre-reg':>8} {'Shopped':>8} {'Flipped':>8}",
        "-" * 78,
    ]
    for r in report.rows:
        flipped_marker = "YES <-" if r.flipped_to_pass else "no"
        lines.append(f"{r.hypothesis:<28} {r.prereg_result:>8} "
                     f"{r.shopped_result[:8]:>8} {flipped_marker:>8}")
        lines.append(f"  pre-reg: {r.prereg_threshold}")
        lines.append(f"  shopped: {r.shopped_threshold}")
        lines.append("")
    lines.append("-" * 78)
    lines.append(report.headline)
    lines.append("=" * 78)
    return "\n".join(lines)


def cli_main(argv: list[str]) -> int:
    import argparse
    import json
    from dataclasses import asdict
    from pathlib import Path

    parser = argparse.ArgumentParser(prog="confpert prereg-ablate")
    parser.add_argument("--prereg", type=Path, required=True)
    parser.add_argument("--results", type=Path, default=Path("baselines/results.json"))
    parser.add_argument("--repo-dir", type=Path, default=Path.cwd())
    parser.add_argument("--out", type=Path, default=None)
    args = parser.parse_args(argv)

    report = run_ablation_panel(
        yaml_path=args.prereg, results_path=args.results, repo_dir=args.repo_dir,
    )
    print(_format_text_report(report))
    if args.out:
        with open(args.out, "w") as fh:
            json.dump({
                "headline": report.headline,
                "n_flipped": report.n_flipped,
                "rows": [asdict(r) for r in report.rows],
            }, fh, indent=2, default=str)
        print(f"\nJSON report written to {args.out}")
    return 0
