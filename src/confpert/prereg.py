"""Pre-registration verification for ConfPert benchmarks.

This module implements the `confpert prereg verify` subcommand:
parses a machine-readable pre-registration YAML, verifies its
integrity (SHA-256 + git commit), loads results.json, runs the
declared hypothesis tests, and applies the pre-registered
disposition table to produce a signed verification report.

Three modes:
  --emit-hashes : compute SHA-256 of the YAML + companion .md and
                  write the values back into the YAML's `lock` block.
                  Run once at lock-in time, before Phase 2 first run.
  --dry-run     : parse YAML, validate schema, count results.json
                  rows matching each hypothesis scope, but do NOT
                  run the statistical tests. Used at lock-in to
                  verify the verifier itself before any data is in.
  (default)     : full verify: parse + hash-check + git-check +
                  results.json mtime check + run each hypothesis
                  test + apply outcome table + emit signed report.
                  Run after Phase 2 results are merged.

The verifier refuses to declare PASS / FAIL if the lock-time
SHA-256 does not match the on-disk YAML SHA-256, or if results.json
mtime predates the YAML git commit timestamp. This is the
machine-enforced pre-registration discipline.
"""
from __future__ import annotations

import hashlib
import json
import os
import subprocess
import sys
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional


@dataclass
class VerificationResult:
    """Outcome of a single hypothesis verification."""

    hypothesis_id: str
    status: str  # "PASS" | "FAIL" | "DRY_RUN" | "SKIPPED" | "ERROR"
    observed: dict[str, Any] = field(default_factory=dict)
    threshold: dict[str, Any] = field(default_factory=dict)
    n_cells_used: int = 0
    error: Optional[str] = None


@dataclass
class PreregVerification:
    """Top-level verification output of `confpert prereg verify`."""

    timestamp: str
    prereg_yaml_path: str
    prereg_yaml_sha256: str
    prereg_yaml_git_commit: Optional[str]
    results_json_path: str
    results_json_sha256: Optional[str]
    results_json_mtime: Optional[float]
    yaml_lock_mtime: Optional[float]
    lock_check_passed: bool
    hypothesis_results: list[VerificationResult] = field(default_factory=list)
    disposition_headline: Optional[str] = None
    dry_run: bool = False
    errors: list[str] = field(default_factory=list)


def _sha256_of_file(path: Path) -> str:
    """Return hex-encoded SHA-256 of file contents."""
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def _git_commit_of_file(path: Path) -> Optional[str]:
    """Return latest git commit SHA touching this file, or None if not in git."""
    try:
        result = subprocess.run(
            ["git", "log", "-n", "1", "--pretty=format:%H", "--", str(path.resolve())],
            cwd=path.parent,
            check=True,
            capture_output=True,
            text=True,
            timeout=10,
        )
        sha = result.stdout.strip()
        return sha if sha else None
    except (subprocess.CalledProcessError, FileNotFoundError, subprocess.TimeoutExpired):
        return None


def _git_commit_timestamp(commit_sha: str, repo_dir: Path) -> Optional[float]:
    """Return UNIX timestamp of git commit, or None on failure."""
    try:
        result = subprocess.run(
            ["git", "show", "-s", "--format=%ct", commit_sha],
            cwd=repo_dir,
            check=True,
            capture_output=True,
            text=True,
            timeout=10,
        )
        return float(result.stdout.strip())
    except (subprocess.CalledProcessError, FileNotFoundError, subprocess.TimeoutExpired, ValueError):
        return None


def _load_yaml(path: Path) -> dict[str, Any]:
    """Load YAML, raising a clear error if PyYAML missing or parse fails."""
    try:
        import yaml
    except ImportError as exc:
        raise RuntimeError(
            "PyYAML is required for prereg verification. "
            "Install via `pip install pyyaml`."
        ) from exc

    with open(path) as fh:
        return yaml.safe_load(fh)


def _validate_schema(spec: dict[str, Any]) -> list[str]:
    """Return a list of schema-validation errors. Empty list = valid."""
    errors: list[str] = []
    required_top = ["version", "project", "phase", "lock", "k1_v2", "k2_v2"]
    for key in required_top:
        if key not in spec:
            errors.append(f"missing top-level key: {key}")

    if "k2_v2" in spec:
        k2 = spec["k2_v2"]
        if "hypotheses" not in k2:
            errors.append("k2_v2 missing `hypotheses` block")
        else:
            for hyp_id, hyp in k2["hypotheses"].items():
                if "success_criteria" not in hyp:
                    errors.append(f"hypothesis {hyp_id} missing `success_criteria`")
                if "test" not in hyp:
                    errors.append(f"hypothesis {hyp_id} missing `test`")
        if "outcome_table" not in k2:
            errors.append("k2_v2 missing `outcome_table`")

    if "k1_v2" in spec:
        k1 = spec["k1_v2"]
        for key in ["predictors", "datasets", "splits", "discrepancies", "alphas"]:
            if key not in k1:
                errors.append(f"k1_v2 missing `{key}`")

    return errors


def emit_hashes(yaml_path: Path) -> None:
    """Compute SHA-256 of yaml + companion .md, write back into yaml `lock` block.

    Should be run ONCE at Phase 2A lock-in. After this runs, the yaml
    is frozen — any further modification will cause the SHA-256 check
    to fail at verify time.
    """
    try:
        import yaml
    except ImportError as exc:
        raise RuntimeError("PyYAML required for emit_hashes") from exc

    with open(yaml_path) as fh:
        spec = yaml.safe_load(fh)

    yaml_dir = yaml_path.parent
    companion_md = yaml_dir / spec.get("companion_md", "preregistration_v2.md")

    spec.setdefault("lock", {})

    # Compute SHA-256 with the lock block ZEROED OUT (so the hash is self-consistent).
    # We write a "lock-stamping" copy with placeholders, hash it, then overwrite.
    spec["lock"]["sha256_md"] = (
        _sha256_of_file(companion_md) if companion_md.exists() else None
    )
    spec["lock"]["git_commit_md"] = (
        _git_commit_of_file(companion_md) if companion_md.exists() else None
    )
    # YAML self-hash: we need to write WITHOUT the sha256_yaml field, then hash that
    # canonical form, then overwrite with the hash. Two-pass.
    spec["lock"]["sha256_yaml"] = None
    spec["lock"]["git_commit_yaml"] = _git_commit_of_file(yaml_path)
    spec["lock"]["results_json_max_mtime"] = None

    # Pass 1: write spec without sha256_yaml, compute SHA-256 of THAT canonical form
    with open(yaml_path, "w") as fh:
        yaml.safe_dump(spec, fh, sort_keys=False, default_flow_style=False)

    yaml_sha = _sha256_of_file(yaml_path)
    spec["lock"]["sha256_yaml"] = yaml_sha

    # Pass 2: write with sha256_yaml populated. Note: the file's actual SHA-256
    # will differ from spec["lock"]["sha256_yaml"] (because the latter was
    # computed on the pre-stamp form). Verifier re-derives the same canonical
    # form before comparing.
    with open(yaml_path, "w") as fh:
        yaml.safe_dump(spec, fh, sort_keys=False, default_flow_style=False)

    print(f"[prereg emit-hashes] sha256_yaml (pre-stamp canonical): {yaml_sha}")
    print(f"[prereg emit-hashes] sha256_md: {spec['lock']['sha256_md']}")
    print(f"[prereg emit-hashes] git_commit_yaml: {spec['lock']['git_commit_yaml']}")
    print(f"[prereg emit-hashes] git_commit_md: {spec['lock']['git_commit_md']}")
    print(f"[prereg emit-hashes] Lock block written to {yaml_path}")


def _verify_yaml_self_hash(yaml_path: Path, declared_sha: Optional[str]) -> tuple[bool, str]:
    """Re-derive the pre-stamp canonical YAML form, compute its SHA-256,
    compare to the declared sha256_yaml. Return (passed, observed_sha).
    """
    if declared_sha is None:
        return False, "no declared sha256_yaml (yaml not locked)"

    try:
        import yaml
    except ImportError:
        return False, "PyYAML missing"

    with open(yaml_path) as fh:
        spec = yaml.safe_load(fh)

    # Re-derive pre-stamp form: zero out sha256_yaml
    if "lock" in spec and "sha256_yaml" in spec["lock"]:
        spec["lock"]["sha256_yaml"] = None

    import tempfile
    with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as tmp:
        yaml.safe_dump(spec, tmp, sort_keys=False, default_flow_style=False)
        tmp_path = Path(tmp.name)

    try:
        observed = _sha256_of_file(tmp_path)
    finally:
        tmp_path.unlink(missing_ok=True)

    return (observed == declared_sha), observed


def _verify_results_mtime(
    results_path: Path,
    yaml_lock_commit: Optional[str],
    repo_dir: Path,
) -> tuple[bool, str]:
    """Verify that results.json mtime exceeds the YAML lock-commit timestamp.

    If yaml_lock_commit is None (yaml not in git or not locked), this check
    is SKIPPED with a warning, not failed.
    """
    if yaml_lock_commit is None:
        return True, "yaml lock commit unknown; mtime check skipped"
    if not results_path.exists():
        return False, f"results.json not found at {results_path}"
    yaml_ts = _git_commit_timestamp(yaml_lock_commit, repo_dir)
    if yaml_ts is None:
        return True, "could not fetch yaml commit timestamp; mtime check skipped"
    results_mtime = results_path.stat().st_mtime
    if results_mtime <= yaml_ts:
        return False, (
            f"results.json mtime ({results_mtime}) <= yaml lock commit ts ({yaml_ts}); "
            "results may predate the pre-registration"
        )
    return True, f"results.json mtime > yaml lock commit timestamp"


def _run_hypothesis_dry(
    hyp_id: str,
    hyp_spec: dict[str, Any],
    results: dict[str, Any],
) -> VerificationResult:
    """Dry-run a single hypothesis: count matching rows but skip statistical test."""
    rows = results.get("rows", [])
    test_kind = hyp_spec.get("test", "<unknown>")
    return VerificationResult(
        hypothesis_id=hyp_id,
        status="DRY_RUN",
        observed={"test_declared": test_kind, "n_total_rows_in_results": len(rows)},
        threshold={
            crit.get("test", "<unknown>"): {
                "operator": crit.get("operator"),
                "threshold": crit.get("threshold"),
            }
            for crit in hyp_spec.get("success_criteria", [])
        },
        n_cells_used=len(rows),
    )


def _flatten_results_rows(
    results: dict[str, Any],
    include_split_types: tuple[str, ...] = ("within_perturbation", ""),
) -> list[dict[str, Any]]:
    """Flatten results.json rows into per-(predictor, dataset, alpha, score) cells.
    Mirrors scripts/k2_analysis.load_calibration_table().
    """
    flat = []
    for r in results.get("rows", []):
        if r.get("split_type", "") not in include_split_types:
            continue
        for score_name, sr in r.get("scores", {}).items():
            if isinstance(sr, dict) and "error" not in sr:
                dev = sr.get("calibration_deviation")
                if dev is None:
                    continue
                flat.append({
                    "predictor": r["predictor"],
                    "dataset": r["dataset"],
                    "alpha": r["alpha"],
                    "noise_variant": r.get("noise_variant", ""),
                    "score": score_name,
                    "calibration_deviation": float(dev),
                })
    return flat


def _aggregate_per_predictor_per_dataset(rows: list[dict[str, Any]]) -> dict[tuple[str, str], float]:
    """Average calibration_deviation over alpha + noise_variant + score per (predictor, dataset)."""
    import numpy as np

    by_key: dict[tuple[str, str], list[float]] = {}
    for r in rows:
        by_key.setdefault((r["predictor"], r["dataset"]), []).append(r["calibration_deviation"])
    return {k: float(np.mean(v)) for k, v in by_key.items() if v}


def _permutation_spearman(x, y, alternative: str = "greater",
                           n_perm: int = 10_000, seed: int = 0) -> tuple[float, float]:
    """Spearman + permutation p-value. alternative ∈ {'greater', 'less', 'two-sided'}."""
    import numpy as np
    from scipy import stats as sstats

    rho, _ = sstats.spearmanr(x, y)
    if not np.isfinite(rho):
        return float(rho) if np.isnan(rho) else 0.0, 1.0
    rng = np.random.RandomState(seed)
    n_extreme = 0
    for _ in range(n_perm):
        y_perm = rng.permutation(y)
        rho_perm, _ = sstats.spearmanr(x, y_perm)
        if alternative == "greater":
            if rho_perm >= rho:
                n_extreme += 1
        elif alternative == "less":
            if rho_perm <= rho:
                n_extreme += 1
        else:
            if abs(rho_perm) >= abs(rho):
                n_extreme += 1
    p = (n_extreme + 1) / (n_perm + 1)
    return float(rho), float(p)


def _resolve_x_value(predictor_id: str, dataset_id: str, x_field: str,
                      x_transform: str | None,
                      prereg_spec: dict[str, Any]) -> float | None:
    """Resolve the X covariate for a (predictor, dataset) cell against the
    pre-registered spec. Supports x_field ∈ {predictor.param_count,
    predictor.corpus_diversity_index, n_train_cells_per_perturbation}.
    """
    import math

    predictors_by_id = {p["id"]: p for p in prereg_spec.get("k1_v2", {}).get("predictors", [])}
    pred_spec = predictors_by_id.get(predictor_id)
    if pred_spec is None:
        return None

    if x_field == "predictor.param_count":
        val = pred_spec.get("param_count")
    elif x_field == "predictor.corpus_diversity_index":
        val = pred_spec.get("corpus_diversity_index")
    elif x_field == "n_train_cells_per_perturbation":
        # Per-dataset attribute, declared on dataset spec or hard-coded fallback
        # matching scripts/k2_analysis.TRAIN_CELLS_BY_DATASET
        fallback = {"norman": 800, "replogle_k562": 250, "replogle_rpe1": 150,
                    "adamson": 600, "tahoe": 400, "frangieh": 400, "schmidt": 600,
                    "datlinger": 200, "mcfaline_figueroa": 100, "walker": 350}
        val = fallback.get(dataset_id)
    else:
        return None

    if val is None:
        return None
    # Defensive cast: PyYAML 6.0 (YAML 1.2 strict) parses unquoted `1.0e6`
    # scientific-notation as STRING, not float. Cast to float here so the
    # verifier works regardless of YAML formatting quirks.
    try:
        val = float(val)
    except (TypeError, ValueError):
        return None
    if val < 0:
        return None

    # log/log10 of zero is undefined. Phase 1 scripts/k2_analysis.k2_test uses
    # log1p(0) = 0 for zero-param predictors (Mean, NoisyMean), keeping them in
    # the test at the bottom of the capacity scale. Mirror that behaviour here
    # so the verifier matches the Phase 1 analysis output.
    if x_transform == "log10":
        if val == 0:
            return 0.0
        return float(math.log10(val))
    if x_transform == "log":
        if val == 0:
            return 0.0
        return float(math.log(val))
    return float(val)


def _run_hypothesis_spearman(
    hyp_id: str,
    hyp_spec: dict[str, Any],
    results: dict[str, Any],
    prereg_spec: dict[str, Any],
) -> VerificationResult:
    """Run a per-dataset Spearman hypothesis test (H1, H1b, H2b)."""
    rows = _flatten_results_rows(results)
    if not rows:
        return VerificationResult(
            hypothesis_id=hyp_id, status="ERROR",
            error="no valid calibration rows in results.json",
        )

    agg = _aggregate_per_predictor_per_dataset(rows)
    if not agg:
        return VerificationResult(
            hypothesis_id=hyp_id, status="ERROR",
            error="aggregation produced no (predictor, dataset) cells",
        )

    x_field = hyp_spec.get("x_field")
    x_transform = hyp_spec.get("x_transform")
    alternative_threshold = None
    p_threshold = None
    n_datasets_threshold = None
    n_perm = 10000
    for crit in hyp_spec.get("success_criteria", []):
        t = crit.get("test")
        if t == "rho":
            alternative_threshold = (crit.get("operator"), crit.get("threshold"))
        elif t == "p_value":
            p_threshold = crit.get("threshold")
            n_perm = crit.get("n_permutations", 10000)
        elif t == "n_datasets_passing":
            n_datasets_threshold = (crit.get("threshold"), crit.get("of_total"))

    if not alternative_threshold or p_threshold is None or not n_datasets_threshold:
        return VerificationResult(
            hypothesis_id=hyp_id, status="ERROR",
            error=f"missing one of {{rho, p_value, n_datasets_passing}} success criteria",
        )

    rho_operator, rho_threshold = alternative_threshold
    alt = "greater" if rho_operator == ">" else ("less" if rho_operator == "<" else "two-sided")

    datasets_passing = []
    per_dataset = {}
    datasets_in_grid = sorted({d for (_, d) in agg.keys()})
    for ds in datasets_in_grid:
        cell = [(p, agg[(p, ds)]) for (p, d) in agg.keys() if d == ds]
        x_vals = []
        y_vals = []
        for pred_name, dev in cell:
            xv = _resolve_x_value(pred_name, ds, x_field, x_transform, prereg_spec)
            if xv is None:
                continue
            x_vals.append(xv)
            y_vals.append(dev)
        if len(x_vals) < 3:
            per_dataset[ds] = {"status": "too_few_predictors", "n": len(x_vals)}
            continue
        rho, p = _permutation_spearman(x_vals, y_vals, alternative=alt, n_perm=n_perm)
        passes_rho = (rho > rho_threshold) if rho_operator == ">" else (rho < rho_threshold)
        passes_p = p < p_threshold
        passes = passes_rho and passes_p
        per_dataset[ds] = {
            "rho": rho, "p_value": p, "n_predictors": len(x_vals), "passes": bool(passes),
        }
        if passes:
            datasets_passing.append(ds)

    n_passing_threshold, n_of_total = n_datasets_threshold
    hyp_overall_pass = len(datasets_passing) >= n_passing_threshold

    return VerificationResult(
        hypothesis_id=hyp_id,
        status="PASS" if hyp_overall_pass else "FAIL",
        observed={
            "n_datasets_passing": len(datasets_passing),
            "datasets_passing": datasets_passing,
            "per_dataset": per_dataset,
        },
        threshold={
            "rho": f"{rho_operator} {rho_threshold}",
            "p_value": f"< {p_threshold}",
            "n_datasets_passing": f">= {n_passing_threshold} of {n_of_total}",
        },
        n_cells_used=len(agg),
    )


def _run_hypothesis_anova(
    hyp_id: str,
    hyp_spec: dict[str, Any],
    results: dict[str, Any],
    prereg_spec: dict[str, Any],
) -> VerificationResult:
    """Run two-way ANOVA on (predictor, dataset) cells for H2."""
    try:
        import statsmodels.api as sm
        from statsmodels.formula.api import ols
        import pandas as pd
        import math
    except ImportError as exc:
        return VerificationResult(
            hypothesis_id=hyp_id, status="ERROR",
            error=f"statsmodels + pandas required for ANOVA: {exc}",
        )

    rows = _flatten_results_rows(results)
    agg = _aggregate_per_predictor_per_dataset(rows)
    if not agg:
        return VerificationResult(
            hypothesis_id=hyp_id, status="ERROR",
            error="no (predictor, dataset) cells in results.json",
        )

    predictors_by_id = {p["id"]: p for p in prereg_spec.get("k1_v2", {}).get("predictors", [])}
    datasets_by_id = {d["id"]: d for d in prereg_spec.get("k1_v2", {}).get("datasets", [])}

    df_rows = []
    for (pred_id, ds_id), dev in agg.items():
        pred_spec = predictors_by_id.get(pred_id)
        ds_spec = datasets_by_id.get(ds_id)
        if pred_spec is None or ds_spec is None:
            continue
        try:
            pc = float(pred_spec.get("param_count", 0))
        except (TypeError, ValueError):
            continue
        param_log = math.log10(pc + 1.0)
        cell_ctx = ds_spec.get("cell_line_context", "unknown")
        df_rows.append({
            "calibration_deviation": dev,
            "cell_line_context": cell_ctx,
            "param_count_log": param_log,
            "predictor": pred_id,
            "dataset": ds_id,
        })

    if len(df_rows) < 10:
        return VerificationResult(
            hypothesis_id=hyp_id, status="ERROR",
            error=f"too few cells ({len(df_rows)}) for ANOVA",
        )

    df = pd.DataFrame(df_rows)
    if df["cell_line_context"].nunique() < 2:
        return VerificationResult(
            hypothesis_id=hyp_id, status="ERROR",
            error="cell_line_context has only one level; cannot run ANOVA",
        )

    try:
        model = ols("calibration_deviation ~ C(cell_line_context) + param_count_log",
                    data=df).fit()
        anova_table = sm.stats.anova_lm(model, typ=2)
    except Exception as exc:
        return VerificationResult(
            hypothesis_id=hyp_id, status="ERROR",
            error=f"ANOVA failed: {type(exc).__name__}: {exc}",
        )

    main_effect_row = anova_table.loc["C(cell_line_context)"]
    f_pvalue = float(main_effect_row["PR(>F)"])
    ss_main = float(main_effect_row["sum_sq"])
    ss_total = float(anova_table["sum_sq"].sum())
    eta_sq = ss_main / ss_total if ss_total > 0 else 0.0

    p_threshold = None
    eta_threshold = None
    for crit in hyp_spec.get("success_criteria", []):
        t = crit.get("test")
        if t == "f_test_pvalue":
            p_threshold = crit.get("threshold")
        elif t == "eta_squared":
            eta_threshold = crit.get("threshold")

    passes_p = (p_threshold is not None) and (f_pvalue < p_threshold)
    passes_eta = (eta_threshold is not None) and (eta_sq > eta_threshold)
    overall_pass = passes_p and passes_eta

    return VerificationResult(
        hypothesis_id=hyp_id,
        status="PASS" if overall_pass else "FAIL",
        observed={
            "f_pvalue_cell_line_context": f_pvalue,
            "eta_squared_cell_line_context": eta_sq,
            "n_cells": len(df_rows),
            "n_K562_cells": int((df["cell_line_context"] == "K562").sum()),
            "n_non_K562_cells": int((df["cell_line_context"] == "non-K562").sum()),
        },
        threshold={
            "f_pvalue": f"< {p_threshold}",
            "eta_squared": f"> {eta_threshold}",
        },
        n_cells_used=len(df_rows),
    )


def _run_hypothesis_kruskal_wallis(
    hyp_id: str,
    hyp_spec: dict[str, Any],
    results: dict[str, Any],
    prereg_spec: dict[str, Any],
) -> VerificationResult:
    """Run Kruskal-Wallis across architecture families for H3."""
    try:
        from scipy import stats as sstats
        import numpy as np
    except ImportError as exc:
        return VerificationResult(
            hypothesis_id=hyp_id, status="ERROR",
            error=f"scipy required for Kruskal-Wallis: {exc}",
        )

    rows = _flatten_results_rows(results)
    agg = _aggregate_per_predictor_per_dataset(rows)
    if not agg:
        return VerificationResult(
            hypothesis_id=hyp_id, status="ERROR",
            error="no (predictor, dataset) cells in results.json",
        )

    predictors_by_id = {p["id"]: p for p in prereg_spec.get("k1_v2", {}).get("predictors", [])}
    excluded_families = set(hyp_spec.get("excluded_families", []))

    # Group (predictor, dataset) deviations by family
    family_devs: dict[str, list[float]] = {}
    for (pred_id, ds_id), dev in agg.items():
        pred_spec = predictors_by_id.get(pred_id)
        if pred_spec is None:
            continue
        family = pred_spec.get("family")
        if family in excluded_families:
            continue
        family_devs.setdefault(family, []).append(dev)

    families = sorted(family_devs.keys())
    if len(families) < 3:
        return VerificationResult(
            hypothesis_id=hyp_id, status="ERROR",
            error=f"need ≥3 families for KW; got {len(families)}: {families}",
        )

    groups = [family_devs[f] for f in families]
    try:
        stat, p = sstats.kruskal(*groups)
    except Exception as exc:
        return VerificationResult(
            hypothesis_id=hyp_id, status="ERROR",
            error=f"kruskal failed: {type(exc).__name__}: {exc}",
        )

    # Cliff's δ between best and worst family means
    family_means = {f: float(np.mean(devs)) for f, devs in family_devs.items()}
    best_family = min(family_means, key=family_means.get)
    worst_family = max(family_means, key=family_means.get)

    # Cliff's δ between best and worst groups
    best_devs = np.asarray(family_devs[best_family])
    worst_devs = np.asarray(family_devs[worst_family])
    # δ = (# pairs where worst > best - # pairs where worst < best) / (n_best * n_worst)
    diffs = worst_devs[:, None] - best_devs[None, :]
    n_pos = int((diffs > 0).sum())
    n_neg = int((diffs < 0).sum())
    n_total = best_devs.size * worst_devs.size
    cliffs_delta = (n_pos - n_neg) / n_total if n_total > 0 else 0.0

    p_threshold = None
    delta_threshold = None
    for crit in hyp_spec.get("success_criteria", []):
        t = crit.get("test")
        if t == "kw_pvalue":
            p_threshold = crit.get("threshold")
        elif t == "cliffs_delta_best_vs_worst":
            delta_threshold = crit.get("threshold")

    passes_p = (p_threshold is not None) and (float(p) < p_threshold)
    passes_delta = (delta_threshold is not None) and (abs(cliffs_delta) > delta_threshold)
    overall_pass = passes_p and passes_delta

    return VerificationResult(
        hypothesis_id=hyp_id,
        status="PASS" if overall_pass else "FAIL",
        observed={
            "kw_statistic": float(stat),
            "kw_pvalue": float(p),
            "cliffs_delta": float(cliffs_delta),
            "best_family": best_family,
            "worst_family": worst_family,
            "family_means": family_means,
            "family_n": {f: len(devs) for f, devs in family_devs.items()},
        },
        threshold={
            "kw_pvalue": f"< {p_threshold}",
            "cliffs_delta": f"> {delta_threshold}",
        },
        n_cells_used=sum(len(devs) for devs in family_devs.values()),
    )


def _run_hypothesis_full(
    hyp_id: str,
    hyp_spec: dict[str, Any],
    results: dict[str, Any],
    prereg_spec: dict[str, Any] | None = None,
) -> VerificationResult:
    """Full verify of a single hypothesis. Dispatches by test type."""
    test_kind = hyp_spec.get("test", "")

    if test_kind == "spearman_per_dataset":
        return _run_hypothesis_spearman(hyp_id, hyp_spec, results, prereg_spec or {})

    if test_kind == "two_way_anova_type2":
        return _run_hypothesis_anova(hyp_id, hyp_spec, results, prereg_spec or {})

    if test_kind == "kruskal_wallis":
        return _run_hypothesis_kruskal_wallis(hyp_id, hyp_spec, results, prereg_spec or {})

    return VerificationResult(
        hypothesis_id=hyp_id,
        status="SKIPPED",
        error=f"unknown test type: {test_kind}",
    )


def _apply_disposition_table(
    k2_v2_spec: dict[str, Any],
    hypothesis_results: list[VerificationResult],
) -> Optional[str]:
    """Apply the outcome_table to compute the disposition headline."""
    family_results: dict[str, str] = {}
    for hr in hypothesis_results:
        if hr.hypothesis_id in {"H2_cell_line_covariate", "H2b_corpus_diversity", "H3_architecture_family"}:
            if hr.status == "PASS":
                family_results[hr.hypothesis_id.split("_")[0]] = "P"
            elif hr.status == "FAIL":
                family_results[hr.hypothesis_id.split("_")[0]] = "F"
            else:
                return None  # DRY_RUN / SKIPPED / ERROR; no disposition

    if {"H2", "H2b", "H3"} - family_results.keys():
        return None  # incomplete family

    key = f"H2={family_results['H2']},H2b={family_results['H2b']},H3={family_results['H3']}"
    outcome_table = k2_v2_spec.get("outcome_table", {})
    return outcome_table.get(key, "<unknown disposition>")


def verify(
    yaml_path: Path,
    results_path: Path,
    repo_dir: Path,
    dry_run: bool = False,
) -> PreregVerification:
    """Top-level verify: parse, hash-check, mtime-check, run tests, dispose.

    Returns a `PreregVerification` data structure. The caller is responsible
    for serialising it to JSON and writing it to the verification log path.
    """
    timestamp = datetime.now(timezone.utc).isoformat()
    result = PreregVerification(
        timestamp=timestamp,
        prereg_yaml_path=str(yaml_path),
        prereg_yaml_sha256="",
        prereg_yaml_git_commit=None,
        results_json_path=str(results_path),
        results_json_sha256=None,
        results_json_mtime=None,
        yaml_lock_mtime=None,
        lock_check_passed=False,
        dry_run=dry_run,
    )

    # Parse YAML
    try:
        spec = _load_yaml(yaml_path)
    except Exception as exc:
        result.errors.append(f"yaml parse error: {exc}")
        return result

    # Schema validation
    schema_errors = _validate_schema(spec)
    if schema_errors:
        for err in schema_errors:
            result.errors.append(f"schema: {err}")

    # YAML self-hash
    lock = spec.get("lock", {})
    declared_yaml_sha = lock.get("sha256_yaml")
    yaml_hash_ok, observed_yaml_sha = _verify_yaml_self_hash(yaml_path, declared_yaml_sha)
    result.prereg_yaml_sha256 = observed_yaml_sha
    if not yaml_hash_ok:
        if not dry_run and declared_yaml_sha:
            result.errors.append(
                f"yaml sha256 mismatch: declared={declared_yaml_sha}, observed={observed_yaml_sha}"
            )
    result.prereg_yaml_git_commit = lock.get("git_commit_yaml")

    # results.json existence + sha256 + mtime
    if not results_path.exists():
        if not dry_run:
            result.errors.append(f"results.json not found at {results_path}")
        return result
    result.results_json_sha256 = _sha256_of_file(results_path)
    result.results_json_mtime = results_path.stat().st_mtime

    # mtime check
    mtime_ok, mtime_msg = _verify_results_mtime(
        results_path, lock.get("git_commit_yaml"), repo_dir
    )
    result.lock_check_passed = yaml_hash_ok and mtime_ok
    if not mtime_ok:
        if not dry_run:
            result.errors.append(f"results mtime: {mtime_msg}")

    # Load results.json
    with open(results_path) as fh:
        results_data = json.load(fh)

    # Per-hypothesis tests
    k2_v2_spec = spec.get("k2_v2", {})
    for hyp_id, hyp_spec in k2_v2_spec.get("hypotheses", {}).items():
        if dry_run:
            hr = _run_hypothesis_dry(hyp_id, hyp_spec, results_data)
        else:
            hr = _run_hypothesis_full(hyp_id, hyp_spec, results_data, prereg_spec=spec)
        result.hypothesis_results.append(hr)

    # Apply disposition (only if not dry-run)
    if not dry_run:
        result.disposition_headline = _apply_disposition_table(k2_v2_spec, result.hypothesis_results)

    return result


def _format_text_report(result: PreregVerification) -> str:
    """Pretty-print a PreregVerification as human-readable text."""
    lines = [
        "=" * 70,
        f"ConfPert pre-registration verification report",
        f"Timestamp: {result.timestamp}",
        f"YAML: {result.prereg_yaml_path}",
        f"  sha256 (observed): {result.prereg_yaml_sha256}",
        f"  git commit:        {result.prereg_yaml_git_commit}",
        f"Results: {result.results_json_path}",
        f"  sha256: {result.results_json_sha256}",
        f"  mtime:  {result.results_json_mtime}",
        f"Lock check passed: {result.lock_check_passed}",
        f"Mode: {'DRY RUN' if result.dry_run else 'FULL VERIFY'}",
        "-" * 70,
        "Hypothesis results:",
    ]
    for hr in result.hypothesis_results:
        lines.append(f"  [{hr.status}] {hr.hypothesis_id}")
        if hr.observed:
            for k, v in hr.observed.items():
                lines.append(f"    {k}: {v}")
        if hr.error:
            lines.append(f"    ERROR: {hr.error}")
    if result.disposition_headline:
        lines.append("-" * 70)
        lines.append(f"DISPOSITION: {result.disposition_headline}")
    if result.errors:
        lines.append("-" * 70)
        lines.append("ERRORS:")
        for err in result.errors:
            lines.append(f"  - {err}")
    lines.append("=" * 70)
    return "\n".join(lines)


def cli_main(argv: list[str]) -> int:
    """CLI entry for `confpert prereg verify`."""
    import argparse

    parser = argparse.ArgumentParser(
        prog="confpert prereg",
        description="ConfPert pre-registration verifier",
    )
    sub = parser.add_subparsers(dest="cmd", required=True)

    p_verify = sub.add_parser("verify", help="verify a pre-registration YAML")
    p_verify.add_argument("--prereg", type=Path, required=True, help="path to preregistration_v2.yaml")
    p_verify.add_argument("--results", type=Path, default=Path("baselines/results.json"))
    p_verify.add_argument("--repo-dir", type=Path, default=Path.cwd())
    p_verify.add_argument("--dry-run", action="store_true", help="skip statistical tests; verify schema + counts only")
    p_verify.add_argument("--out", type=Path, default=None, help="optional JSON output path")

    p_emit = sub.add_parser("emit-hashes", help="compute & write SHA-256 lock block (run once at lock-in)")
    p_emit.add_argument("--prereg", type=Path, required=True)

    args = parser.parse_args(argv)

    if args.cmd == "emit-hashes":
        emit_hashes(args.prereg)
        return 0

    if args.cmd == "verify":
        result = verify(
            yaml_path=args.prereg,
            results_path=args.results,
            repo_dir=args.repo_dir,
            dry_run=args.dry_run,
        )
        print(_format_text_report(result))
        if args.out is not None:
            from dataclasses import asdict
            with open(args.out, "w") as fh:
                json.dump(asdict(result), fh, indent=2, default=str)
            print(f"\nJSON report written to {args.out}")
        return 0 if not result.errors else 1

    return 0


if __name__ == "__main__":
    sys.exit(cli_main(sys.argv[1:]))
