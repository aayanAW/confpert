"""ConfPert top-level CLI dispatcher.

Subcommands:
  prereg verify        -- verify a pre-registration YAML against results.json
  prereg emit-hashes   -- compute & write SHA-256 lock block (once at lock-in)
  prereg-ablate        -- run pre-registration ablation panel (Pipeline A vs B)
  phase2d              -- full Phase 2D pipeline: verify + ablate + power + bootstrap
  style-audit          -- LLM-style audit on a .tex / .md file (em-dash, bold, etc.)
"""
from __future__ import annotations

import sys


def main() -> int:
    if len(sys.argv) < 2 or sys.argv[1] in {"-h", "--help"}:
        print(
            "confpert — Calibrated distributional uncertainty for single-cell "
            "perturbation predictors.\n\n"
            "Usage:\n"
            "  confpert prereg verify --prereg <yaml> [--results <json>] [--dry-run] [--out <json>]\n"
            "  confpert prereg emit-hashes --prereg <yaml>\n"
            "  confpert prereg-ablate --prereg <yaml> --results <json> [--out <json>]\n"
            "  confpert phase2d --prereg <yaml> --results <json> [--out <json>] [--skip-ablation]\n"
            "  confpert style-audit <path.tex|path.md> [--strict]\n"
            "\n"
            "See confpert/PHASE_2_PLAN.md and paper_neurips_dnb_2026/preregistration_v2.md.\n"
        )
        return 0

    sub = sys.argv[1]
    rest = sys.argv[2:]

    if sub == "prereg":
        from confpert import prereg
        return prereg.cli_main(rest)

    if sub == "prereg-ablate":
        from confpert import prereg_ablation
        return prereg_ablation.cli_main(rest)

    if sub == "phase2d":
        from confpert import phase2d_runner
        return phase2d_runner.cli_main(rest)

    if sub == "style-audit":
        from confpert import style_audit
        return style_audit.cli_main(rest)

    print(f"confpert: unknown subcommand '{sub}'. Run `confpert --help` for usage.")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
