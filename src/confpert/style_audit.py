r"""LLM-style audit on .tex / .md files per PHASE_2_PLAN.md §6 writing-repair rules.

Enforces concrete writing-style rules (no judgment calls — every rule is a
counted regex match):

  - Em-dash budget: ≤ 2 per estimated-page
  - Rule-of-three lists: occurrences of "(i)... (ii)... (iii)..." or
    "three findings: ..." etc.
  - Bold-shouting in prose: `\textbf{PASS}` / `\textbf{FAIL}` / `\textbf{EXACT}`
    outside table-cell context
  - Paragraph headers in short sections: `\paragraph{...}` count per page
  - Self-praising adjectives: "rigorous", "principled", "comprehensive",
    "extensive" used unquantified
  - "we introduce" repetition: >1 per abstract
  - Inflated phrases: load-bearing, strictly better, explicitly not pre-reg,
    meets-in-our-framework, etc.

Output is a JSON report + a per-finding diff-style commentary suitable for
direct review.

Usage:
  confpert style-audit paper_neurips_dnb_2026/tex/main.tex
  confpert style-audit paper_neurips_dnb_2026/preregistration_v2.md --strict
"""
from __future__ import annotations

import json
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any


# Approximate page-length heuristic: NeurIPS D&B template runs about 600
# words per main-text page, ~3500 characters per page including LaTeX
# noise. Override via --chars-per-page CLI flag.
CHARS_PER_PAGE = 3500


# ----------------------------------------------------------------------------
# Rules. Each rule is (id, regex, severity, budget_per_page, fix_hint).
# severity: "block" (counts toward strict-mode exit code) | "warn"
# budget_per_page: max occurrences per estimated page; None = absolute count
# ----------------------------------------------------------------------------

@dataclass
class StyleRule:
    id: str
    regex: re.Pattern
    severity: str
    budget_per_page: float | None
    absolute_max: int | None
    fix_hint: str
    # If True, drop matches that occur inside a LaTeX comment line (% ...)
    # or inside \begin{tabular}…\end{tabular} or \begin{table}…\end{table}.
    skip_in_table_or_comment: bool = False


RULES: list[StyleRule] = [
    StyleRule(
        id="em_dash_budget",
        # Match en/em dashes (—, --); LaTeX --- in prose renders as em-dash.
        # Exclude Markdown horizontal-rule lines (lines consisting of 3+ dashes
        # surrounded by whitespace, optionally end-of-string).
        regex=re.compile(r"(?<!\n)(?<![\-])(?:---|—|–)(?![\-])(?!\n)"),
        severity="block",
        budget_per_page=2.0,
        absolute_max=None,
        fix_hint="Replace excess em-dashes with period + new sentence OR comma.",
    ),
    StyleRule(
        id="rule_of_three_list",
        regex=re.compile(r"\(i\).{1,200}\(ii\).{1,200}\(iii\)", re.DOTALL),
        severity="block",
        budget_per_page=None,
        absolute_max=2,
        fix_hint="Rewrite (i)/(ii)/(iii) as three separate sentences in prose.",
    ),
    StyleRule(
        id="bold_shouting",
        # PASS/FAIL/EXACT/TRUE/FALSE wrapped in \textbf{} in prose. Matches
        # inside \begin{tabular}…\end{tabular}, \begin{table}…\end{table}, or
        # LaTeX comment lines are filtered (see skip_in_table_or_comment).
        regex=re.compile(r"\\textbf\{(PASS|FAIL|EXACT|TRUE|FALSE|YES|NO)\}"),
        severity="block",
        budget_per_page=None,
        absolute_max=0,
        fix_hint=(
            "Reserve bold for table cells only; in prose, replace with "
            "'passes', 'fails', etc."
        ),
        skip_in_table_or_comment=True,
    ),
    StyleRule(
        id="paragraph_header_in_short_section",
        # \paragraph{...} is a strong style smell when overused.
        regex=re.compile(r"\\paragraph\{[^}]+\}"),
        severity="warn",
        budget_per_page=1.0,
        absolute_max=None,
        fix_hint=(
            "Drop \\paragraph{} in any section under 1/2 page; inline as "
            "regular paragraph break."
        ),
    ),
    StyleRule(
        id="self_praising_adjective",
        regex=re.compile(
            r"\b(rigorous|principled|comprehensive|extensive|robust|seamless|"
            r"holistic|elegant|state-of-the-art|cutting-edge|novel|innovative)\b",
            re.IGNORECASE,
        ),
        severity="warn",
        budget_per_page=1.0,
        absolute_max=None,
        fix_hint=(
            "Remove unless quantified. 'rigorous' -> drop; 'comprehensive' -> "
            "drop; 'novel' -> show, don't tell."
        ),
    ),
    StyleRule(
        id="we_introduce_repetition",
        regex=re.compile(r"\bwe introduce\b", re.IGNORECASE),
        severity="warn",
        budget_per_page=None,
        absolute_max=1,
        fix_hint=(
            "Use 'we introduce' ONCE in the abstract. Elsewhere use "
            "'ConfPert reports', 'Each head calibrates', etc."
        ),
    ),
    StyleRule(
        id="inflated_phrase_load_bearing",
        regex=re.compile(r"\bload[-\s]bearing\b", re.IGNORECASE),
        severity="warn",
        budget_per_page=None,
        absolute_max=2,
        fix_hint="Replace 'load-bearing' with 'the main' or 'the strongest'.",
    ),
    StyleRule(
        id="inflated_phrase_strictly_better",
        regex=re.compile(r"\bstrictly better\b", re.IGNORECASE),
        severity="warn",
        budget_per_page=None,
        absolute_max=0,
        fix_hint="Replace 'strictly better' with 'lower' or 'smaller deviation'.",
    ),
    StyleRule(
        id="inflated_phrase_meets_in_our_framework",
        regex=re.compile(r"\bmeets\b.{1,40}\bin our framework\b", re.IGNORECASE),
        severity="warn",
        budget_per_page=None,
        absolute_max=0,
        fix_hint="Drop 'X meets Y in our framework' entirely.",
    ),
    StyleRule(
        id="filler_just_really_basically_actually",
        regex=re.compile(r"\b(just|really|basically|actually|simply)\b",
                         re.IGNORECASE),
        severity="warn",
        budget_per_page=2.0,
        absolute_max=None,
        fix_hint="Drop filler words. They add no information.",
    ),
    StyleRule(
        id="pleasantry_in_text",
        regex=re.compile(
            r"\b(of course|happy to|certainly|naturally)\b", re.IGNORECASE),
        severity="warn",
        budget_per_page=None,
        absolute_max=0,
        fix_hint="Drop pleasantries.",
    ),
]


@dataclass
class StyleFinding:
    rule_id: str
    severity: str
    count: int
    budget: float | int | None
    over_budget: bool
    examples: list[str] = field(default_factory=list)
    fix_hint: str = ""


@dataclass
class StyleReport:
    file: str
    n_chars: int
    n_estimated_pages: float
    findings: list[StyleFinding]
    over_budget_count: int
    n_block_violations: int


_TABLE_ENV_RE = re.compile(
    r"\\begin\{(?:tabular|table|tabularx|longtable)\b[^}]*\}.*?"
    r"\\end\{(?:tabular|table|tabularx|longtable)\}",
    re.DOTALL,
)


def _table_or_comment_ranges(text: str) -> list[tuple[int, int]]:
    """Return char-index ranges (start, end) inside which matches should be
    ignored for rules with skip_in_table_or_comment=True. Covers
    \begin{tabular}…\end{tabular}, \begin{table}…\end{table}, and any LaTeX
    line whose first non-whitespace character is '%' (comment).
    """
    ranges: list[tuple[int, int]] = []
    for m in _TABLE_ENV_RE.finditer(text):
        ranges.append((m.start(), m.end()))
    # Comment lines: scan line-by-line
    offset = 0
    for line in text.split("\n"):
        ls = line.lstrip()
        if ls.startswith("%"):
            ranges.append((offset, offset + len(line)))
        offset += len(line) + 1  # +1 for the newline
    return ranges


def _in_any_range(pos: int, ranges: list[tuple[int, int]]) -> bool:
    for s, e in ranges:
        if s <= pos < e:
            return True
    return False


def audit_text(text: str, chars_per_page: int = CHARS_PER_PAGE,
                file_label: str = "<text>") -> StyleReport:
    n_chars = len(text)
    n_pages = max(1.0, n_chars / chars_per_page)

    findings: list[StyleFinding] = []
    block_count = 0
    over_budget = 0
    skip_ranges = _table_or_comment_ranges(text)
    for rule in RULES:
        matches = list(rule.regex.finditer(text))
        if rule.skip_in_table_or_comment and skip_ranges:
            matches = [m for m in matches
                        if not _in_any_range(m.start(), skip_ranges)]
        n = len(matches)
        examples = []
        for m in matches[:3]:
            start = max(0, m.start() - 30)
            end = min(len(text), m.end() + 30)
            excerpt = text[start:end].replace("\n", "  ")
            examples.append(excerpt)

        over = False
        if rule.absolute_max is not None and n > rule.absolute_max:
            over = True
        if rule.budget_per_page is not None and n > rule.budget_per_page * n_pages:
            over = True

        if over:
            over_budget += 1
            if rule.severity == "block":
                block_count += 1

        findings.append(StyleFinding(
            rule_id=rule.id,
            severity=rule.severity,
            count=n,
            budget=(rule.absolute_max if rule.absolute_max is not None
                    else (rule.budget_per_page * n_pages
                          if rule.budget_per_page is not None else None)),
            over_budget=over,
            examples=examples,
            fix_hint=rule.fix_hint,
        ))

    return StyleReport(
        file=file_label,
        n_chars=n_chars,
        n_estimated_pages=n_pages,
        findings=findings,
        over_budget_count=over_budget,
        n_block_violations=block_count,
    )


def cli_main(argv: list[str]) -> int:
    import argparse
    p = argparse.ArgumentParser(description="ConfPert LLM-style audit")
    p.add_argument("path", help=".tex / .md / .txt file to audit")
    p.add_argument("--strict", action="store_true",
                   help="Exit non-zero if any 'block'-severity rule is over budget")
    p.add_argument("--json", action="store_true",
                   help="Emit JSON instead of human-readable text")
    p.add_argument("--chars-per-page", type=int, default=CHARS_PER_PAGE)
    args = p.parse_args(argv)

    path = Path(args.path)
    if not path.exists():
        print(f"style-audit: file not found: {path}")
        return 2

    text = path.read_text()
    report = audit_text(text, chars_per_page=args.chars_per_page,
                         file_label=str(path))

    if args.json:
        out = {
            "file": report.file,
            "n_chars": report.n_chars,
            "n_estimated_pages": report.n_estimated_pages,
            "over_budget_count": report.over_budget_count,
            "n_block_violations": report.n_block_violations,
            "findings": [
                {
                    "rule_id": f.rule_id, "severity": f.severity,
                    "count": f.count, "budget": f.budget,
                    "over_budget": f.over_budget, "fix_hint": f.fix_hint,
                    "examples": f.examples,
                }
                for f in report.findings
            ],
        }
        print(json.dumps(out, indent=2))
    else:
        print("=" * 72)
        print(f"LLM-style audit: {report.file}")
        print(f"  {report.n_chars} chars  ~{report.n_estimated_pages:.1f} pages")
        print(f"  Block violations: {report.n_block_violations}")
        print(f"  Over-budget rules: {report.over_budget_count}")
        print("=" * 72)
        for f in report.findings:
            marker = "[FAIL]" if f.over_budget else "[ ok ]"
            sev = f"({f.severity})"
            budget = f.budget if f.budget is not None else "∞"
            print(f"{marker} {sev:7} {f.rule_id:42} {f.count:4d} (budget ~{budget})")
            if f.over_budget:
                print(f"        -> {f.fix_hint}")
                for ex in f.examples:
                    print(f"        ex: ...{ex}...")
        print("=" * 72)

    if args.strict and report.n_block_violations > 0:
        return 1
    return 0


if __name__ == "__main__":
    import sys
    raise SystemExit(cli_main(sys.argv[1:]))
