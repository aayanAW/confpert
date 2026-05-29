"""Tests for the LLM-style audit module."""
from __future__ import annotations

import pytest

from confpert.style_audit import audit_text, RULES


def test_em_dash_budget_flagged_over_budget():
    text = "First — second — third — fourth — fifth." * 50
    report = audit_text(text, chars_per_page=200)
    em_dash_finding = next(f for f in report.findings
                            if f.rule_id == "em_dash_budget")
    assert em_dash_finding.count > 0
    assert em_dash_finding.over_budget is True


def test_rule_of_three_caught():
    text = (
        "We find (i) calibration improves; (ii) capacity does not predict "
        "calibration; (iii) cell-line context dominates the picture."
    )
    report = audit_text(text)
    f = next(f for f in report.findings if f.rule_id == "rule_of_three_list")
    assert f.count >= 1


def test_bold_shouting_blocked():
    text = r"The result is \textbf{PASS} on every dataset."
    report = audit_text(text)
    f = next(f for f in report.findings if f.rule_id == "bold_shouting")
    assert f.count == 1
    assert f.over_budget is True
    assert report.n_block_violations >= 1


def test_self_praising_adjective_caught():
    text = "Our novel rigorous comprehensive framework " * 30
    report = audit_text(text, chars_per_page=100)
    f = next(f for f in report.findings if f.rule_id == "self_praising_adjective")
    assert f.count > 0


def test_we_introduce_repetition():
    text = "We introduce ConfPert. We introduce a benchmark. We introduce a verifier."
    report = audit_text(text)
    f = next(f for f in report.findings if f.rule_id == "we_introduce_repetition")
    assert f.count == 3
    assert f.over_budget is True


def test_strictly_better_blocked():
    text = "Method A is strictly better than method B."
    report = audit_text(text)
    f = next(f for f in report.findings
             if f.rule_id == "inflated_phrase_strictly_better")
    assert f.count == 1
    assert f.over_budget is True


def test_clean_prose_no_violations():
    text = (
        "ConfPert reports calibration deviation per discrepancy. Each head "
        "calibrates against a held-out perturbation set. The bilinear "
        "baseline matches transformer performance on five of six metrics."
    )
    report = audit_text(text)
    block_findings = [f for f in report.findings
                      if f.severity == "block" and f.over_budget]
    assert len(block_findings) == 0


def test_all_rules_have_fix_hint():
    for r in RULES:
        assert r.fix_hint.strip(), f"rule {r.id} missing fix_hint"


def test_load_bearing_under_budget_for_2_uses():
    text = "The cell-line covariate is load-bearing for H2 and load-bearing for K2."
    report = audit_text(text)
    f = next(f for f in report.findings
             if f.rule_id == "inflated_phrase_load_bearing")
    assert f.count == 2
    assert f.over_budget is False  # budget = 2


def test_pleasantry_blocked():
    text = "Of course, the result is naturally certainly correct."
    report = audit_text(text)
    f = next(f for f in report.findings if f.rule_id == "pleasantry_in_text")
    assert f.count >= 2
    assert f.over_budget is True
