"""K3 top-K signature-size robustness figure (2-panel PDF).

Reads `results/k3_tahoe_predictor_derived[_top{K}].json` for K in {25,50,100,200}
and emits a 2-panel chart showing (left) ConfPert vs uncalibrated signature counts
per K, (right) the ConfPert/uncalibrated ratio per K.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"
OUT = ROOT / "paper" / "fig_k3_robustness.pdf"


def _load(K: int) -> dict:
    name = "k3_tahoe_predictor_derived.json" if K == 50 else f"k3_tahoe_predictor_derived_top{K}.json"
    with open(RESULTS / name) as f:
        d = json.load(f)
    sk = d["k3_success_check"]
    return {
        "K": K,
        "n_confpert": sk["n_confpert_signatures"],
        "n_uncal_misses": sk["n_uncalibrated_misses_on_confpert_signatures"],
        "passed": bool(sk["success"]),
    }


def main() -> None:
    Ks = [25, 50, 100, 200]
    rows = [_load(K) for K in Ks]
    n_cf = [r["n_confpert"] for r in rows]
    n_un = [r["n_uncal_misses"] for r in rows]
    ratio = [c / u for c, u in zip(n_cf, n_un)]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 3.2))

    x = list(range(len(Ks)))
    width = 0.36
    ax1.bar([xi - width / 2 for xi in x], n_cf, width, label="ConfPert sigs",
            color="#2D6CDF")
    ax1.bar([xi + width / 2 for xi in x], n_un, width, label="uncal. misses",
            color="#C44E52")
    for xi, c, u in zip(x, n_cf, n_un):
        ax1.text(xi - width / 2, c + 15, str(c), ha="center", fontsize=8)
        ax1.text(xi + width / 2, u + 15, str(u), ha="center", fontsize=8)
    ax1.set_xticks(x)
    ax1.set_xticklabels([f"top-{K}" for K in Ks])
    ax1.set_ylabel("BH-FDR signatures at q < 0.05")
    ax1.set_title("(a) Per-K signature counts")
    ax1.legend(loc="upper left", fontsize=9, frameon=False)
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)

    ax2.plot(Ks, ratio, "o-", color="#2D6CDF", linewidth=2, markersize=8,
             label="ConfPert / uncalibrated")
    ax2.axhline(1.0, color="gray", linestyle="--", linewidth=1, label="parity")
    for K, r in zip(Ks, ratio):
        ax2.annotate(f"{r:.2f}×", (K, r), textcoords="offset points",
                     xytext=(6, 6), fontsize=9)
    ax2.set_xscale("log")
    ax2.set_xticks(Ks)
    ax2.set_xticklabels([str(K) for K in Ks])
    ax2.set_xlabel("top-K signature size (genes)")
    ax2.set_ylabel("ratio")
    ax2.set_title("(b) ConfPert / uncalibrated ratio")
    ax2.set_ylim(0.95, max(ratio) + 0.2)
    ax2.legend(loc="upper left", fontsize=9, frameon=False)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)

    plt.tight_layout()
    fig.savefig(OUT, bbox_inches="tight")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
