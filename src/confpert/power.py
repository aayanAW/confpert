"""Power analysis for K2 v2 hypothesis family (H2 ANOVA, H2b Spearman, H3 Kruskal-Wallis).

Implements minimum-detectable-effect computations and post-hoc power checks
referenced in `preregistration_v2.md` §2.7 and `PHASE_2_PLAN.md` §14.4.

Each function returns a `PowerResult` dataclass with:
  - test name
  - sample size used
  - threshold / effect size targeted
  - estimated power (probability of detecting the effect at α)
  - minimum detectable effect at 80% power
  - notes (e.g. underpowering warning, computation method)

Two approaches:
  - Analytical (via `statsmodels.stats.power`) where available (F-test / ANOVA)
  - Monte-Carlo simulation otherwise (Spearman, Kruskal-Wallis)
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional


@dataclass
class PowerResult:
    test_name: str
    n: int
    alpha: float
    target_effect: float
    estimated_power: float
    min_detectable_effect_at_80: float
    notes: str = ""


def power_two_way_anova(
    n_cells: int,
    n_factors: int = 2,
    target_eta_sq: float = 0.10,
    alpha: float = 0.0167,
) -> PowerResult:
    """Power of detecting a main-effect η² in two-way ANOVA via F-test.

    Uses statsmodels.stats.power.FTestPower. Cohen's f² = η² / (1 - η²).
    """
    try:
        from statsmodels.stats.power import FTestAnovaPower
    except ImportError:
        return PowerResult(
            test_name="two_way_anova",
            n=n_cells,
            alpha=alpha,
            target_effect=target_eta_sq,
            estimated_power=float("nan"),
            min_detectable_effect_at_80=float("nan"),
            notes="statsmodels not installed; install via `pip install statsmodels`.",
        )

    # Cohen's f for effect size (used by FTestAnovaPower in `effect_size` argument)
    # f = sqrt(eta_sq / (1 - eta_sq))
    f = (target_eta_sq / (1 - target_eta_sq)) ** 0.5
    power = FTestAnovaPower().power(
        effect_size=f, nobs=n_cells, alpha=alpha, k_groups=n_factors + 1
    )

    # Find min detectable effect at 80% power: binary search on eta_sq
    lo, hi = 1e-4, 0.5
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        f_mid = (mid / (1 - mid)) ** 0.5
        p_mid = FTestAnovaPower().power(
            effect_size=f_mid, nobs=n_cells, alpha=alpha, k_groups=n_factors + 1
        )
        if p_mid < 0.80:
            lo = mid
        else:
            hi = mid
    mde = 0.5 * (lo + hi)

    return PowerResult(
        test_name="two_way_anova",
        n=n_cells,
        alpha=alpha,
        target_effect=target_eta_sq,
        estimated_power=float(power),
        min_detectable_effect_at_80=float(mde),
        notes=f"FTestAnovaPower with k_groups={n_factors+1}",
    )


def power_spearman_per_dataset(
    n_predictors: int,
    target_rho: float = 0.5,
    alpha: float = 0.0167,
    n_simulations: int = 2000,
    seed: int = 42,
) -> PowerResult:
    """Monte-Carlo estimate of Spearman-correlation test power.

    Simulates n_simulations draws of n pairs (X, Y) where Y is correlated with X
    at the target rho via a Gaussian copula. Counts fraction with Spearman p < alpha.
    """
    import numpy as np
    from scipy import stats as sstats

    rng = np.random.default_rng(seed)
    target_abs = abs(target_rho)

    n_hits = 0
    for _ in range(n_simulations):
        x = rng.standard_normal(n_predictors)
        eps = rng.standard_normal(n_predictors)
        # Generate y with target Pearson correlation; Spearman ~= Pearson for
        # joint-Gaussian data
        y = target_abs * x + ((1 - target_abs**2) ** 0.5) * eps
        # Use absolute observed rho — pre-reg uses one-sided test, but absolute
        # detection captures both H2b (negative) and H1 (positive) symmetrically
        rho, p = sstats.spearmanr(x, y)
        if p < alpha:
            n_hits += 1
    power = n_hits / n_simulations

    # MDE: find rho that gives 80% power at this n
    # Bisection on positive rho only (symmetric for spearman)
    lo, hi = 0.01, 0.99
    for _ in range(15):  # 2^15 ≈ 32768 simulations total — keep budget reasonable
        mid = 0.5 * (lo + hi)
        n_hits_mid = 0
        for _ in range(max(200, n_simulations // 10)):
            x = rng.standard_normal(n_predictors)
            eps = rng.standard_normal(n_predictors)
            y = mid * x + ((1 - mid**2) ** 0.5) * eps
            _, p = sstats.spearmanr(x, y)
            if p < alpha:
                n_hits_mid += 1
        p_mid = n_hits_mid / max(200, n_simulations // 10)
        if p_mid < 0.80:
            lo = mid
        else:
            hi = mid
    mde = 0.5 * (lo + hi)

    notes = f"Monte-Carlo n_simulations={n_simulations}, seed={seed}"
    if power < 0.80:
        notes += f"; UNDERPOWERED: estimated power {power:.2f} < 0.80 at target rho {target_rho}"
    return PowerResult(
        test_name="spearman",
        n=n_predictors,
        alpha=alpha,
        target_effect=target_rho,
        estimated_power=float(power),
        min_detectable_effect_at_80=float(mde),
        notes=notes,
    )


def power_kruskal_wallis(
    n_per_group: int,
    k_groups: int,
    target_cliffs_delta: float = 0.30,
    alpha: float = 0.0167,
    n_simulations: int = 2000,
    seed: int = 42,
) -> PowerResult:
    """Monte-Carlo estimate of Kruskal-Wallis test power.

    Simulates k_groups groups of n_per_group samples each. Group means
    separated by an offset chosen so Cliff's δ between extreme groups equals
    `target_cliffs_delta`. (Approximation: δ ≈ Φ(d/√2) - Φ(-d/√2) where d is
    Cohen's d for standard-normal groups.)
    """
    import numpy as np
    from scipy import stats as sstats

    rng = np.random.default_rng(seed)

    # Approximate offset for target Cliff's δ between extreme groups
    # δ ≈ 2 * Φ(d/√2) - 1 where d is Cohen's d
    # → d = √2 * Φ⁻¹((δ+1)/2)
    from scipy.stats import norm
    cohen_d = 2**0.5 * norm.ppf((target_cliffs_delta + 1) / 2)

    n_hits = 0
    for _ in range(n_simulations):
        groups = []
        for i in range(k_groups):
            shift = cohen_d * i / (k_groups - 1) if k_groups > 1 else 0.0
            groups.append(rng.standard_normal(n_per_group) + shift)
        stat, p = sstats.kruskal(*groups)
        if p < alpha:
            n_hits += 1
    power = n_hits / n_simulations

    # MDE: bisection on cliffs_delta
    lo, hi = 0.05, 0.95
    for _ in range(15):
        mid = 0.5 * (lo + hi)
        cohen_d_mid = 2**0.5 * norm.ppf((mid + 1) / 2)
        n_hits_mid = 0
        for _ in range(max(200, n_simulations // 10)):
            groups = []
            for i in range(k_groups):
                shift = cohen_d_mid * i / (k_groups - 1) if k_groups > 1 else 0.0
                groups.append(rng.standard_normal(n_per_group) + shift)
            _, p = sstats.kruskal(*groups)
            if p < alpha:
                n_hits_mid += 1
        p_mid = n_hits_mid / max(200, n_simulations // 10)
        if p_mid < 0.80:
            lo = mid
        else:
            hi = mid
    mde = 0.5 * (lo + hi)

    notes = f"Monte-Carlo n_simulations={n_simulations}, k_groups={k_groups}, seed={seed}"
    if power < 0.80:
        notes += f"; UNDERPOWERED"
    return PowerResult(
        test_name="kruskal_wallis",
        n=n_per_group * k_groups,
        alpha=alpha,
        target_effect=target_cliffs_delta,
        estimated_power=float(power),
        min_detectable_effect_at_80=float(mde),
        notes=notes,
    )


def k2_v2_power_report() -> dict:
    """Compute the full K2 v2 power-analysis table per `preregistration_v2.md` §2.7.

    Returns a dict suitable for JSON-dumping into the verification report.
    """
    return {
        "H2_anova": power_two_way_anova(
            n_cells=120, n_factors=2, target_eta_sq=0.10, alpha=0.0167
        ),
        # H2b: n=14 predictors per dataset, ρ < -0.5
        "H2b_spearman_per_dataset": power_spearman_per_dataset(
            n_predictors=14, target_rho=-0.5, alpha=0.0167, n_simulations=2000
        ),
        # H3: pre-reg v2.0 final lock uses per-(predictor, dataset) cells as
        # sample unit (dependent_aggregation = mean_over_alphas_and_discrepancies).
        # F1 baselines (6) + F2 latent-delta (4) + F3 transformer (5) = 15 predictors
        # x 10 datasets = 150 cells. Per-family n: F1=60, F2=40, F3=50. Worst-case
        # group n ~= 40 for KW.
        "H3_kruskal_wallis": power_kruskal_wallis(
            n_per_group=40, k_groups=3, target_cliffs_delta=0.50, alpha=0.0167,
            n_simulations=2000,
        ),
    }


def k1_v2_h1_h1b_power_report() -> dict:
    """Power analysis for H1 (capacity) and H1b (data) tests re-run on the
    expanded Phase 2 grid (14 predictors x 10 datasets).

    H1/H1b carry over from Phase 1 with the same alpha=0.05 threshold and the
    same rho>0.5/<-0.5 thresholds. They are NOT part of the K2 v2 family-wise
    correction (which only covers H2/H2b/H3). Sample sizes at Phase 2:

      H1 per-dataset Spearman: n=14 predictors (or 15 with NoisyMean
        control), at alpha=0.05. Bigger than Phase 1's n=8-9, so power should
        be substantially better.

      H1b per-dataset Spearman: same n=14, but X is constant per dataset
        (training-cell-count doesn't vary across predictors WITHIN a dataset),
        so the per-dataset test will produce ρ=NaN. The meaningful H1b test is
        cross-dataset: one mean-deviation per dataset (n=10), one X per dataset
        (log10 training cells). Power at n=10 for ρ<-0.5 is computed too.
    """
    return {
        "H1_capacity_per_dataset_phase2": power_spearman_per_dataset(
            n_predictors=14, target_rho=0.5, alpha=0.05, n_simulations=2000,
        ),
        "H1b_data_cross_dataset_phase2": power_spearman_per_dataset(
            n_predictors=10, target_rho=-0.5, alpha=0.05, n_simulations=2000,
        ),
    }


def full_power_report() -> dict:
    """Combined K1+K2 v2 power report for Phase 2A.6 / Phase 2D analyses."""
    out = {}
    out.update(k1_v2_h1_h1b_power_report())
    out.update(k2_v2_power_report())
    return out


if __name__ == "__main__":
    import json
    from dataclasses import asdict
    report = k2_v2_power_report()
    print(json.dumps({k: asdict(v) for k, v in report.items()}, indent=2))
