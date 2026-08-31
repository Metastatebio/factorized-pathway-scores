"""Carrier feasibility and power approximations for interaction studies."""

from __future__ import annotations

from collections.abc import Iterable, Sequence

import numpy as np
import pandas as pd
from scipy.stats import norm


def joint_signature_prevalence(component_frequencies: Sequence[float]) -> float:
    """Return joint prevalence under an explicit independence assumption."""
    frequencies = np.asarray(component_frequencies, dtype=float)
    if frequencies.size == 0:
        raise ValueError("At least one component frequency is required.")
    if np.any((frequencies < 0) | (frequencies > 1)):
        raise ValueError("Component frequencies must lie between 0 and 1.")
    return float(np.prod(frequencies))


def interaction_power_normal(
    sample_size: int,
    signature_prevalence: float,
    interaction_effect_sd: float,
    residual_sd: float,
    alpha: float,
) -> float:
    """Approximate two-sided power for A×G under standardized independent A and binary G.

    Residualizing A×G against A and G gives approximate variance p(1-p), so the
    interaction z statistic has noncentrality beta*sqrt(n*p*(1-p))/sigma.
    """
    if sample_size <= 0:
        raise ValueError("sample_size must be positive.")
    if not 0 < signature_prevalence < 1:
        raise ValueError("signature_prevalence must lie strictly between 0 and 1.")
    if residual_sd <= 0:
        raise ValueError("residual_sd must be positive.")
    if not 0 < alpha < 1:
        raise ValueError("alpha must lie strictly between 0 and 1.")

    noncentrality = (
        abs(interaction_effect_sd)
        * np.sqrt(sample_size * signature_prevalence * (1 - signature_prevalence))
        / residual_sd
    )
    critical = norm.ppf(1 - alpha / 2)
    return float(norm.sf(critical - noncentrality) + norm.cdf(-critical - noncentrality))


def randomized_preconditioning_interaction_power(
    sample_size: int,
    genotype_by_preconditioning_effect_sd: float,
    within_person_difference_sd: float = 1.0,
    alpha: float = 0.05,
) -> float:
    """Power for a two-genotype difference-in-differences challenge design.

    Each participant is measured under both preconditioning states, and the
    within-person response difference is compared between equally sized genotype
    strata. The target effect and difference standard deviation must use the same
    endpoint scale.
    """
    if sample_size < 4:
        raise ValueError("sample_size must be at least 4.")
    if within_person_difference_sd <= 0:
        raise ValueError("within_person_difference_sd must be positive.")
    if not 0 < alpha < 1:
        raise ValueError("alpha must lie strictly between 0 and 1.")
    group_zero = sample_size // 2
    group_one = sample_size - group_zero
    standard_error = within_person_difference_sd * np.sqrt(1 / group_zero + 1 / group_one)
    noncentrality = abs(genotype_by_preconditioning_effect_sd) / standard_error
    critical = norm.ppf(1 - alpha / 2)
    return float(norm.sf(critical - noncentrality) + norm.cdf(-critical - noncentrality))


def required_randomized_preconditioning_sample_size(
    genotype_by_preconditioning_effect_sd: float,
    within_person_difference_sd: float = 1.0,
    alpha: float = 0.05,
    target_power: float = 0.8,
    attrition_fraction: float = 0.15,
) -> int:
    """Return an even recruitment target for the difference-in-differences test."""
    if genotype_by_preconditioning_effect_sd == 0:
        raise ValueError("The target interaction effect must be non-zero.")
    if within_person_difference_sd <= 0:
        raise ValueError("within_person_difference_sd must be positive.")
    if not 0 < alpha < 1 or not 0 < target_power < 1:
        raise ValueError("alpha and target_power must lie strictly between 0 and 1.")
    if not 0 <= attrition_fraction < 1:
        raise ValueError("attrition_fraction must lie in [0, 1).")
    critical = norm.ppf(1 - alpha / 2)
    target = norm.ppf(target_power)
    analyzable = (
        4
        * within_person_difference_sd**2
        * (critical + target) ** 2
        / genotype_by_preconditioning_effect_sd**2
    )
    recruitment = int(np.ceil(analyzable / (1 - attrition_fraction)))
    return recruitment + recruitment % 2


def monte_carlo_power_normal(
    noncentrality: float, alpha: float, replicates: int, rng: np.random.Generator
) -> tuple[float, float]:
    """Simulate normal test statistics and return rejection-rate estimate and SE."""
    if replicates <= 0:
        raise ValueError("replicates must be positive.")
    critical = norm.ppf(1 - alpha / 2)
    test_statistics = rng.normal(loc=noncentrality, scale=1.0, size=replicates)
    rejected = np.abs(test_statistics) > critical
    estimate = float(rejected.mean())
    standard_error = float(np.sqrt(estimate * (1 - estimate) / replicates))
    return estimate, standard_error


def make_power_grid(
    signature_scenarios: Sequence[dict[str, object]],
    cohort_sizes: Iterable[int],
    interaction_effects_sd: Iterable[float],
    residual_sd: float,
    discovery_alpha: float,
    replication_alpha: float,
    monte_carlo_replicates: int,
    discovery_fraction: float,
    minimum_carriers_discovery: int,
    minimum_carriers_replication: int,
    random_seed: int,
) -> pd.DataFrame:
    """Build the registered feasibility grid from configuration values."""
    if not 0 < discovery_fraction <= 1:
        raise ValueError("discovery_fraction must lie in (0, 1].")

    rng = np.random.default_rng(random_seed)
    records: list[dict[str, object]] = []
    for scenario in signature_scenarios:
        scenario_name = str(scenario["name"])
        component_frequencies = [float(value) for value in scenario["component_frequencies"]]
        prevalence = joint_signature_prevalence(component_frequencies)
        for total_size in cohort_sizes:
            discovery_size = int(np.floor(int(total_size) * discovery_fraction))
            replication_size = int(total_size) - discovery_size
            expected_discovery_carriers = discovery_size * prevalence
            expected_replication_carriers = replication_size * prevalence
            for effect in interaction_effects_sd:
                discovery_power = interaction_power_normal(
                    sample_size=discovery_size,
                    signature_prevalence=prevalence,
                    interaction_effect_sd=float(effect),
                    residual_sd=residual_sd,
                    alpha=discovery_alpha,
                )
                replication_power = interaction_power_normal(
                    sample_size=replication_size,
                    signature_prevalence=prevalence,
                    interaction_effect_sd=float(effect),
                    residual_sd=residual_sd,
                    alpha=replication_alpha,
                )
                discovery_noncentrality = (
                    abs(float(effect))
                    * np.sqrt(discovery_size * prevalence * (1 - prevalence))
                    / residual_sd
                )
                replication_noncentrality = (
                    abs(float(effect))
                    * np.sqrt(replication_size * prevalence * (1 - prevalence))
                    / residual_sd
                )
                discovery_mc_power, discovery_mc_se = monte_carlo_power_normal(
                    discovery_noncentrality,
                    discovery_alpha,
                    monte_carlo_replicates,
                    rng,
                )
                replication_mc_power, replication_mc_se = monte_carlo_power_normal(
                    replication_noncentrality,
                    replication_alpha,
                    monte_carlo_replicates,
                    rng,
                )
                passes_carrier_gate = (
                    expected_discovery_carriers >= minimum_carriers_discovery
                    and expected_replication_carriers >= minimum_carriers_replication
                )
                records.append(
                    {
                        "signature_scenario": scenario_name,
                        "component_frequencies": ";".join(map(str, component_frequencies)),
                        "joint_prevalence": prevalence,
                        "total_cohort_size": int(total_size),
                        "discovery_size": discovery_size,
                        "replication_size": replication_size,
                        "expected_discovery_carriers": expected_discovery_carriers,
                        "expected_replication_carriers": expected_replication_carriers,
                        "interaction_effect_sd": float(effect),
                        "discovery_alpha": discovery_alpha,
                        "replication_alpha": replication_alpha,
                        "discovery_power": discovery_power,
                        "replication_power": replication_power,
                        "discovery_mc_power": discovery_mc_power,
                        "discovery_mc_se": discovery_mc_se,
                        "replication_mc_power": replication_mc_power,
                        "replication_mc_se": replication_mc_se,
                        "passes_carrier_gate": passes_carrier_gate,
                        "passes_power_80": discovery_power >= 0.80 and replication_power >= 0.80,
                    }
                )
    return pd.DataFrame.from_records(records)
