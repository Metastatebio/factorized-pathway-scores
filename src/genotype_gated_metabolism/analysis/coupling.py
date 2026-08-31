"""Robust inference for continuous cross-omic coupling modifiers."""

from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np
import pandas as pd
import statsmodels.api as sm


@dataclass(frozen=True)
class CouplingResult:
    n: int
    groups: int
    estimate: float
    standard_error: float
    p_value: float
    confidence_interval_low: float
    confidence_interval_high: float
    slope_at_low_modifier: float
    slope_at_high_modifier: float
    covariance_type: str

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class InfluenceRefitResult:
    """Interaction refit after removing the largest Cook's-distance values."""

    removed_n: int
    maximum_cooks_distance: float
    estimate: float
    p_value: float

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


def _design_matrix(
    frame: pd.DataFrame,
    exposure: str,
    modifier: str,
    group: str | None,
    control_group_slopes: bool,
) -> tuple[pd.DataFrame, str]:
    interaction = f"{exposure}:x:{modifier}"
    design = pd.DataFrame(index=frame.index)
    design[exposure] = frame[exposure]
    design[modifier] = frame[modifier]
    design[interaction] = frame[exposure] * frame[modifier]

    if group is not None:
        dummies = pd.get_dummies(frame[group].astype(str), prefix=group, drop_first=True)
        dummies = dummies.astype(float)
        design = pd.concat([design, dummies], axis=1)
        if control_group_slopes:
            for column in dummies:
                design[f"{exposure}:x:{column}"] = frame[exposure] * dummies[column]
                design[f"{modifier}:x:{column}"] = frame[modifier] * dummies[column]
    return sm.add_constant(design, has_constant="add"), interaction


def fit_coupling(
    data: pd.DataFrame,
    outcome: str,
    exposure: str,
    modifier: str,
    group: str | None = None,
    covariance_type: str = "HC3",
    control_group_slopes: bool = True,
) -> CouplingResult:
    """Fit B ~ A + G + A:G with optional group effects and lower-order slopes."""
    required = [outcome, exposure, modifier, *([group] if group else [])]
    missing = sorted(set(required) - set(data.columns))
    if missing:
        raise ValueError(f"Missing columns: {', '.join(missing)}")
    frame = data.loc[:, required].replace([np.inf, -np.inf], np.nan).dropna().copy()
    if len(frame) < 20:
        raise ValueError("At least 20 complete samples are required.")
    for column in (outcome, exposure, modifier):
        if frame[column].nunique() < 2:
            raise ValueError(f"Column {column} has no usable variation.")

    design, interaction = _design_matrix(frame, exposure, modifier, group, control_group_slopes)
    if design.shape[0] <= design.shape[1] + 5:
        raise ValueError("The interaction design has insufficient residual degrees of freedom.")
    fitted = sm.OLS(frame[outcome], design).fit(cov_type=covariance_type)
    interval = fitted.conf_int().loc[interaction]
    main_slope = float(fitted.params[exposure])
    interaction_estimate = float(fitted.params[interaction])
    return CouplingResult(
        n=len(frame),
        groups=int(frame[group].nunique()) if group else 1,
        estimate=interaction_estimate,
        standard_error=float(fitted.bse[interaction]),
        p_value=float(fitted.pvalues[interaction]),
        confidence_interval_low=float(interval.iloc[0]),
        confidence_interval_high=float(interval.iloc[1]),
        slope_at_low_modifier=main_slope - interaction_estimate,
        slope_at_high_modifier=main_slope + interaction_estimate,
        covariance_type=covariance_type,
    )


def refit_without_most_influential(
    data: pd.DataFrame,
    outcome: str,
    exposure: str,
    modifier: str,
    group: str | None = None,
    exclusion_fraction: float = 0.01,
    covariance_type: str = "HC3",
    control_group_slopes: bool = True,
) -> InfluenceRefitResult:
    """Remove the largest Cook's distances and refit the interaction model."""
    if not 0 < exclusion_fraction < 0.5:
        raise ValueError("exclusion_fraction must be between 0 and 0.5.")
    required = [outcome, exposure, modifier, *([group] if group else [])]
    frame = data.loc[:, required].replace([np.inf, -np.inf], np.nan).dropna().copy()
    design, _ = _design_matrix(frame, exposure, modifier, group, control_group_slopes)
    fitted = sm.OLS(frame[outcome], design).fit()
    cooks_distance = np.asarray(fitted.get_influence().cooks_distance[0])
    removed_n = max(1, int(np.ceil(len(frame) * exclusion_fraction)))
    retained = np.ones(len(frame), dtype=bool)
    retained[np.argsort(cooks_distance)[-removed_n:]] = False
    refit = fit_coupling(
        frame.iloc[np.flatnonzero(retained)],
        outcome,
        exposure,
        modifier,
        group=group,
        covariance_type=covariance_type,
        control_group_slopes=control_group_slopes,
    )
    return InfluenceRefitResult(
        removed_n=removed_n,
        maximum_cooks_distance=float(np.nanmax(cooks_distance)),
        estimate=refit.estimate,
        p_value=refit.p_value,
    )


def permutation_p_value(
    data: pd.DataFrame,
    outcome: str,
    exposure: str,
    modifier: str,
    observed_estimate: float,
    group: str | None,
    permutations: int,
    random_seed: int,
    control_group_slopes: bool = True,
) -> float:
    """Compute a two-sided permutation P value, shuffling within groups."""
    if permutations <= 0:
        raise ValueError("permutations must be positive.")
    rng = np.random.default_rng(random_seed)
    required = [outcome, exposure, modifier, *([group] if group else [])]
    frame = data.loc[:, required].replace([np.inf, -np.inf], np.nan).dropna().copy()
    exposure_values = frame[exposure].to_numpy(dtype=float)
    modifier_values = frame[modifier].to_numpy(dtype=float)
    outcome_values = frame[outcome].to_numpy(dtype=float)

    if group is None:
        dummies = np.empty((len(frame), 0), dtype=float)
        group_indices = [np.arange(len(frame))]
    else:
        encoded = pd.get_dummies(frame[group].astype(str), drop_first=True).astype(float)
        dummies = encoded.to_numpy()
        group_indices = [
            np.flatnonzero(frame[group].to_numpy() == value) for value in frame[group].unique()
        ]

    fixed_columns = [np.ones(len(frame)), exposure_values]
    if dummies.shape[1]:
        fixed_columns.extend(dummies.T)
        if control_group_slopes:
            fixed_columns.extend((exposure_values[:, None] * dummies).T)
    fixed = np.column_stack(fixed_columns)

    def interaction_estimate(shuffled_modifier: np.ndarray) -> float:
        variable_columns = [shuffled_modifier, exposure_values * shuffled_modifier]
        interaction_index = fixed.shape[1] + 1
        if dummies.shape[1] and control_group_slopes:
            variable_columns.extend((shuffled_modifier[:, None] * dummies).T)
        design = np.column_stack([fixed, *variable_columns])
        coefficients = np.linalg.lstsq(design, outcome_values, rcond=None)[0]
        return float(coefficients[interaction_index])

    extreme = 0
    for _ in range(permutations):
        shuffled_modifier = modifier_values.copy()
        for indices in group_indices:
            shuffled_modifier[indices] = rng.permutation(modifier_values[indices])
        estimate = interaction_estimate(shuffled_modifier)
        extreme += abs(estimate) >= abs(observed_estimate)
    return (extreme + 1) / (permutations + 1)
