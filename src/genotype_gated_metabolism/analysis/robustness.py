"""Reusable sensitivity analyses for cross-omic interaction models."""

from __future__ import annotations

import numpy as np
import pandas as pd

from .coupling import fit_coupling


def winsorize_within_groups(
    data: pd.DataFrame,
    columns: list[str],
    group: str,
    tail_fraction: float = 0.01,
) -> pd.DataFrame:
    """Clip numeric columns to group-specific empirical tail quantiles."""
    if not 0 <= tail_fraction < 0.5:
        raise ValueError("tail_fraction must be in [0, 0.5).")
    frame = data.copy()
    for column in columns:
        frame[column] = frame.groupby(group, observed=True)[column].transform(
            lambda values: values.clip(
                lower=values.quantile(tail_fraction),
                upper=values.quantile(1 - tail_fraction),
            )
        )
    return frame


def leave_one_group_out_summary(
    data: pd.DataFrame,
    outcome: str,
    exposure: str,
    modifier: str,
    group: str,
    reference_estimate: float,
    covariance_type: str = "HC3",
    control_group_slopes: bool = True,
) -> dict[str, float | int]:
    """Summarize interaction stability after excluding each group in turn."""
    estimates: list[float] = []
    for value in sorted(data[group].dropna().astype(str).unique()):
        subset = data.loc[data[group].astype(str).ne(value)]
        result = fit_coupling(
            subset,
            outcome,
            exposure,
            modifier,
            group=group,
            covariance_type=covariance_type,
            control_group_slopes=control_group_slopes,
        )
        estimates.append(result.estimate)
    estimate_array = np.asarray(estimates)
    return {
        "leave_one_group_out_runs": len(estimates),
        "leave_one_group_out_sign_agreement": float(
            np.mean(estimate_array * reference_estimate > 0)
        ),
        "leave_one_group_out_min_estimate": float(estimate_array.min()),
        "leave_one_group_out_max_estimate": float(estimate_array.max()),
    }
