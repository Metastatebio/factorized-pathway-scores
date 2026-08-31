"""Model-agnostic grouped validation contracts for cross-omic prediction."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class GroupedValidationSpec:
    """Serializable contract for repeated, group-isolated validation."""

    outer_folds: int = 5
    repeats: int = 1
    seed: int = 20260830

    def __post_init__(self) -> None:
        if self.outer_folds < 2:
            raise ValueError("Grouped validation requires at least two folds.")
        if self.repeats < 1:
            raise ValueError("Grouped validation requires at least one repeat.")


def _balanced_group_assignment(
    groups: np.ndarray, *, folds: int, seed: int
) -> dict[object, int]:
    unique, counts = np.unique(groups, return_counts=True)
    if len(unique) < folds:
        raise ValueError("Fewer unique groups than requested folds.")
    rng = np.random.default_rng(seed)
    tie_break = rng.random(len(unique))
    order = np.lexsort((tie_break, -counts))
    fold_sizes = np.zeros(folds, dtype=int)
    assignment: dict[object, int] = {}
    for index in order:
        smallest = np.flatnonzero(fold_sizes == fold_sizes.min())
        fold = int(rng.choice(smallest))
        assignment[unique[index]] = fold
        fold_sizes[fold] += int(counts[index])
    return assignment


def repeated_balanced_group_splits(
    groups: pd.Series | np.ndarray | list[object],
    spec: GroupedValidationSpec,
) -> list[tuple[np.ndarray, np.ndarray, int, int]]:
    """Return reproducible splits that keep every group wholly inside one outer fold."""
    values = np.asarray(groups, dtype=object)
    if values.ndim != 1 or len(values) == 0:
        raise ValueError("groups must be a non-empty one-dimensional vector.")
    if pd.isna(values).any():
        raise ValueError("groups cannot contain missing identifiers.")
    splits: list[tuple[np.ndarray, np.ndarray, int, int]] = []
    for repeat in range(spec.repeats):
        assignment = _balanced_group_assignment(
            values, folds=spec.outer_folds, seed=spec.seed + repeat
        )
        fold_ids = np.asarray([assignment[value] for value in values], dtype=int)
        for fold in range(spec.outer_folds):
            test = np.flatnonzero(fold_ids == fold)
            train = np.flatnonzero(fold_ids != fold)
            splits.append((train, test, repeat, fold))
    assert_group_isolation(values, splits)
    return splits


def assert_group_isolation(
    groups: pd.Series | np.ndarray | list[object],
    splits: list[tuple[np.ndarray, np.ndarray, int, int]],
) -> None:
    """Raise if any split leaks a group or fails to cover every row exactly once per repeat."""
    values = np.asarray(groups, dtype=object)
    repeats: dict[int, np.ndarray] = {}
    for train, test, repeat, _ in splits:
        if set(values[train]) & set(values[test]):
            raise ValueError("Group leakage detected between training and test rows.")
        repeats.setdefault(repeat, np.zeros(len(values), dtype=int))[test] += 1
    if not repeats or any(not np.all(coverage == 1) for coverage in repeats.values()):
        raise ValueError("Every row must appear in exactly one test fold per repeat.")


def regression_metrics(
    observed: np.ndarray,
    predicted: np.ndarray,
    *,
    groups: pd.Series | np.ndarray | list[object],
) -> dict[str, float | int]:
    """Report row-weighted and equal-group-weighted regression estimands separately."""
    truth = np.asarray(observed, dtype=float)
    estimate = np.asarray(predicted, dtype=float)
    group_values = np.asarray(groups, dtype=object)
    if truth.shape != estimate.shape or truth.shape[0] != len(group_values):
        raise ValueError("Observed, predicted and group rows must align.")
    if truth.ndim == 1:
        truth = truth[:, None]
        estimate = estimate[:, None]
    finite = np.isfinite(truth) & np.isfinite(estimate)
    if not finite.any() or not np.all(finite.any(axis=1)):
        raise ValueError("Every row must contain at least one finite observed-predicted pair.")
    row_mse = np.divide(
        np.where(finite, (truth - estimate) ** 2, 0.0).sum(axis=1),
        finite.sum(axis=1),
    )
    group_mse = pd.Series(row_mse).groupby(group_values).mean()
    centered = truth - np.nanmean(truth, axis=0)
    sse = float(np.where(finite, (truth - estimate) ** 2, 0.0).sum())
    sst = float(np.where(finite, centered**2, 0.0).sum())
    return {
        "rows": len(truth),
        "groups": len(group_mse),
        "row_weighted_rmse": float(np.sqrt(row_mse.mean())),
        "equal_group_weighted_rmse": float(np.sqrt(group_mse.mean())),
        "row_weighted_r2": float(1.0 - sse / sst) if sst > 0 else np.nan,
    }
