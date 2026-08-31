"""Mechanism-aware feature builders for one- and multi-gene signatures."""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
import pandas as pd


def parse_gene_signature(value: str) -> tuple[str, ...]:
    """Parse the semicolon-delimited canonical signature representation."""
    genes = tuple(sorted({part.strip() for part in value.split(";") if part.strip()}))
    if not genes:
        raise ValueError("A gene signature must contain at least one gene.")
    return genes


def build_expression_signature(
    expression: pd.DataFrame,
    genes: Sequence[str],
    aggregation: str = "limiting_subunit",
) -> pd.Series:
    """Aggregate log-expression into a reaction-state proxy.

    Candidate generation emits one row for each GPR disjunct. Multiple genes in
    a row therefore represent an enzyme complex (AND logic); the default proxy
    uses the least-expressed required subunit. Single genes pass through.
    """
    canonical = tuple(sorted(set(genes)))
    if not canonical:
        raise ValueError("At least one gene is required.")
    missing = sorted(set(canonical) - set(expression.columns))
    if missing:
        raise KeyError(f"Missing expression genes: {', '.join(missing)}")

    values = expression.loc[:, canonical]
    if len(canonical) == 1:
        signature = values.iloc[:, 0]
    elif aggregation == "limiting_subunit":
        signature = values.min(axis=1)
    elif aggregation == "mean":
        signature = values.mean(axis=1)
    elif aggregation == "sum":
        signature = values.sum(axis=1)
    else:
        raise ValueError(f"Unsupported expression-signature aggregation: {aggregation}")
    signature = signature.astype(float).copy()
    signature.name = ";".join(canonical)
    return signature


def zscore_within_groups(values: pd.Series, groups: pd.Series) -> pd.Series:
    """Standardize a feature within context groups, retaining missing values."""
    if not values.index.equals(groups.index):
        groups = groups.reindex(values.index)

    def standardize(group: pd.Series) -> pd.Series:
        scale = group.std(ddof=0)
        if not np.isfinite(scale) or scale == 0:
            return pd.Series(np.nan, index=group.index, dtype=float)
        return (group - group.mean()) / scale

    result = values.groupby(groups, observed=True, group_keys=False).apply(standardize)
    result.name = values.name
    return result
