"""Reusable leakage-control, validation and falsification primitives."""

from .validation import (
    GroupedValidationSpec,
    assert_group_isolation,
    regression_metrics,
    repeated_balanced_group_splits,
)

__all__ = [
    "GroupedValidationSpec",
    "assert_group_isolation",
    "regression_metrics",
    "repeated_balanced_group_splits",
]
