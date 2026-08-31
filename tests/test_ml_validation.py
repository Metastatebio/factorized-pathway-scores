import numpy as np

from genotype_gated_metabolism.ml.validation import (
    GroupedValidationSpec,
    regression_metrics,
    repeated_balanced_group_splits,
)


def test_repeated_group_splits_are_isolated_and_deterministic() -> None:
    groups = np.repeat(["a", "b", "c", "d", "e", "f"], [8, 2, 7, 3, 6, 4])
    spec = GroupedValidationSpec(outer_folds=3, repeats=2, seed=19)
    first = repeated_balanced_group_splits(groups, spec)
    second = repeated_balanced_group_splits(groups, spec)
    assert len(first) == 6
    for left, right in zip(first, second, strict=True):
        assert np.array_equal(left[0], right[0])
        assert np.array_equal(left[1], right[1])
        assert left[2:] == right[2:]
        assert not (set(groups[left[0]]) & set(groups[left[1]]))


def test_regression_metrics_separate_row_and_group_estimands() -> None:
    observed = np.array([0.0, 0.0, 0.0, 10.0])
    predicted = np.array([0.0, 0.0, 0.0, 0.0])
    metrics = regression_metrics(observed, predicted, groups=["many", "many", "many", "one"])
    assert metrics["rows"] == 4
    assert metrics["groups"] == 2
    assert metrics["equal_group_weighted_rmse"] > metrics["row_weighted_rmse"]
