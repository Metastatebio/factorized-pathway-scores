from __future__ import annotations

import numpy as np
import pandas as pd

from genotype_gated_metabolism.analysis.pathway_scores import (
    bootstrap_subject_model_difference,
    disjoint_family_masks,
    eligible_lipid_families,
    masked_panel_benchmark,
)


def test_disjoint_masks_cover_every_eligible_feature_once() -> None:
    features = pd.Index(
        [*[f"PC({index})" for index in range(10)], *[f"TAG{index}" for index in range(10)]]
    )
    families = eligible_lipid_families(features, minimum_family_features=5)
    masks = disjoint_family_masks(families, masks=5, seed=7)
    assert set(masks["feature"]) == set(features)
    assert masks["feature"].value_counts().eq(1).all()
    assert masks.groupby(["family", "mask"]).size().eq(2).all()


def test_masked_panel_benchmark_is_subject_isolated_and_finite() -> None:
    rng = np.random.default_rng(11)
    subjects = np.repeat([f"p{index}" for index in range(10)], 2)
    latent = np.repeat(rng.normal(size=(10, 2)), 2, axis=0)
    values = pd.DataFrame(
        {
            **{f"PC({index})": latent[:, 0] + rng.normal(scale=0.1, size=20) for index in range(5)},
            **{f"TAG{index}": latent[:, 1] + rng.normal(scale=0.1, size=20) for index in range(5)},
        },
        index=[f"s{index}" for index in range(20)],
    )
    families = eligible_lipid_families(values.columns, minimum_family_features=5)
    masks = disjoint_family_masks(families, masks=5, seed=13)
    metrics, samples, targets = masked_panel_benchmark(
        values,
        pd.Series(subjects, index=values.index),
        families,
        masks,
        folds=5,
        ridge_alpha=1.0,
        random_group_seed=17,
        priority_count=1,
        true_priority_fraction=0.2,
    )
    assert set(metrics["model"]) == {
        "population_mean",
        "random_group_score_ridge",
        "family_score_ridge",
        "all_visible_ridge",
    }
    assert samples.groupby(["subject_id", "mask"])["fold"].nunique().eq(1).all()
    assert targets.groupby(["feature", "model"]).size().eq(1).all()
    assert np.isfinite(metrics["row_weighted_rmse_sd"]).all()


def test_subject_bootstrap_returns_paired_intervals() -> None:
    rows = []
    for subject in range(12):
        rows.extend(
            [
                {
                    "subject_id": f"p{subject}",
                    "model": "reference",
                    "mean_squared_error": 1.0,
                    "precision_at_k": 0.2,
                },
                {
                    "subject_id": f"p{subject}",
                    "model": "challenger",
                    "mean_squared_error": 0.25,
                    "precision_at_k": 0.6,
                },
            ]
        )
    result = bootstrap_subject_model_difference(
        pd.DataFrame(rows),
        reference="reference",
        challenger="challenger",
        draws=200,
        seed=19,
    )
    assert result["rmse_ci_lower"] > 0
    assert result["precision_ci_lower"] > 0
