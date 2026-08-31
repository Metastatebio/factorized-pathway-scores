from __future__ import annotations

import pandas as pd

from genotype_gated_metabolism.pipelines.pathway_score_ccle_null_ensemble import (
    adjudicate_null_ensemble,
)


def test_null_ensemble_requires_seed_and_expected_null_robustness() -> None:
    seeds = pd.DataFrame(
        {
            "improvement_sd": [0.10, 0.12, 0.08],
            "ci_lower": [0.02, 0.03, 0.01],
            "factorized_rmse_sd": [0.8, 0.8, 0.8],
            "targets": [60, 60, 60],
            "target_digest": ["same", "same", "same"],
        }
    )
    decision, audit = adjudicate_null_ensemble(
        seeds,
        ensemble_ci_lower=0.02,
        minimum_seed_ci_pass_rate=0.9,
        improvement_threshold=0.0,
        invariant_tolerance=1e-12,
    )
    assert decision == "CCLE_RANDOM_FEATURE_NULL_ENSEMBLE_ROBUST"
    assert audit["seed_ci_pass_rate"] == 1.0

    seeds.loc[0, "improvement_sd"] = -0.01
    decision, _ = adjudicate_null_ensemble(
        seeds,
        ensemble_ci_lower=0.02,
        minimum_seed_ci_pass_rate=0.9,
        improvement_threshold=0.0,
        invariant_tolerance=1e-12,
    )
    assert decision == "NULL_ENSEMBLE_WARNING"
