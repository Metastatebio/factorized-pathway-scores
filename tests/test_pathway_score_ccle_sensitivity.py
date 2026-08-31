from __future__ import annotations

import pandas as pd

from genotype_gated_metabolism.pipelines.pathway_score_ccle_sensitivity import (
    adjudicate_ccle_sensitivity,
)


def test_ccle_sensitivity_requires_positive_random_control_results() -> None:
    summary = pd.DataFrame(
        {
            "rmse_improvement_vs_random_sd": [0.1, 0.2],
            "random_ci_lower": [0.02, 0.04],
            "median_factorized_target_r2": [0.2, 0.1],
            "rmse_improvement_vs_additive_sd": [-0.01, 0.01],
        }
    )
    decision, result = adjudicate_ccle_sensitivity(
        summary, point_threshold=0.0, minimum_ci_pass_rate=0.8
    )
    assert decision == "CCLE_REACTION_TOPOLOGY_SENSITIVITY_ROBUST"
    assert result["interaction_beats_additive_rate"] == 0.5
    summary.loc[0, "rmse_improvement_vs_random_sd"] = -0.01
    decision, _ = adjudicate_ccle_sensitivity(
        summary, point_threshold=0.0, minimum_ci_pass_rate=0.8
    )
    assert decision == "SENSITIVITY_WARNING"
