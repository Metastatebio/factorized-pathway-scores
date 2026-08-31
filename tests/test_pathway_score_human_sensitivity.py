from __future__ import annotations

import pandas as pd

from genotype_gated_metabolism.pipelines.pathway_score_human_sensitivity import (
    adjudicate_human_sensitivity,
)


def test_human_sensitivity_requires_direction_and_ci_pass_rate() -> None:
    results = pd.DataFrame(
        {
            "dataset": ["a", "a", "b", "b"],
            "row_degrees_preserved": [True] * 4,
            "column_degrees_preserved": [True] * 4,
            "rmse_improvement_sd": [0.1] * 4,
            "precision_improvement": [0.1] * 4,
            "rmse_ci_lower": [0.01] * 4,
            "precision_ci_lower": [0.01] * 4,
        }
    )
    decision, _ = adjudicate_human_sensitivity(
        results, point_threshold=0.0, minimum_ci_pass_rate=0.8
    )
    assert decision == "HUMAN_STRUCTURAL_SENSITIVITY_ROBUST"
    results.loc[0, "precision_improvement"] = -0.01
    decision, _ = adjudicate_human_sensitivity(
        results, point_threshold=0.0, minimum_ci_pass_rate=0.8
    )
    assert decision == "SENSITIVITY_WARNING"
