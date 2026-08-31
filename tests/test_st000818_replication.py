from __future__ import annotations

import pandas as pd

from genotype_gated_metabolism.pipelines.pathway_score_st000818_replication import (
    adjudicate_external_replication,
)


def test_external_replication_requires_all_nulls_and_eligibility() -> None:
    results = pd.DataFrame(
        {
            "row_degrees_preserved": [True, True],
            "column_degrees_preserved": [True, True],
            "rmse_ci_lower": [0.05, 0.02],
            "precision_ci_lower": [0.03, 0.01],
        }
    )
    decision, gates = adjudicate_external_replication(
        results,
        eligible_features=120,
        validation_groups=15,
        minimum_features=100,
        minimum_groups=10,
        threshold=0.0,
    )
    assert decision == "EXTERNAL_STRUCTURAL_REPLICATION"
    assert all(gates.values())

    results.loc[1, "precision_ci_lower"] = -0.01
    decision, _ = adjudicate_external_replication(
        results,
        eligible_features=120,
        validation_groups=15,
        minimum_features=100,
        minimum_groups=10,
        threshold=0.0,
    )
    assert decision == "EXTERNAL_STRUCTURAL_RECONSTRUCTION_ONLY"
