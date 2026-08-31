from __future__ import annotations

import pandas as pd

from genotype_gated_metabolism.pipelines.pathway_score_st002081_structural_null_sensitivity import (
    adjudicate_null_ensemble,
)


def test_null_ensemble_requires_every_graph_to_pass_both_endpoints() -> None:
    passing = pd.DataFrame(
        {
            "row_degrees_preserved": [True, True],
            "column_degrees_preserved": [True, True],
            "rmse_ci_lower": [0.08, 0.04],
            "precision_ci_lower": [0.1, 0.02],
        }
    )
    decision, gates = adjudicate_null_ensemble(passing, threshold=0.0)
    assert decision == "ADAPTIVE_STRUCTURE_SIGNAL_ROBUST_ACROSS_GRAPH_NULLS"
    assert all(gates.values())

    failing = passing.copy()
    failing.loc[1, "precision_ci_lower"] = -0.01
    decision, gates = adjudicate_null_ensemble(failing, threshold=0.0)
    assert decision == "ADAPTIVE_STRUCTURE_GRAPH_NULL_SENSITIVITY_FAILED"
    assert not gates["precision_ci_passes_all_nulls"]
