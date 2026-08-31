from __future__ import annotations

import pandas as pd

from genotype_gated_metabolism.pipelines.pathway_score_human_graph_mixing import (
    adjudicate_graph_mixing,
)


def test_graph_mixing_requires_degrees_distance_and_uniqueness() -> None:
    nulls = pd.DataFrame(
        {
            "row_degrees_preserved": [True, True],
            "column_degrees_preserved": [True, True],
            "edge_replacement_fraction": [0.7, 0.8],
        }
    )
    masks = pd.DataFrame(
        {"unique_nulls": [20, 20], "expected_nulls": [20, 20], "maximum_pairwise_jaccard": [0.3, 0.4]}
    )
    decision, audit = adjudicate_graph_mixing(
        nulls,
        masks,
        minimum_edge_replacement_fraction=0.5,
        maximum_pairwise_null_jaccard=0.8,
    )
    assert decision == "HUMAN_GRAPH_NULL_MIXING_ADEQUATE"
    assert audit["all_nulls_unique"]

    nulls.loc[0, "edge_replacement_fraction"] = 0.2
    decision, _ = adjudicate_graph_mixing(
        nulls,
        masks,
        minimum_edge_replacement_fraction=0.5,
        maximum_pairwise_null_jaccard=0.8,
    )
    assert decision == "GRAPH_NULL_MIXING_WARNING"
