from __future__ import annotations

from genotype_gated_metabolism.pipelines.pathway_score_ccle_reaction_cluster import (
    adjudicate_cluster_robustness,
)


def test_cluster_adjudication_requires_interval_and_leave_one_out() -> None:
    decision, audit = adjudicate_cluster_robustness(
        clusters=12,
        ci_lower=0.04,
        minimum_leave_one_cluster_effect=0.03,
        minimum_clusters=10,
        threshold=0.0,
    )
    assert decision == "CCLE_REACTION_CLUSTER_ROBUST"
    assert audit["minimum_clusters_pass"]

    decision, _ = adjudicate_cluster_robustness(
        clusters=12,
        ci_lower=-0.01,
        minimum_leave_one_cluster_effect=0.03,
        minimum_clusters=10,
        threshold=0.0,
    )
    assert decision == "REACTION_CLUSTER_WARNING"
