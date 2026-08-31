from __future__ import annotations

from genotype_gated_metabolism.pipelines.pathway_score_mapping_supplement import (
    feature_exclusion_reason,
)


def test_feature_exclusion_reason_is_mutually_exclusive() -> None:
    assert (
        feature_exclusion_reason(
            complete=True, accepted=True, recognized_family=True
        )
        == "accepted"
    )
    assert (
        feature_exclusion_reason(
            complete=False, accepted=False, recognized_family=True
        )
        == "incomplete_across_eligible_samples"
    )
    assert (
        feature_exclusion_reason(
            complete=True, accepted=False, recognized_family=False
        )
        == "unrecognized_or_undersized_lipid_family"
    )
