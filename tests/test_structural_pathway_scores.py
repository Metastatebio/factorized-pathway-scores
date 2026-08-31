from __future__ import annotations

import numpy as np
import pandas as pd

from genotype_gated_metabolism.analysis.pathway_scores import (
    degree_preserving_descriptor_null,
    disjoint_family_masks,
    eligible_lipid_families,
    lipid_structure_descriptors,
    masked_structural_descriptor_benchmark,
    structural_descriptor_incidence,
)


def test_lipid_descriptors_parse_chain_and_total_composition() -> None:
    pc = set(lipid_structure_descriptors("PC(16:0_18:2)"))
    tag = set(lipid_structure_descriptors("TAG52:4-FA18:1"))
    workbench_dg = set(lipid_structure_descriptors("DG 34:2; [M+NH4]+@7.52"))
    workbench_lyso = set(lipid_structure_descriptors("lysoPC 18:1; [M+H]+@1.04"))
    assert {"family:PC", "total_carbon:34", "total_unsaturation:2", "acyl:18:2"} <= pc
    assert {"family:TAG", "total_carbon:52", "total_unsaturation:4", "acyl:18:1"} <= tag
    assert {"family:DAG", "total_carbon:34", "total_unsaturation:2"} <= workbench_dg
    assert {"family:LPC", "total_carbon:18", "total_unsaturation:1"} <= workbench_lyso


def test_degree_preserving_null_keeps_both_degree_sequences() -> None:
    incidence = pd.DataFrame(
        [
            [1, 1, 0, 0],
            [1, 0, 1, 0],
            [0, 1, 0, 1],
            [0, 0, 1, 1],
        ],
        index=list("abcd"),
        columns=list("wxyz"),
        dtype=bool,
    )
    randomized = degree_preserving_descriptor_null(incidence, swaps_per_edge=50, seed=7)
    assert np.array_equal(incidence.sum(axis=0), randomized.sum(axis=0))
    assert np.array_equal(incidence.sum(axis=1), randomized.sum(axis=1))


def test_structural_benchmark_emits_degree_audit() -> None:
    rng = np.random.default_rng(9)
    features = pd.Index(
        [
            *[f"PC(16:0_18:{index % 3})" for index in range(10)],
            *[f"TAG{48 + index % 4}:{index % 3}-FA18:{index % 2}" for index in range(10)],
        ]
    )
    # Make names unique while retaining parseable structure.
    features = pd.Index([f"{value}_{index}" for index, value in enumerate(features)])
    values = pd.DataFrame(rng.normal(size=(30, 20)), columns=features)
    subjects = pd.Series(np.repeat([f"p{i}" for i in range(10)], 3), index=values.index)
    families = eligible_lipid_families(features, minimum_family_features=5)
    masks = disjoint_family_masks(families, masks=5, seed=11)
    incidence = structural_descriptor_incidence(
        features, minimum_features=2, maximum_feature_fraction=0.95
    )
    metrics, samples, audit = masked_structural_descriptor_benchmark(
        values,
        subjects,
        masks,
        incidence,
        folds=5,
        ridge_alpha=1.0,
        split_seed=13,
        null_seed=17,
        swaps_per_edge=5,
        priority_count=2,
        true_priority_fraction=0.2,
    )
    assert set(metrics["model"]) == {
        "population_mean",
        "degree_preserving_random_structural_ridge",
        "structural_descriptor_ridge",
    }
    assert audit["row_degrees_preserved"].all()
    assert audit["column_degrees_preserved"].all()
    assert samples["subject_id"].nunique() == 10
