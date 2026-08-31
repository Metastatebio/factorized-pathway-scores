from __future__ import annotations

import numpy as np
import pandas as pd

from genotype_gated_metabolism.analysis.ccle_pathway_prediction import (
    bootstrap_target_model_difference,
    build_signature_matrix,
    ccle_pathway_prediction_benchmark,
    propensity_matched_target_feature_sets,
    target_feature_sets,
)


def test_target_feature_sets_preserve_random_dimensions() -> None:
    candidates = pd.DataFrame(
        {
            "metabolite_a": ["a", "a", "b"],
            "metabolite_b": ["b", "c", "d"],
            "genes": ["G1", "G2", "G3"],
        }
    )
    result = target_feature_sets(
        candidates,
        mapped_metabolites=["a", "b", "c", "d", "e", "f"],
        available_signatures=["G1", "G2", "G3", "G4", "G5"],
        seed=3,
    )
    target = next(item for item in result if item.target == "a")
    assert len(target.network_metabolites) == len(target.random_metabolites)
    assert len(target.network_signatures) == len(target.random_signatures)
    assert len(target.network_interactions) == len(target.random_interactions)


def test_propensity_matched_features_preserve_dimensions_and_exclusions() -> None:
    candidates = pd.DataFrame(
        {
            "metabolite_a": ["a", "a", "b", "c", "d", "e"],
            "metabolite_b": ["b", "c", "d", "e", "f", "f"],
            "genes": ["G1", "G2", "G3", "G4", "G5", "G6"],
        }
    )
    metabolites = list("abcdefghi")
    result = propensity_matched_target_feature_sets(
        candidates,
        mapped_metabolites=metabolites,
        available_signatures=[f"G{index}" for index in range(1, 9)],
        metabolite_coverage=pd.Series(1.0, index=metabolites),
        seed=19,
    )
    target = next(item for item in result if item.target == "a")
    assert len(target.network_metabolites) == len(target.random_metabolites)
    assert len(target.network_signatures) == len(target.random_signatures)
    assert len(target.network_interactions) == len(target.random_interactions)
    assert not set(target.network_metabolites) & set(target.random_metabolites)
    assert not set(target.network_signatures) & set(target.random_signatures)


def test_ccle_benchmark_keeps_lineages_isolated_and_emits_all_models() -> None:
    rng = np.random.default_rng(5)
    lineages = np.repeat([f"l{index}" for index in range(6)], 10)
    latent = rng.normal(size=60)
    metabolomics = pd.DataFrame(
        {
            "a": latent + rng.normal(scale=0.1, size=60),
            "b": latent + rng.normal(scale=0.1, size=60),
            "c": rng.normal(size=60),
            "d": rng.normal(size=60),
        },
        index=[f"s{index}" for index in range(60)],
    )
    expression = pd.DataFrame(
        {"G1": latent + rng.normal(scale=0.2, size=60), "G2": rng.normal(size=60)},
        index=metabolomics.index,
    )
    signatures = build_signature_matrix(expression, ["G1", "G2"], aggregation="mean")
    candidates = pd.DataFrame(
        {"metabolite_a": ["a"], "metabolite_b": ["b"], "genes": ["G1"]}
    )
    feature_sets = target_feature_sets(
        candidates,
        mapped_metabolites=list(metabolomics.columns),
        available_signatures=list(signatures.columns),
        seed=7,
    )
    metrics, targets, features = ccle_pathway_prediction_benchmark(
        metabolomics,
        signatures,
        pd.Series(lineages, index=metabolomics.index),
        feature_sets,
        minimum_target_samples=40,
        minimum_lineage_samples=5,
        folds=3,
        repeats=1,
        ridge_alpha=1.0,
        seed=11,
    )
    assert set(metrics["model"]) == {
        "population_mean",
        "random_factorized_ridge",
        "network_metabolites_ridge",
        "network_additive_ridge",
        "factorized_interaction_ridge",
        "all_metabolites_ridge",
    }
    assert targets["target"].nunique() == 2
    assert features["lineages"].eq(6).all()


def test_target_bootstrap_detects_uniform_improvement() -> None:
    rows = []
    for target in range(20):
        rows.extend(
            [
                {
                    "target": f"m{target}",
                    "model": "reference",
                    "equal_lineage_weighted_rmse_sd": 1.0,
                },
                {
                    "target": f"m{target}",
                    "model": "challenger",
                    "equal_lineage_weighted_rmse_sd": 0.7,
                },
            ]
        )
    result = bootstrap_target_model_difference(
        pd.DataFrame(rows),
        reference="reference",
        challenger="challenger",
        draws=200,
        seed=13,
    )
    assert result["ci_lower"] > 0
