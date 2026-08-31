"""Lineage-isolated prediction from direct HumanGEM and GPR feature sets."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from sklearn.linear_model import Ridge

from ..features.signatures import build_expression_signature, parse_gene_signature
from ..ml.validation import GroupedValidationSpec, repeated_balanced_group_splits


@dataclass(frozen=True)
class TargetFeatureSet:
    """Declared mechanistic and size-matched random features for one metabolite."""

    target: str
    network_metabolites: tuple[str, ...]
    network_signatures: tuple[str, ...]
    network_interactions: tuple[tuple[str, str], ...]
    random_metabolites: tuple[str, ...]
    random_signatures: tuple[str, ...]
    random_interactions: tuple[tuple[str, str], ...]


def build_signature_matrix(
    expression: pd.DataFrame,
    signatures: list[str],
    *,
    aggregation: str,
) -> pd.DataFrame:
    """Build every complete canonical GPR signature exactly once."""
    columns: dict[str, pd.Series] = {}
    for signature in sorted(set(signatures)):
        genes = parse_gene_signature(signature)
        if not set(genes).issubset(expression.columns):
            continue
        columns[signature] = build_expression_signature(
            expression, genes, aggregation=aggregation
        )
    return pd.DataFrame(columns, index=expression.index, dtype=float)


def target_feature_sets(
    candidates: pd.DataFrame,
    *,
    mapped_metabolites: list[str],
    available_signatures: list[str],
    seed: int,
) -> list[TargetFeatureSet]:
    """Create direct-reaction features and dimension-matched random controls."""
    metabolite_set = set(mapped_metabolites)
    signature_set = set(available_signatures)
    edges_by_target: dict[str, set[tuple[str, str]]] = {}
    for row in candidates.to_dict(orient="records"):
        left = str(row["metabolite_a"])
        right = str(row["metabolite_b"])
        signature = str(row["genes"])
        if left not in metabolite_set or right not in metabolite_set:
            continue
        if signature not in signature_set:
            continue
        edges_by_target.setdefault(left, set()).add((right, signature))
        edges_by_target.setdefault(right, set()).add((left, signature))

    results: list[TargetFeatureSet] = []
    for target_index, target in enumerate(sorted(edges_by_target)):
        edges = sorted(edges_by_target[target])
        network_metabolites = tuple(sorted({metabolite for metabolite, _ in edges}))
        network_signatures = tuple(sorted({signature for _, signature in edges}))
        if not network_metabolites or not network_signatures:
            continue
        rng = np.random.default_rng(seed + target_index)
        metabolite_pool = sorted(
            metabolite_set - {target} - set(network_metabolites)
        )
        if len(metabolite_pool) < len(network_metabolites):
            metabolite_pool = sorted(metabolite_set - {target})
        signature_pool = sorted(signature_set - set(network_signatures))
        if len(signature_pool) < len(network_signatures):
            signature_pool = sorted(signature_set)
        random_metabolites = tuple(
            sorted(
                rng.choice(
                    metabolite_pool, size=len(network_metabolites), replace=False
                ).tolist()
            )
        )
        random_signatures = tuple(
            sorted(
                rng.choice(
                    signature_pool, size=len(network_signatures), replace=False
                ).tolist()
            )
        )
        combinations = [
            (metabolite, signature)
            for metabolite in random_metabolites
            for signature in random_signatures
        ]
        selected = rng.choice(
            len(combinations), size=len(edges), replace=False
        )
        random_interactions = tuple(sorted(combinations[int(index)] for index in selected))
        results.append(
            TargetFeatureSet(
                target=target,
                network_metabolites=network_metabolites,
                network_signatures=network_signatures,
                network_interactions=tuple(edges),
                random_metabolites=random_metabolites,
                random_signatures=random_signatures,
                random_interactions=random_interactions,
            )
        )
    return results


def propensity_matched_target_feature_sets(
    candidates: pd.DataFrame,
    *,
    mapped_metabolites: list[str],
    available_signatures: list[str],
    metabolite_coverage: pd.Series,
    seed: int,
    assignment_jitter: float = 0.05,
) -> list[TargetFeatureSet]:
    """Create dimension-, network-degree- and assay-coverage-matched null features.

    Matching uses no target abundance values or target associations. Metabolites are matched on
    direct-network degree and assay coverage; GPR signatures are matched on candidate-network
    degree. A small seeded assignment jitter varies equally plausible matches across null draws.
    """
    metabolite_set = set(mapped_metabolites)
    signature_set = set(available_signatures)
    edges_by_target: dict[str, set[tuple[str, str]]] = {}
    signature_degree: dict[str, int] = {signature: 0 for signature in signature_set}
    for row in candidates.to_dict(orient="records"):
        left = str(row["metabolite_a"])
        right = str(row["metabolite_b"])
        signature = str(row["genes"])
        if left not in metabolite_set or right not in metabolite_set:
            continue
        if signature not in signature_set:
            continue
        edges_by_target.setdefault(left, set()).add((right, signature))
        edges_by_target.setdefault(right, set()).add((left, signature))
        signature_degree[signature] = signature_degree.get(signature, 0) + 1
    metabolite_degree = {
        metabolite: len({neighbor for neighbor, _ in edges_by_target.get(metabolite, set())})
        for metabolite in metabolite_set
    }

    def assignment(
        requested: tuple[str, ...],
        pool: list[str],
        profiles: dict[str, tuple[float, ...]],
        rng: np.random.Generator,
    ) -> tuple[str, ...]:
        if len(pool) < len(requested):
            raise ValueError("Matched random-feature pool is smaller than the requested set.")
        requested_matrix = np.asarray([profiles[value] for value in requested], dtype=float)
        pool_matrix = np.asarray([profiles[value] for value in pool], dtype=float)
        combined = np.vstack([requested_matrix, pool_matrix])
        scale = np.std(combined, axis=0, ddof=0)
        scale = np.where(np.isfinite(scale) & (scale > 1e-12), scale, 1.0)
        costs = np.abs(
            requested_matrix[:, None, :] - pool_matrix[None, :, :]
        ) / scale[None, None, :]
        costs = costs.sum(axis=2)
        costs += rng.uniform(0.0, assignment_jitter, size=costs.shape)
        rows, columns = linear_sum_assignment(costs)
        if len(rows) != len(requested):
            raise ValueError("Matched assignment did not cover every requested feature.")
        return tuple(sorted(pool[int(column)] for column in columns))

    metabolite_profiles = {
        metabolite: (
            float(np.log1p(metabolite_degree.get(metabolite, 0))),
            float(metabolite_coverage.get(metabolite, 0.0)),
        )
        for metabolite in metabolite_set
    }
    signature_profiles = {
        signature: (float(np.log1p(signature_degree.get(signature, 0))),)
        for signature in signature_set
    }
    results: list[TargetFeatureSet] = []
    for target_index, target in enumerate(sorted(edges_by_target)):
        edges = sorted(edges_by_target[target])
        network_metabolites = tuple(sorted({metabolite for metabolite, _ in edges}))
        network_signatures = tuple(sorted({signature for _, signature in edges}))
        if not network_metabolites or not network_signatures:
            continue
        rng = np.random.default_rng(seed + target_index)
        metabolite_pool = sorted(
            metabolite_set - {target} - set(network_metabolites)
        )
        signature_pool = sorted(signature_set - set(network_signatures))
        random_metabolites = assignment(
            network_metabolites, metabolite_pool, metabolite_profiles, rng
        )
        random_signatures = assignment(
            network_signatures, signature_pool, signature_profiles, rng
        )
        combinations = [
            (metabolite, signature)
            for metabolite in random_metabolites
            for signature in random_signatures
        ]
        selected = rng.choice(len(combinations), size=len(edges), replace=False)
        random_interactions = tuple(
            sorted(combinations[int(index)] for index in selected)
        )
        results.append(
            TargetFeatureSet(
                target=target,
                network_metabolites=network_metabolites,
                network_signatures=network_signatures,
                network_interactions=tuple(edges),
                random_metabolites=random_metabolites,
                random_signatures=random_signatures,
                random_interactions=random_interactions,
            )
        )
    return results


def _standardize(
    train: np.ndarray, test: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    finite_count = np.isfinite(train).sum(axis=0)
    center = np.divide(
        np.nansum(train, axis=0),
        finite_count,
        out=np.zeros(train.shape[1], dtype=float),
        where=finite_count > 0,
    )
    centered = np.where(np.isfinite(train), train - center, 0.0)
    scale = np.sqrt(np.mean(centered**2, axis=0))
    scale = np.where(np.isfinite(scale) & (scale > 1e-12), scale, 1.0)
    train_filled = np.where(np.isfinite(train), train, center)
    test_filled = np.where(np.isfinite(test), test, center)
    return (train_filled - center) / scale, (test_filled - center) / scale, center, scale


def _predict_feature_model(
    metabolomics: pd.DataFrame,
    signatures: pd.DataFrame,
    train_index: np.ndarray,
    test_index: np.ndarray,
    train_y: np.ndarray,
    *,
    metabolites: tuple[str, ...],
    expression_signatures: tuple[str, ...] = (),
    interactions: tuple[tuple[str, str], ...] = (),
    alpha: float,
) -> np.ndarray:
    blocks = []
    names: list[str] = []
    if metabolites:
        blocks.append(metabolomics.loc[:, list(metabolites)].to_numpy(dtype=float))
        names.extend(f"M::{value}" for value in metabolites)
    if expression_signatures:
        blocks.append(signatures.loc[:, list(expression_signatures)].to_numpy(dtype=float))
        names.extend(f"T::{value}" for value in expression_signatures)
    if not blocks:
        raise ValueError("A feature model requires at least one base feature.")
    raw = np.column_stack(blocks)
    train_base, test_base, _, _ = _standardize(raw[train_index], raw[test_index])
    name_index = {name: index for index, name in enumerate(names)}
    train_columns = [train_base]
    test_columns = [test_base]
    if interactions:
        train_products = []
        test_products = []
        for metabolite, signature in interactions:
            metabolite_index = name_index[f"M::{metabolite}"]
            signature_index = name_index[f"T::{signature}"]
            train_products.append(
                train_base[:, metabolite_index] * train_base[:, signature_index]
            )
            test_products.append(
                test_base[:, metabolite_index] * test_base[:, signature_index]
            )
        train_interactions = np.column_stack(train_products)
        test_interactions = np.column_stack(test_products)
        train_interactions, test_interactions, _, _ = _standardize(
            train_interactions, test_interactions
        )
        train_columns.append(train_interactions)
        test_columns.append(test_interactions)
    model = Ridge(alpha=alpha, solver="lsqr")
    model.fit(np.column_stack(train_columns), train_y)
    return np.asarray(model.predict(np.column_stack(test_columns)), dtype=float)


def ccle_pathway_prediction_benchmark(
    metabolomics: pd.DataFrame,
    signatures: pd.DataFrame,
    lineages: pd.Series,
    feature_sets: list[TargetFeatureSet],
    *,
    minimum_target_samples: int,
    minimum_lineage_samples: int,
    folds: int,
    repeats: int,
    ridge_alpha: float,
    seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Compare compact mechanism-aware features in lineage-isolated folds."""
    model_names = (
        "population_mean",
        "random_factorized_ridge",
        "network_metabolites_ridge",
        "network_additive_ridge",
        "factorized_interaction_ridge",
        "all_metabolites_ridge",
    )
    all_metabolites = tuple(sorted(metabolomics.columns))
    prediction_records: list[dict[str, object]] = []
    feature_records: list[dict[str, object]] = []

    for feature_set in feature_sets:
        target = feature_set.target
        target_values = metabolomics[target]
        eligible = target_values.notna() & lineages.notna()
        lineage_counts = lineages.loc[eligible].astype(str).value_counts()
        retained_lineages = set(
            lineage_counts.loc[lineage_counts >= minimum_lineage_samples].index
        )
        eligible &= lineages.astype(str).isin(retained_lineages)
        indexes = np.flatnonzero(eligible.to_numpy())
        if len(indexes) < minimum_target_samples or len(retained_lineages) < folds:
            continue
        local_metabolomics = metabolomics.iloc[indexes]
        local_signatures = signatures.iloc[indexes]
        local_lineages = lineages.iloc[indexes].astype(str)
        target_y = target_values.iloc[indexes].to_numpy(dtype=float)
        splits = repeated_balanced_group_splits(
            local_lineages,
            GroupedValidationSpec(outer_folds=folds, repeats=repeats, seed=seed),
        )
        feature_records.append(
            {
                "target": target,
                "samples": len(indexes),
                "lineages": len(retained_lineages),
                "network_metabolites": ";".join(feature_set.network_metabolites),
                "network_signatures": ";".join(feature_set.network_signatures),
                "network_interactions": ";".join(
                    f"{metabolite}*{signature}"
                    for metabolite, signature in feature_set.network_interactions
                ),
                "random_metabolites": ";".join(feature_set.random_metabolites),
                "random_signatures": ";".join(feature_set.random_signatures),
                "random_interactions": ";".join(
                    f"{metabolite}*{signature}"
                    for metabolite, signature in feature_set.random_interactions
                ),
                "network_metabolite_count": len(feature_set.network_metabolites),
                "network_signature_count": len(feature_set.network_signatures),
                "network_interaction_count": len(feature_set.network_interactions),
            }
        )
        for train, test, repeat, fold in splits:
            train_center = float(np.mean(target_y[train]))
            train_scale = float(np.std(target_y[train], ddof=0))
            if not np.isfinite(train_scale) or train_scale <= 1e-12:
                continue
            train_y = (target_y[train] - train_center) / train_scale
            test_y = (target_y[test] - train_center) / train_scale
            predictions = {
                "population_mean": np.zeros(len(test), dtype=float),
                "random_factorized_ridge": _predict_feature_model(
                    local_metabolomics,
                    local_signatures,
                    train,
                    test,
                    train_y,
                    metabolites=feature_set.random_metabolites,
                    expression_signatures=feature_set.random_signatures,
                    interactions=feature_set.random_interactions,
                    alpha=ridge_alpha,
                ),
                "network_metabolites_ridge": _predict_feature_model(
                    local_metabolomics,
                    local_signatures,
                    train,
                    test,
                    train_y,
                    metabolites=feature_set.network_metabolites,
                    alpha=ridge_alpha,
                ),
                "network_additive_ridge": _predict_feature_model(
                    local_metabolomics,
                    local_signatures,
                    train,
                    test,
                    train_y,
                    metabolites=feature_set.network_metabolites,
                    expression_signatures=feature_set.network_signatures,
                    alpha=ridge_alpha,
                ),
                "factorized_interaction_ridge": _predict_feature_model(
                    local_metabolomics,
                    local_signatures,
                    train,
                    test,
                    train_y,
                    metabolites=feature_set.network_metabolites,
                    expression_signatures=feature_set.network_signatures,
                    interactions=feature_set.network_interactions,
                    alpha=ridge_alpha,
                ),
                "all_metabolites_ridge": _predict_feature_model(
                    local_metabolomics,
                    local_signatures,
                    train,
                    test,
                    train_y,
                    metabolites=tuple(value for value in all_metabolites if value != target),
                    alpha=ridge_alpha,
                ),
            }
            test_ids = local_metabolomics.index[test]
            test_lineages = local_lineages.iloc[test].to_numpy()
            for model_name in model_names:
                prediction = predictions[model_name]
                for row, sample_id in enumerate(test_ids):
                    prediction_records.append(
                        {
                            "target": target,
                            "sample_id": str(sample_id),
                            "lineage": str(test_lineages[row]),
                            "repeat": int(repeat),
                            "fold": int(fold),
                            "model": model_name,
                            "observed_sd": float(test_y[row]),
                            "predicted_sd": float(prediction[row]),
                            "squared_error": float((test_y[row] - prediction[row]) ** 2),
                            "truth_energy": float(test_y[row] ** 2),
                        }
                    )

    predictions = pd.DataFrame.from_records(prediction_records)
    if predictions.empty:
        raise ValueError("No CCLE target passed the configured eligibility gates.")
    lineage_mse = (
        predictions.groupby(["target", "model", "lineage"], as_index=False)
        .agg(mean_squared_error=("squared_error", "mean"))
    )
    equal_lineage = (
        lineage_mse.groupby(["target", "model"], as_index=False)
        .agg(equal_lineage_mse=("mean_squared_error", "mean"))
    )
    sums = (
        predictions.groupby(["target", "model"], as_index=False)
        .agg(
            rows=("squared_error", "size"),
            squared_error_sum=("squared_error", "sum"),
            truth_energy_sum=("truth_energy", "sum"),
        )
        .merge(equal_lineage, on=["target", "model"], validate="one_to_one")
    )
    sums["row_weighted_rmse_sd"] = np.sqrt(sums["squared_error_sum"] / sums["rows"])
    sums["equal_lineage_weighted_rmse_sd"] = np.sqrt(sums["equal_lineage_mse"])
    sums["row_weighted_r2"] = 1.0 - sums["squared_error_sum"] / sums["truth_energy_sum"]
    target_metrics = sums.drop(columns="equal_lineage_mse")
    model_metrics = (
        target_metrics.groupby("model", as_index=False)
        .agg(
            targets=("target", "nunique"),
            mean_equal_lineage_rmse_sd=("equal_lineage_weighted_rmse_sd", "mean"),
            median_target_r2=("row_weighted_r2", "median"),
            positive_r2_targets=("row_weighted_r2", lambda value: int((value > 0).sum())),
        )
    )
    return model_metrics, target_metrics, pd.DataFrame.from_records(feature_records)


def bootstrap_target_model_difference(
    target_metrics: pd.DataFrame,
    *,
    reference: str,
    challenger: str,
    draws: int,
    seed: int,
) -> dict[str, float | int | str]:
    """Bootstrap paired target-level equal-lineage RMSE differences."""
    pivot = target_metrics.pivot(
        index="target", columns="model", values="equal_lineage_weighted_rmse_sd"
    )
    required = [reference, challenger]
    if any(model not in pivot for model in required):
        raise ValueError("Requested target comparison model is absent.")
    values = pivot[required].dropna().to_numpy(dtype=float)
    rng = np.random.default_rng(seed)
    indexes = rng.integers(0, len(values), size=(draws, len(values)))
    differences = values[indexes, 0].mean(axis=1) - values[indexes, 1].mean(axis=1)
    point = float(values[:, 0].mean() - values[:, 1].mean())
    return {
        "reference": reference,
        "challenger": challenger,
        "targets": len(values),
        "draws": draws,
        "rmse_improvement_sd": point,
        "ci_lower": float(np.quantile(differences, 0.025)),
        "ci_upper": float(np.quantile(differences, 0.975)),
    }
