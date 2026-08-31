"""Leakage-safe benchmarks for compact biochemical pathway representations."""

from __future__ import annotations

import re
from collections.abc import Mapping

import numpy as np
import pandas as pd
from sklearn.linear_model import Ridge

from ..ml.validation import GroupedValidationSpec, repeated_balanced_group_splits

_LIPID_FAMILY_ALIASES = {
    "lysoPC": "LPC",
    "lysoPE": "LPE",
    "TG": "TAG",
    "DG": "DAG",
}
_LIPID_FAMILY_PATTERN = re.compile(
    r"^(lysoPC|lysoPE|HexCer|LacCer|LPC|LPE|TAG|DAG|Cer|SM|PC|PE|PG|PI|PS|CL|CE|TG|DG)"
    r"(?=\(|\s|\d)",
    flags=re.IGNORECASE,
)


def lipid_family(feature: str) -> str:
    """Return the declared ST002081 lipid family from an assay feature name."""
    name = str(feature)
    known = _LIPID_FAMILY_PATTERN.match(name)
    if known:
        observed = known.group(1)
        canonical = next(
            (
                family
                for family in (
                    "lysoPC",
                    "lysoPE",
                    "HexCer",
                    "LacCer",
                    "LPC",
                    "LPE",
                    "TAG",
                    "DAG",
                    "Cer",
                    "SM",
                    "PC",
                    "PE",
                    "PG",
                    "PI",
                    "PS",
                    "CL",
                    "CE",
                    "TG",
                    "DG",
                )
                if family.lower() == observed.lower()
            ),
            observed,
        )
        return _LIPID_FAMILY_ALIASES.get(canonical, canonical)
    match = re.match(r"^([^\(]+)", name)
    return match.group(1).strip() if match else name


def lipid_structure_descriptors(feature: str) -> tuple[str, ...]:
    """Parse outcome-independent lipid headgroup and chain descriptors."""
    name = str(feature)
    family = lipid_family(name)
    descriptors = {f"family:{family}"}
    if family == "TAG":
        total_match = re.match(r"^(?:TAG|TG)\s*(\d+):(\d+)", name)
        chain_matches = re.findall(r"FA(\d+):(\d+)", name)
        if total_match:
            total_carbon, total_unsaturation = map(int, total_match.groups())
            descriptors.update(
                {
                    f"total_carbon:{total_carbon}",
                    f"total_unsaturation:{total_unsaturation}",
                    f"family_unsaturation:{family}:{total_unsaturation}",
                }
            )
    else:
        inside = re.search(r"\(([^)]+)\)", name)
        if inside:
            chain_matches = re.findall(r"(\d+):(\d+)", inside.group(1))
        else:
            composition = re.match(
                r"^[A-Za-z]+(?:Cer)?\s+[dmtOP-]*(\d+):(\d+)", name
            )
            chain_matches = [composition.groups()] if composition else []
        if chain_matches:
            total_carbon = sum(int(carbon) for carbon, _ in chain_matches)
            total_unsaturation = sum(int(double) for _, double in chain_matches)
            descriptors.update(
                {
                    f"total_carbon:{total_carbon}",
                    f"total_unsaturation:{total_unsaturation}",
                    f"family_unsaturation:{family}:{total_unsaturation}",
                }
            )
    for carbon, double in chain_matches:
        descriptors.update(
            {
                f"acyl:{carbon}:{double}",
                f"chain_carbon:{carbon}",
                f"chain_unsaturation:{double}",
            }
        )
    return tuple(sorted(descriptors))


def structural_descriptor_incidence(
    features: pd.Index,
    *,
    minimum_features: int,
    maximum_feature_fraction: float,
) -> pd.DataFrame:
    """Return a filtered feature-by-structural-descriptor incidence matrix."""
    records = {
        str(feature): set(lipid_structure_descriptors(str(feature))) for feature in features
    }
    descriptors = sorted(set().union(*records.values()))
    incidence = pd.DataFrame(
        {
            descriptor: [descriptor in records[str(feature)] for feature in features]
            for descriptor in descriptors
        },
        index=pd.Index([str(feature) for feature in features], name="feature"),
        dtype=bool,
    )
    counts = incidence.sum(axis=0)
    keep = (counts >= minimum_features) & (
        counts <= maximum_feature_fraction * len(incidence)
    )
    return incidence.loc[:, keep].copy()


def degree_preserving_descriptor_null(
    incidence: pd.DataFrame, *, swaps_per_edge: int, seed: int
) -> pd.DataFrame:
    """Randomize a bipartite incidence graph while preserving both degree sequences."""
    matrix = incidence.to_numpy(dtype=bool, copy=True)
    rows, columns = np.nonzero(matrix)
    if len(rows) < 2:
        raise ValueError("Descriptor incidence has fewer than two edges.")
    rng = np.random.default_rng(seed)
    attempts = swaps_per_edge * len(rows)
    for _ in range(attempts):
        left, right = rng.integers(0, len(rows), size=2)
        row_left, column_left = int(rows[left]), int(columns[left])
        row_right, column_right = int(rows[right]), int(columns[right])
        if row_left == row_right or column_left == column_right:
            continue
        if matrix[row_left, column_right] or matrix[row_right, column_left]:
            continue
        matrix[row_left, column_left] = False
        matrix[row_right, column_right] = False
        matrix[row_left, column_right] = True
        matrix[row_right, column_left] = True
        columns[left], columns[right] = column_right, column_left
    return pd.DataFrame(matrix, index=incidence.index, columns=incidence.columns)


def _descriptor_scores(
    standardized: np.ndarray,
    visible_features: list[str],
    incidence: pd.DataFrame,
) -> np.ndarray:
    columns = []
    for descriptor in incidence.columns:
        members = incidence.index[incidence[descriptor]].intersection(visible_features)
        indices = [visible_features.index(str(member)) for member in members]
        if not indices:
            raise ValueError(f"Descriptor {descriptor!r} has no visible members.")
        columns.append(np.median(standardized[:, indices], axis=1))
    return np.column_stack(columns)


def masked_structural_descriptor_benchmark(
    values: pd.DataFrame,
    subjects: pd.Series,
    masks: pd.DataFrame,
    incidence: pd.DataFrame,
    *,
    folds: int,
    ridge_alpha: float,
    split_seed: int,
    null_seed: int,
    swaps_per_edge: int,
    priority_count: int,
    true_priority_fraction: float,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Compare structural descriptors with a bipartite degree-preserving null."""
    feature_order = list(incidence.index)
    matrix = values.loc[:, feature_order].astype(float)
    if matrix.isna().any().any():
        raise ValueError("The structural benchmark requires complete eligible features.")
    aligned_subjects = subjects.reindex(matrix.index).astype(str)
    splits = repeated_balanced_group_splits(
        aligned_subjects,
        GroupedValidationSpec(outer_folds=folds, repeats=1, seed=split_seed),
    )
    model_names = (
        "population_mean",
        "degree_preserving_random_structural_ridge",
        "structural_descriptor_ridge",
    )
    sample_records: list[dict[str, object]] = []
    null_records: list[dict[str, object]] = []
    for mask_id in sorted(masks["mask"].unique()):
        hidden = sorted(masks.loc[masks["mask"].eq(mask_id), "feature"].astype(str))
        visible = sorted(set(feature_order) - set(hidden))
        visible_incidence = incidence.loc[visible]
        usable = visible_incidence.columns[visible_incidence.sum(axis=0).ge(2)]
        visible_incidence = visible_incidence.loc[:, usable]
        null_incidence = degree_preserving_descriptor_null(
            visible_incidence,
            swaps_per_edge=swaps_per_edge,
            seed=null_seed + int(mask_id),
        )
        null_records.append(
            {
                "mask": int(mask_id),
                "visible_features": len(visible),
                "hidden_features": len(hidden),
                "descriptors": len(usable),
                "incidence_edges": int(visible_incidence.to_numpy().sum()),
                "row_degrees_preserved": bool(
                    np.array_equal(
                        visible_incidence.sum(axis=1).to_numpy(),
                        null_incidence.sum(axis=1).to_numpy(),
                    )
                ),
                "column_degrees_preserved": bool(
                    np.array_equal(
                        visible_incidence.sum(axis=0).to_numpy(),
                        null_incidence.sum(axis=0).to_numpy(),
                    )
                ),
            }
        )
        for train_index, test_index, repeat, fold in splits:
            train_visible = matrix.iloc[train_index].loc[:, visible].to_numpy(dtype=float)
            test_visible = matrix.iloc[test_index].loc[:, visible].to_numpy(dtype=float)
            train_visible_z, test_visible_z, _, _ = _standardize(train_visible, test_visible)
            train_y = matrix.iloc[train_index].loc[:, hidden].to_numpy(dtype=float)
            test_y = matrix.iloc[test_index].loc[:, hidden].to_numpy(dtype=float)
            train_y_z, test_y_z, _, _ = _standardize(train_y, test_y)
            train_structural = _descriptor_scores(
                train_visible_z, visible, visible_incidence
            )
            test_structural = _descriptor_scores(test_visible_z, visible, visible_incidence)
            train_random = _descriptor_scores(train_visible_z, visible, null_incidence)
            test_random = _descriptor_scores(test_visible_z, visible, null_incidence)
            predictions = {
                "population_mean": np.zeros_like(test_y_z),
                "degree_preserving_random_structural_ridge": _ridge_prediction(
                    train_random, test_random, train_y_z, alpha=ridge_alpha
                ),
                "structural_descriptor_ridge": _ridge_prediction(
                    train_structural, test_structural, train_y_z, alpha=ridge_alpha
                ),
            }
            test_ids = matrix.index[test_index]
            test_subjects = aligned_subjects.iloc[test_index].to_numpy()
            for model_name in model_names:
                prediction = predictions[model_name]
                squared = (test_y_z - prediction) ** 2
                precision, ndcg = _priority_metrics(
                    test_y_z,
                    prediction,
                    priority_count=priority_count,
                    true_priority_fraction=true_priority_fraction,
                )
                for row, sample_id in enumerate(test_ids):
                    sample_records.append(
                        {
                            "sample_id": str(sample_id),
                            "subject_id": str(test_subjects[row]),
                            "mask": int(mask_id),
                            "repeat": int(repeat),
                            "fold": int(fold),
                            "model": model_name,
                            "hidden_features": len(hidden),
                            "mean_squared_error": float(squared[row].mean()),
                            "truth_energy": float((test_y_z[row] ** 2).mean()),
                            "precision_at_k": float(precision[row]),
                            "ndcg_at_k": float(ndcg[row]),
                        }
                    )
    sample_metrics = pd.DataFrame.from_records(sample_records)
    model_records = []
    for model_name, frame in sample_metrics.groupby("model", sort=True):
        pairs = frame["hidden_features"].sum()
        sse = float((frame["mean_squared_error"] * frame["hidden_features"]).sum())
        sst = float((frame["truth_energy"] * frame["hidden_features"]).sum())
        subject_mse = frame.groupby("subject_id")["mean_squared_error"].mean()
        model_records.append(
            {
                "model": model_name,
                "samples": frame["sample_id"].nunique(),
                "subjects": frame["subject_id"].nunique(),
                "prediction_pairs": int(pairs),
                "row_weighted_rmse_sd": float(np.sqrt(sse / pairs)),
                "equal_subject_weighted_rmse_sd": float(np.sqrt(subject_mse.mean())),
                "pooled_r2": float(1.0 - sse / sst),
                "precision_at_k": float(frame["precision_at_k"].mean()),
                "ndcg_at_k": float(frame["ndcg_at_k"].mean()),
            }
        )
    return (
        pd.DataFrame.from_records(model_records),
        sample_metrics,
        pd.DataFrame.from_records(null_records),
    )


def eligible_lipid_families(
    features: pd.Index, *, minimum_family_features: int
) -> pd.Series:
    """Map features to deterministic families and remove undersized families."""
    mapping = pd.Series(
        {str(feature): lipid_family(str(feature)) for feature in features},
        name="family",
        dtype="object",
    )
    counts = mapping.value_counts()
    eligible = set(counts.loc[counts >= minimum_family_features].index)
    return mapping.loc[mapping.isin(eligible)].sort_index()


def disjoint_family_masks(
    family_by_feature: pd.Series, *, masks: int, seed: int
) -> pd.DataFrame:
    """Assign every feature to exactly one seeded, family-stratified mask."""
    if masks < 2:
        raise ValueError("At least two masks are required.")
    rng = np.random.default_rng(seed)
    records: list[dict[str, object]] = []
    for family, members in family_by_feature.groupby(family_by_feature, sort=True):
        ordered = np.asarray(sorted(members.index), dtype=object)
        shuffled = ordered[rng.permutation(len(ordered))]
        for index, feature in enumerate(shuffled):
            records.append(
                {"feature": str(feature), "family": str(family), "mask": index % masks}
            )
    result = pd.DataFrame.from_records(records).sort_values(["mask", "family", "feature"])
    if result["feature"].duplicated().any():
        raise ValueError("A feature was assigned to more than one mask.")
    return result.reset_index(drop=True)


def _standardize(
    train: np.ndarray, test: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    center = np.nanmean(train, axis=0)
    scale = np.nanstd(train, axis=0, ddof=0)
    scale = np.where(np.isfinite(scale) & (scale > 1e-12), scale, 1.0)
    train_filled = np.where(np.isfinite(train), train, center)
    test_filled = np.where(np.isfinite(test), test, center)
    return (train_filled - center) / scale, (test_filled - center) / scale, center, scale


def _score_matrix(
    standardized: np.ndarray,
    features: list[str],
    group_by_feature: Mapping[str, str],
    group_order: list[str],
) -> np.ndarray:
    columns = []
    feature_index = {feature: index for index, feature in enumerate(features)}
    for group in group_order:
        indices = [
            feature_index[feature]
            for feature in features
            if group_by_feature[feature] == group
        ]
        if not indices:
            raise ValueError(f"Score group {group!r} has no visible features.")
        columns.append(np.median(standardized[:, indices], axis=1))
    return np.column_stack(columns)


def _random_visible_groups(
    visible: list[str], family_by_feature: pd.Series, *, seed: int
) -> dict[str, str]:
    """Permute visible labels while preserving every visible family size."""
    labels = np.asarray([str(family_by_feature.loc[feature]) for feature in visible], dtype=object)
    rng = np.random.default_rng(seed)
    shuffled = labels[rng.permutation(len(labels))]
    return dict(zip(visible, shuffled, strict=True))


def _ridge_prediction(
    train_x: np.ndarray,
    test_x: np.ndarray,
    train_y_standardized: np.ndarray,
    *,
    alpha: float,
) -> np.ndarray:
    train_scaled, test_scaled, _, _ = _standardize(train_x, test_x)
    model = Ridge(alpha=alpha, solver="lsqr")
    model.fit(train_scaled, train_y_standardized)
    return np.asarray(model.predict(test_scaled), dtype=float)


def _priority_metrics(
    truth: np.ndarray,
    prediction: np.ndarray,
    *,
    priority_count: int,
    true_priority_fraction: float,
) -> tuple[np.ndarray, np.ndarray]:
    features = truth.shape[1]
    k = min(priority_count, features)
    truth_count = max(k, int(np.ceil(features * true_priority_fraction)))
    precision = np.empty(len(truth), dtype=float)
    ndcg = np.empty(len(truth), dtype=float)
    discounts = 1.0 / np.log2(np.arange(k) + 2.0)
    for row in range(len(truth)):
        relevance = np.abs(truth[row])
        predicted_order = np.argsort(-np.abs(prediction[row]), kind="stable")[:k]
        truth_set = set(np.argsort(-relevance, kind="stable")[:truth_count])
        precision[row] = len(set(predicted_order) & truth_set) / k
        dcg = float(np.sum(relevance[predicted_order] * discounts))
        ideal_order = np.argsort(-relevance, kind="stable")[:k]
        ideal = float(np.sum(relevance[ideal_order] * discounts))
        ndcg[row] = dcg / ideal if ideal > 0 else 0.0
    return precision, ndcg


def masked_panel_benchmark(
    values: pd.DataFrame,
    subjects: pd.Series,
    family_by_feature: pd.Series,
    masks: pd.DataFrame,
    *,
    folds: int,
    ridge_alpha: float,
    random_group_seed: int,
    priority_count: int,
    true_priority_fraction: float,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Evaluate disjoint hidden panels with participant-isolated outer folds."""
    feature_order = list(family_by_feature.index)
    matrix = values.loc[:, feature_order].astype(float)
    if matrix.isna().any().any():
        raise ValueError("The masked-panel benchmark requires complete eligible features.")
    aligned_subjects = subjects.reindex(matrix.index).astype(str)
    split_spec = GroupedValidationSpec(outer_folds=folds, repeats=1, seed=random_group_seed)
    splits = repeated_balanced_group_splits(aligned_subjects, split_spec)
    model_names = (
        "population_mean",
        "random_group_score_ridge",
        "family_score_ridge",
        "all_visible_ridge",
    )
    sample_records: list[dict[str, object]] = []
    target_records: list[dict[str, object]] = []

    for mask_id in sorted(masks["mask"].unique()):
        hidden = sorted(masks.loc[masks["mask"].eq(mask_id), "feature"].astype(str))
        visible = sorted(set(feature_order) - set(hidden))
        group_order = sorted(set(family_by_feature.loc[visible].astype(str)))
        random_groups = _random_visible_groups(
            visible,
            family_by_feature,
            seed=random_group_seed + int(mask_id),
        )
        family_groups = {feature: str(family_by_feature.loc[feature]) for feature in visible}

        for train_index, test_index, repeat, fold in splits:
            train_visible = matrix.iloc[train_index].loc[:, visible].to_numpy(dtype=float)
            test_visible = matrix.iloc[test_index].loc[:, visible].to_numpy(dtype=float)
            train_visible_z, test_visible_z, _, _ = _standardize(train_visible, test_visible)
            train_y = matrix.iloc[train_index].loc[:, hidden].to_numpy(dtype=float)
            test_y = matrix.iloc[test_index].loc[:, hidden].to_numpy(dtype=float)
            train_y_z, test_y_z, _, _ = _standardize(train_y, test_y)

            train_family = _score_matrix(
                train_visible_z, visible, family_groups, group_order
            )
            test_family = _score_matrix(test_visible_z, visible, family_groups, group_order)
            train_random = _score_matrix(
                train_visible_z, visible, random_groups, group_order
            )
            test_random = _score_matrix(test_visible_z, visible, random_groups, group_order)
            predictions = {
                "population_mean": np.zeros_like(test_y_z),
                "random_group_score_ridge": _ridge_prediction(
                    train_random,
                    test_random,
                    train_y_z,
                    alpha=ridge_alpha,
                ),
                "family_score_ridge": _ridge_prediction(
                    train_family,
                    test_family,
                    train_y_z,
                    alpha=ridge_alpha,
                ),
                "all_visible_ridge": _ridge_prediction(
                    train_visible_z,
                    test_visible_z,
                    train_y_z,
                    alpha=ridge_alpha,
                ),
            }

            test_ids = matrix.index[test_index]
            test_subjects = aligned_subjects.iloc[test_index].to_numpy()
            for model_name in model_names:
                prediction = predictions[model_name]
                squared = (test_y_z - prediction) ** 2
                precision, ndcg = _priority_metrics(
                    test_y_z,
                    prediction,
                    priority_count=priority_count,
                    true_priority_fraction=true_priority_fraction,
                )
                for row, sample_id in enumerate(test_ids):
                    sample_records.append(
                        {
                            "sample_id": str(sample_id),
                            "subject_id": str(test_subjects[row]),
                            "mask": int(mask_id),
                            "repeat": int(repeat),
                            "fold": int(fold),
                            "model": model_name,
                            "hidden_features": len(hidden),
                            "mean_squared_error": float(squared[row].mean()),
                            "truth_energy": float((test_y_z[row] ** 2).mean()),
                            "precision_at_k": float(precision[row]),
                            "ndcg_at_k": float(ndcg[row]),
                        }
                    )
                for column, feature in enumerate(hidden):
                    target_records.append(
                        {
                            "feature": feature,
                            "family": str(family_by_feature.loc[feature]),
                            "mask": int(mask_id),
                            "repeat": int(repeat),
                            "fold": int(fold),
                            "model": model_name,
                            "n": len(test_index),
                            "squared_error_sum": float(squared[:, column].sum()),
                            "truth_energy_sum": float((test_y_z[:, column] ** 2).sum()),
                        }
                    )

    sample_metrics = pd.DataFrame.from_records(sample_records)
    target_parts = pd.DataFrame.from_records(target_records)
    target_metrics = (
        target_parts.groupby(["feature", "family", "mask", "model"], as_index=False)
        .agg(
            n=("n", "sum"),
            squared_error_sum=("squared_error_sum", "sum"),
            truth_energy_sum=("truth_energy_sum", "sum"),
        )
    )
    target_metrics["rmse_sd"] = np.sqrt(
        target_metrics["squared_error_sum"] / target_metrics["n"]
    )
    target_metrics["r2"] = 1.0 - (
        target_metrics["squared_error_sum"] / target_metrics["truth_energy_sum"]
    )

    metric_records = []
    for model_name, frame in sample_metrics.groupby("model", sort=True):
        pairs = frame["hidden_features"].sum()
        sse = float((frame["mean_squared_error"] * frame["hidden_features"]).sum())
        sst = float((frame["truth_energy"] * frame["hidden_features"]).sum())
        subject_mse = frame.groupby("subject_id")["mean_squared_error"].mean()
        subject_precision = frame.groupby("subject_id")["precision_at_k"].mean()
        metric_records.append(
            {
                "model": model_name,
                "samples": frame["sample_id"].nunique(),
                "subjects": frame["subject_id"].nunique(),
                "prediction_pairs": int(pairs),
                "row_weighted_rmse_sd": float(np.sqrt(sse / pairs)),
                "equal_subject_weighted_rmse_sd": float(np.sqrt(subject_mse.mean())),
                "pooled_r2": float(1.0 - sse / sst),
                "precision_at_k": float(frame["precision_at_k"].mean()),
                "equal_subject_precision_at_k": float(subject_precision.mean()),
                "ndcg_at_k": float(frame["ndcg_at_k"].mean()),
            }
        )
    model_metrics = pd.DataFrame.from_records(metric_records)
    return model_metrics, sample_metrics, target_metrics


def bootstrap_subject_model_difference(
    sample_metrics: pd.DataFrame,
    *,
    reference: str,
    challenger: str,
    draws: int,
    seed: int,
) -> dict[str, float | int | str]:
    """Bootstrap participant-paired RMSE and precision differences."""
    subject_model = (
        sample_metrics.groupby(["subject_id", "model"], as_index=False)
        .agg(
            mean_squared_error=("mean_squared_error", "mean"),
            precision_at_k=("precision_at_k", "mean"),
        )
    )
    mse = subject_model.pivot(index="subject_id", columns="model", values="mean_squared_error")
    precision = subject_model.pivot(
        index="subject_id", columns="model", values="precision_at_k"
    )
    required = [reference, challenger]
    if any(model not in mse or model not in precision for model in required):
        raise ValueError("Requested comparison model is absent.")
    shared = mse[required].dropna().index.intersection(precision[required].dropna().index)
    mse_values = mse.loc[shared, required].to_numpy(dtype=float)
    precision_values = precision.loc[shared, required].to_numpy(dtype=float)
    rng = np.random.default_rng(seed)
    indices = rng.integers(0, len(shared), size=(draws, len(shared)))
    sampled_mse = mse_values[indices]
    rmse_difference = np.sqrt(sampled_mse[:, :, 0].mean(axis=1)) - np.sqrt(
        sampled_mse[:, :, 1].mean(axis=1)
    )
    sampled_precision = precision_values[indices]
    precision_difference = (
        sampled_precision[:, :, 1].mean(axis=1)
        - sampled_precision[:, :, 0].mean(axis=1)
    )
    point_rmse = float(np.sqrt(mse_values[:, 0].mean()) - np.sqrt(mse_values[:, 1].mean()))
    point_precision = float(precision_values[:, 1].mean() - precision_values[:, 0].mean())
    return {
        "reference": reference,
        "challenger": challenger,
        "subjects": len(shared),
        "draws": draws,
        "rmse_improvement_sd": point_rmse,
        "rmse_ci_lower": float(np.quantile(rmse_difference, 0.025)),
        "rmse_ci_upper": float(np.quantile(rmse_difference, 0.975)),
        "precision_improvement": point_precision,
        "precision_ci_lower": float(np.quantile(precision_difference, 0.025)),
        "precision_ci_upper": float(np.quantile(precision_difference, 0.975)),
    }
