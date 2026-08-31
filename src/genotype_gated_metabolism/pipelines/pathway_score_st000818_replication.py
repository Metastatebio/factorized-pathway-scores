"""Run the locked ST000818 external structural pathway-score replication."""

from __future__ import annotations

import argparse
import hashlib
import inspect
import json
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import yaml
from joblib import Parallel, delayed, parallel_backend

from ..analysis.intervention_registry import parse_factor_string
from ..analysis.pathway_scores import (
    bootstrap_subject_model_difference,
    disjoint_family_masks,
    eligible_lipid_families,
    masked_panel_benchmark,
    masked_structural_descriptor_benchmark,
    structural_descriptor_incidence,
)
from ..datasets.metabolomics_workbench import (
    DEFAULT_BASE_URL,
    MetabolomicsWorkbenchClient,
)
from ..reporting import dataframe_to_markdown
from .pathway_score_st002081_structural_null_sensitivity import _run_null_seed


def _resolve(config_path: Path, value: str) -> Path:
    path = Path(value)
    return path if path.is_absolute() else (config_path.parent / path).resolve()


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _cache_path(cache_dir: Path, url: str) -> Path:
    return cache_dir / f"{hashlib.sha256(url.encode()).hexdigest()}.txt"


def implementation_sha256() -> str:
    digest = hashlib.sha256()
    for implementation in (
        run,
        masked_panel_benchmark,
        masked_structural_descriptor_benchmark,
    ):
        path = Path(inspect.getfile(implementation))
        digest.update(path.name.encode())
        digest.update(_sha256(path).encode())
    return digest.hexdigest()


def _annotate_groups(factors: pd.DataFrame, factor_name: str) -> pd.Series:
    parsed = factors["factors"].map(parse_factor_string)
    groups = parsed.map(lambda values: values.get(factor_name))
    return pd.Series(
        groups.to_numpy(),
        index=factors["local_sample_id"].astype(str),
        name="validation_group",
        dtype="object",
    )


def adjudicate_external_replication(
    null_results: pd.DataFrame,
    *,
    eligible_features: int,
    validation_groups: int,
    minimum_features: int,
    minimum_groups: int,
    threshold: float,
) -> tuple[str, dict[str, bool]]:
    """Apply the locked all-null-realizations external replication gates."""
    gates = {
        "minimum_features": eligible_features >= minimum_features,
        "minimum_validation_groups": validation_groups >= minimum_groups,
        "all_degrees_preserved": bool(
            null_results["row_degrees_preserved"].all()
            and null_results["column_degrees_preserved"].all()
        ),
        "rmse_ci_passes_all_nulls": bool(
            null_results["rmse_ci_lower"].gt(threshold).all()
        ),
        "precision_ci_passes_all_nulls": bool(
            null_results["precision_ci_lower"].gt(threshold).all()
        ),
    }
    common = (
        gates["minimum_features"]
        and gates["minimum_validation_groups"]
        and gates["all_degrees_preserved"]
    )
    if common and gates["rmse_ci_passes_all_nulls"] and gates[
        "precision_ci_passes_all_nulls"
    ]:
        decision = "EXTERNAL_STRUCTURAL_REPLICATION"
    elif common and gates["rmse_ci_passes_all_nulls"]:
        decision = "EXTERNAL_STRUCTURAL_RECONSTRUCTION_ONLY"
    else:
        decision = "EXTERNAL_STRUCTURAL_REPLICATION_FAILED"
    return decision, gates


def _render_report(
    manifest: dict[str, Any],
    metrics: pd.DataFrame,
    comparisons: pd.DataFrame,
    null_results: pd.DataFrame,
) -> str:
    return "\n".join(
        [
            "# ST000818 external structural pathway-score replication",
            "",
            f"**Decision:** `{manifest['decision']}`",
            "",
            "## Cohort",
            "",
            (
                f"The locked analysis retained {manifest['eligible_features']} complete lipids "
                f"and {manifest['descriptors']} descriptors across {manifest['samples']} "
                f"participants from {manifest['validation_groups']} population categories. "
                "Every outer fold excluded complete categories."
            ),
            "",
            "## Model performance",
            "",
            dataframe_to_markdown(metrics.round(4)),
            "",
            "## Primary paired comparisons",
            "",
            dataframe_to_markdown(comparisons.round(4)),
            "",
            "## Repeated degree-preserving nulls",
            "",
            (
                f"Across {manifest['null_repeats']} graph realizations, random RMSE ranged from "
                f"{manifest['random_rmse_range'][0]:.4f} to "
                f"{manifest['random_rmse_range'][1]:.4f} SD. Structural RMSE improvement ranged "
                f"from {manifest['rmse_improvement_range'][0]:.4f} to "
                f"{manifest['rmse_improvement_range'][1]:.4f} SD."
            ),
            "",
            dataframe_to_markdown(
                null_results.loc[
                    :,
                    [
                        "null_seed",
                        "random_rmse_sd",
                        "rmse_improvement_sd",
                        "rmse_ci_lower",
                        "random_precision_at_k",
                        "precision_improvement",
                        "precision_ci_lower",
                    ],
                ].round(4)
            ),
            "",
            "## Boundary",
            "",
            manifest["claim_boundary"],
            "",
        ]
    )


def run(config_path: Path, output_dir: Path) -> dict[str, Any]:
    config_path = config_path.resolve()
    config = yaml.safe_load(config_path.read_text())
    protocol_path = _resolve(config_path, str(config["protocol"]))
    cache_dir = _resolve(config_path, str(config["source"]["cache_dir"]))
    study_id = str(config["source"]["study_id"])
    analysis_id = str(config["source"]["analysis_id"])
    client = MetabolomicsWorkbenchClient(cache_dir=cache_dir)

    factors = client.factors(study_id)
    groups = _annotate_groups(
        factors, str(config["columns"]["validation_group_factor"])
    ).dropna()
    measurements = client.measurements(study_id)
    selected = measurements.loc[measurements["analysis_id"].eq(analysis_id)].copy()
    if selected.empty:
        raise ValueError(f"No measurements were returned for {study_id}/{analysis_id}.")
    matrix = selected.pivot_table(
        index="local_sample_id",
        columns="metabolite_name",
        values="value",
        aggfunc="mean",
    )
    shared = matrix.index.astype(str).intersection(groups.index)
    matrix = matrix.loc[shared]
    groups = groups.loc[shared]
    complete = matrix.columns[matrix.notna().all()]
    families = eligible_lipid_families(
        complete,
        minimum_family_features=int(
            config["eligibility"]["minimum_family_features"]
        ),
    )
    values = matrix.loc[:, families.index]
    minimum_features = int(config["eligibility"]["minimum_total_features"])
    minimum_groups = int(config["eligibility"]["minimum_validation_groups"])
    if len(families) < minimum_features:
        raise ValueError(
            f"Only {len(families)} complete eligible features; {minimum_features} required."
        )
    if groups.nunique() < minimum_groups:
        raise ValueError(
            f"Only {groups.nunique()} validation groups; {minimum_groups} required."
        )

    masks = disjoint_family_masks(
        families,
        masks=int(config["masking"]["masks"]),
        seed=int(config["masking"]["seed"]),
    )
    incidence = structural_descriptor_incidence(
        families.index,
        minimum_features=int(config["descriptors"]["minimum_features"]),
        maximum_feature_fraction=float(
            config["descriptors"]["maximum_feature_fraction"]
        ),
    )
    coarse_metrics, coarse_samples, _ = masked_panel_benchmark(
        values,
        groups,
        families,
        masks,
        folds=int(config["model"]["folds"]),
        ridge_alpha=float(config["model"]["ridge_alpha"]),
        random_group_seed=int(config["model"]["random_group_seed"]),
        priority_count=int(config["priority"]["priority_count"]),
        true_priority_fraction=float(config["priority"]["true_priority_fraction"]),
    )
    primary_seed = int(config["null_ensemble"]["first_seed"])
    structural_metrics, structural_samples, primary_audit = (
        masked_structural_descriptor_benchmark(
            values,
            groups,
            masks,
            incidence,
            folds=int(config["model"]["folds"]),
            ridge_alpha=float(config["model"]["ridge_alpha"]),
            split_seed=int(config["model"]["split_seed"]),
            null_seed=primary_seed,
            swaps_per_edge=int(config["descriptors"]["swaps_per_edge"]),
            priority_count=int(config["priority"]["priority_count"]),
            true_priority_fraction=float(config["priority"]["true_priority_fraction"]),
        )
    )
    settings: dict[str, int | float] = {
        "folds": int(config["model"]["folds"]),
        "ridge_alpha": float(config["model"]["ridge_alpha"]),
        "split_seed": int(config["model"]["split_seed"]),
        "swaps_per_edge": int(config["descriptors"]["swaps_per_edge"]),
        "priority_count": int(config["priority"]["priority_count"]),
        "true_priority_fraction": float(
            config["priority"]["true_priority_fraction"]
        ),
        "bootstrap_draws": int(config["inference"]["bootstrap_draws"]),
        "bootstrap_seed": int(config["inference"]["seed"]),
    }
    repeats = int(config["null_ensemble"]["repeats"])
    jobs = min(int(config["null_ensemble"]["parallel_jobs"]), repeats)
    with parallel_backend("loky", inner_max_num_threads=1):
        null_records = Parallel(n_jobs=jobs)(
            delayed(_run_null_seed)(
                index,
                primary_seed + index,
                values,
                groups,
                masks,
                incidence,
                settings,
            )
            for index in range(repeats)
        )
    null_results = pd.DataFrame.from_records(null_records).sort_values("null_index")

    comparison_samples = pd.concat([coarse_samples, structural_samples], ignore_index=True)
    comparison_records = []
    for index, reference in enumerate(
        ["population_mean", "family_score_ridge", "all_visible_ridge"]
    ):
        comparison_records.append(
            bootstrap_subject_model_difference(
                comparison_samples,
                reference=reference,
                challenger="structural_descriptor_ridge",
                draws=int(config["inference"]["bootstrap_draws"]),
                seed=int(config["inference"]["seed"]) + repeats + index,
            )
        )
    comparisons = pd.DataFrame.from_records(comparison_records)
    threshold = float(config["gates"]["bootstrap_ci_lower_must_exceed"])
    decision, gates = adjudicate_external_replication(
        null_results,
        eligible_features=len(families),
        validation_groups=int(groups.nunique()),
        minimum_features=minimum_features,
        minimum_groups=minimum_groups,
        threshold=threshold,
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    feature_path = output_dir / "eligible-features.csv"
    metrics_path = output_dir / "model-metrics.csv"
    comparisons_path = output_dir / "bootstrap-comparisons.csv"
    null_path = output_dir / "null-ensemble-results.csv"
    audit_path = output_dir / "primary-null-audit.csv"
    report_path = output_dir / "report.md"
    manifest_path = output_dir / "manifest.json"
    pd.DataFrame(
        {"feature": families.index, "family": families.to_numpy()}
    ).to_csv(feature_path, index=False)
    structural_nonpopulation = structural_metrics.loc[
        ~structural_metrics["model"].eq("population_mean")
    ]
    combined_metrics = pd.concat(
        [coarse_metrics, structural_nonpopulation], ignore_index=True
    )
    combined_metrics.to_csv(metrics_path, index=False)
    comparisons.to_csv(comparisons_path, index=False)
    null_results.to_csv(null_path, index=False)
    primary_audit.to_csv(audit_path, index=False)

    factors_url = f"{DEFAULT_BASE_URL}/study/study_id/{study_id}/factors/json"
    data_url = f"{DEFAULT_BASE_URL}/study/study_id/{study_id}/data/json"
    source_paths = {
        "factors": _cache_path(cache_dir, factors_url),
        "measurements": _cache_path(cache_dir, data_url),
    }
    if any(not path.exists() for path in source_paths.values()):
        raise ValueError("Expected checksum-pinned Workbench cache files are absent.")
    manifest: dict[str, Any] = {
        "analysis_id": str(config["analysis_id"]),
        "completed_at": datetime.now(UTC).isoformat(),
        "decision": decision,
        "independent_public_cohort": True,
        "config": str(config_path),
        "config_sha256": _sha256(config_path),
        "protocol": str(protocol_path),
        "protocol_sha256": _sha256(protocol_path),
        "implementation_sha256": implementation_sha256(),
        "source_sha256": {key: _sha256(path) for key, path in source_paths.items()},
        "study_id": study_id,
        "study_analysis_id": analysis_id,
        "samples": len(values),
        "validation_groups": int(groups.nunique()),
        "eligible_features": len(families),
        "eligible_families": int(families.nunique()),
        "descriptors": incidence.shape[1],
        "null_repeats": repeats,
        "parallel_jobs": jobs,
        "model_metrics": combined_metrics.replace({np.nan: None}).to_dict(
            orient="records"
        ),
        "bootstrap_comparisons": comparison_records,
        "random_rmse_range": [
            float(null_results["random_rmse_sd"].min()),
            float(null_results["random_rmse_sd"].max()),
        ],
        "rmse_improvement_range": [
            float(null_results["rmse_improvement_sd"].min()),
            float(null_results["rmse_improvement_sd"].max()),
        ],
        "precision_improvement_range": [
            float(null_results["precision_improvement"].min()),
            float(null_results["precision_improvement"].max()),
        ],
        "minimum_rmse_ci_lower": float(null_results["rmse_ci_lower"].min()),
        "minimum_precision_ci_lower": float(
            null_results["precision_ci_lower"].min()
        ),
        "gate_results": gates,
        "claim_boundary": str(config["claim_boundary"]),
    }
    report_path.write_text(
        _render_report(manifest, combined_metrics, comparisons, null_results)
    )
    outputs = [
        feature_path,
        metrics_path,
        comparisons_path,
        null_path,
        audit_path,
        report_path,
    ]
    manifest["output_sha256"] = {path.name: _sha256(path) for path in outputs}
    manifest_path.write_text(json.dumps(manifest, indent=2, allow_nan=False) + "\n")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        type=Path,
        default=Path("config/pathway-score-st000818-replication.yaml"),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("artifacts/pathway-score-st000818-replication"),
    )
    args = parser.parse_args()
    print(json.dumps(run(args.config, args.output), indent=2))


if __name__ == "__main__":
    main()
