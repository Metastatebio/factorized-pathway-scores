"""Run repeated graph-null sensitivity for adaptive ST002081 structural scores."""

from __future__ import annotations

import argparse
import hashlib
import inspect
import json
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import pandas as pd
import yaml
from joblib import Parallel, delayed, parallel_backend

from ..analysis.pathway_scores import (
    bootstrap_subject_model_difference,
    disjoint_family_masks,
    eligible_lipid_families,
    masked_structural_descriptor_benchmark,
    structural_descriptor_incidence,
)
from ..analysis.publication_readiness import load_verified_manifest
from ..datasets.mwtab import load_mwtab
from ..reporting import dataframe_to_markdown


def _resolve(config_path: Path, value: str) -> Path:
    path = Path(value)
    return path if path.is_absolute() else (config_path.parent / path).resolve()


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def implementation_sha256() -> str:
    digest = hashlib.sha256()
    for implementation in (run, masked_structural_descriptor_benchmark):
        path = Path(inspect.getfile(implementation))
        digest.update(path.name.encode())
        digest.update(_sha256(path).encode())
    return digest.hexdigest()


def _run_null_seed(
    null_index: int,
    null_seed: int,
    values: pd.DataFrame,
    subjects: pd.Series,
    masks: pd.DataFrame,
    incidence: pd.DataFrame,
    settings: dict[str, int | float],
) -> dict[str, Any]:
    metrics, samples, audit = masked_structural_descriptor_benchmark(
        values,
        subjects,
        masks,
        incidence,
        folds=int(settings["folds"]),
        ridge_alpha=float(settings["ridge_alpha"]),
        split_seed=int(settings["split_seed"]),
        null_seed=null_seed,
        swaps_per_edge=int(settings["swaps_per_edge"]),
        priority_count=int(settings["priority_count"]),
        true_priority_fraction=float(settings["true_priority_fraction"]),
    )
    comparison = bootstrap_subject_model_difference(
        samples,
        reference="degree_preserving_random_structural_ridge",
        challenger="structural_descriptor_ridge",
        draws=int(settings["bootstrap_draws"]),
        seed=int(settings["bootstrap_seed"]) + null_index,
    )
    by_model = metrics.set_index("model")
    random_row = by_model.loc["degree_preserving_random_structural_ridge"]
    structural_row = by_model.loc["structural_descriptor_ridge"]
    return {
        "null_index": null_index,
        "null_seed": null_seed,
        "random_rmse_sd": float(random_row["row_weighted_rmse_sd"]),
        "structural_rmse_sd": float(structural_row["row_weighted_rmse_sd"]),
        "random_precision_at_k": float(random_row["precision_at_k"]),
        "structural_precision_at_k": float(structural_row["precision_at_k"]),
        **comparison,
        "row_degrees_preserved": bool(audit["row_degrees_preserved"].all()),
        "column_degrees_preserved": bool(audit["column_degrees_preserved"].all()),
    }


def adjudicate_null_ensemble(
    results: pd.DataFrame, *, threshold: float
) -> tuple[str, dict[str, bool]]:
    """Apply the post-result all-null-realizations robustness rule."""
    gates = {
        "all_degrees_preserved": bool(
            results["row_degrees_preserved"].all()
            and results["column_degrees_preserved"].all()
        ),
        "rmse_ci_passes_all_nulls": bool(results["rmse_ci_lower"].gt(threshold).all()),
        "precision_ci_passes_all_nulls": bool(
            results["precision_ci_lower"].gt(threshold).all()
        ),
    }
    decision = (
        "ADAPTIVE_STRUCTURE_SIGNAL_ROBUST_ACROSS_GRAPH_NULLS"
        if all(gates.values())
        else "ADAPTIVE_STRUCTURE_GRAPH_NULL_SENSITIVITY_FAILED"
    )
    return decision, gates


def _render_report(manifest: dict[str, Any], results: pd.DataFrame) -> str:
    columns = [
        "null_seed",
        "random_rmse_sd",
        "rmse_improvement_sd",
        "rmse_ci_lower",
        "rmse_ci_upper",
        "random_precision_at_k",
        "precision_improvement",
        "precision_ci_lower",
        "precision_ci_upper",
    ]
    return "\n".join(
        [
            "# ST002081 repeated graph-null sensitivity",
            "",
            f"**Decision:** `{manifest['decision']}`  ",
            "**Role:** post-result adaptive sensitivity; independent replication: false",
            "",
            "## Summary",
            "",
            (
                f"All {manifest['null_repeats']} degree-preserving graph realizations were "
                f"evaluated. Random-graph RMSE ranged from "
                f"{manifest['random_rmse_range'][0]:.4f} to "
                f"{manifest['random_rmse_range'][1]:.4f} SD. The worst participant-bootstrap "
                f"lower bound for structural RMSE improvement was "
                f"{manifest['minimum_rmse_ci_lower']:.4f} SD; the worst lower bound for "
                f"precision improvement was {manifest['minimum_precision_ci_lower']:.4f}."
            ),
            "",
            "## Complete null ensemble",
            "",
            dataframe_to_markdown(results.loc[:, columns].round(4)),
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
    base_config_path = _resolve(config_path, str(config["base_config"]))
    base_config = yaml.safe_load(base_config_path.read_text())
    base_artifact = _resolve(config_path, str(config["base_artifact"]))
    _, base_integrity = load_verified_manifest(base_artifact)
    if base_integrity != "verified":
        raise ValueError(f"Base structural artifact integrity failed: {base_integrity}")

    source_path = _resolve(base_config_path, str(base_config["source"]["mwtab"]))
    dataset = load_mwtab(source_path)
    subject_column = str(base_config["columns"]["subject"])
    excluded = {
        str(value) for value in base_config["eligibility"]["exclude_subject_values"]
    }
    metadata = dataset.sample_metadata.copy()
    metadata = metadata.loc[~metadata[subject_column].astype(str).isin(excluded)]
    metabolomics = dataset.blocks["metabolomics"].loc[metadata.index]
    complete = metabolomics.columns[metabolomics.notna().all()]
    families = eligible_lipid_families(
        complete,
        minimum_family_features=int(
            base_config["eligibility"]["minimum_family_features"]
        ),
    )
    values = metabolomics.loc[:, families.index]
    masks = disjoint_family_masks(
        families,
        masks=int(base_config["masking"]["masks"]),
        seed=int(base_config["masking"]["seed"]),
    )
    incidence = structural_descriptor_incidence(
        families.index,
        minimum_features=int(base_config["descriptors"]["minimum_features"]),
        maximum_feature_fraction=float(
            base_config["descriptors"]["maximum_feature_fraction"]
        ),
    )
    settings: dict[str, int | float] = {
        "folds": int(base_config["model"]["folds"]),
        "ridge_alpha": float(base_config["model"]["ridge_alpha"]),
        "split_seed": int(base_config["model"]["split_seed"]),
        "swaps_per_edge": int(base_config["descriptors"]["swaps_per_edge"]),
        "priority_count": int(base_config["priority"]["priority_count"]),
        "true_priority_fraction": float(
            base_config["priority"]["true_priority_fraction"]
        ),
        "bootstrap_draws": int(config["inference"]["bootstrap_draws"]),
        "bootstrap_seed": int(config["inference"]["seed"]),
    }
    repeats = int(config["null_ensemble"]["repeats"])
    first_seed = int(config["null_ensemble"]["first_seed"])
    jobs = min(int(config["null_ensemble"]["parallel_jobs"]), repeats)
    with parallel_backend("loky", inner_max_num_threads=1):
        records = Parallel(n_jobs=jobs)(
            delayed(_run_null_seed)(
                index,
                first_seed + index,
                values,
                metadata[subject_column],
                masks,
                incidence,
                settings,
            )
            for index in range(repeats)
        )
    results = pd.DataFrame.from_records(records).sort_values("null_index")
    threshold = float(config["gates"]["bootstrap_ci_lower_must_exceed"])
    decision, gates = adjudicate_null_ensemble(results, threshold=threshold)

    output_dir.mkdir(parents=True, exist_ok=True)
    results_path = output_dir / "null-ensemble-results.csv"
    report_path = output_dir / "report.md"
    manifest_path = output_dir / "manifest.json"
    results.to_csv(results_path, index=False)
    manifest: dict[str, Any] = {
        "analysis_id": str(config["analysis_id"]),
        "completed_at": datetime.now(UTC).isoformat(),
        "decision": decision,
        "adaptive_post_result_sensitivity": True,
        "config": str(config_path),
        "config_sha256": _sha256(config_path),
        "base_config_sha256": _sha256(base_config_path),
        "protocol": str(protocol_path),
        "protocol_sha256": _sha256(protocol_path),
        "implementation_sha256": implementation_sha256(),
        "source_sha256": _sha256(source_path),
        "base_artifact_integrity": base_integrity,
        "null_repeats": repeats,
        "parallel_jobs": jobs,
        "random_rmse_range": [
            float(results["random_rmse_sd"].min()),
            float(results["random_rmse_sd"].max()),
        ],
        "rmse_improvement_range": [
            float(results["rmse_improvement_sd"].min()),
            float(results["rmse_improvement_sd"].max()),
        ],
        "precision_improvement_range": [
            float(results["precision_improvement"].min()),
            float(results["precision_improvement"].max()),
        ],
        "minimum_rmse_ci_lower": float(results["rmse_ci_lower"].min()),
        "minimum_precision_ci_lower": float(results["precision_ci_lower"].min()),
        "gate_results": gates,
        "claim_boundary": str(config["claim_boundary"]),
    }
    report_path.write_text(_render_report(manifest, results))
    manifest["output_sha256"] = {
        path.name: _sha256(path) for path in (results_path, report_path)
    }
    manifest_path.write_text(json.dumps(manifest, indent=2, allow_nan=False) + "\n")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        type=Path,
        default=Path("config/pathway-score-st002081-structural-null-sensitivity.yaml"),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("artifacts/pathway-score-st002081-structural-null-sensitivity"),
    )
    args = parser.parse_args()
    print(json.dumps(run(args.config, args.output), indent=2))


if __name__ == "__main__":
    main()
