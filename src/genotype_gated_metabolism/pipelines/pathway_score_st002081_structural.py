"""Run the adaptive structure-resolved ST002081 pathway-score benchmark."""

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


def _render_report(manifest: dict[str, Any], metrics: pd.DataFrame) -> str:
    return "\n".join(
        [
            "# Adaptive ST002081 structure-resolved pathway-score benchmark",
            "",
            f"**Decision:** `{manifest['decision']}`  ",
            "**Independent replication:** false",
            "",
            "## Result",
            "",
            (
                f"The adaptive score used {manifest['descriptors']} filtered structural "
                f"descriptors across {manifest['eligible_features']} lipids. Every null preserved "
                "the visible feature and descriptor degree sequences exactly."
            ),
            "",
            dataframe_to_markdown(metrics.round(4)),
            "",
            "## Paired participant bootstrap",
            "",
            *[
                (
                    f"- `{item['challenger']}` versus `{item['reference']}`: RMSE improvement "
                    f"{item['rmse_improvement_sd']:.4f} SD "
                    f"({item['rmse_ci_lower']:.4f} to {item['rmse_ci_upper']:.4f}); "
                    f"precision improvement {item['precision_improvement']:.4f} "
                    f"({item['precision_ci_lower']:.4f} to {item['precision_ci_upper']:.4f})."
                )
                for item in manifest["bootstrap_comparisons"]
            ],
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
    source_path = _resolve(config_path, str(config["source"]["mwtab"]))
    coarse_dir = _resolve(config_path, str(config["source"]["coarse_artifact"]))
    _, coarse_integrity = load_verified_manifest(coarse_dir)
    if coarse_integrity != "verified":
        raise ValueError(f"Coarse benchmark integrity failed: {coarse_integrity}")
    dataset = load_mwtab(source_path)
    subject_column = str(config["columns"]["subject"])
    excluded = {str(value) for value in config["eligibility"]["exclude_subject_values"]}
    metadata = dataset.sample_metadata.copy()
    keep = ~metadata[subject_column].astype(str).isin(excluded)
    metadata = metadata.loc[keep]
    metabolomics = dataset.blocks["metabolomics"].loc[metadata.index]
    complete = metabolomics.columns[metabolomics.notna().all()]
    families = eligible_lipid_families(
        complete,
        minimum_family_features=int(config["eligibility"]["minimum_family_features"]),
    )
    values = metabolomics.loc[:, families.index]
    masks = disjoint_family_masks(
        families,
        masks=int(config["masking"]["masks"]),
        seed=int(config["masking"]["seed"]),
    )
    incidence = structural_descriptor_incidence(
        families.index,
        minimum_features=int(config["descriptors"]["minimum_features"]),
        maximum_feature_fraction=float(config["descriptors"]["maximum_feature_fraction"]),
    )
    metrics, sample_metrics, null_audit = masked_structural_descriptor_benchmark(
        values,
        metadata[subject_column],
        masks,
        incidence,
        folds=int(config["model"]["folds"]),
        ridge_alpha=float(config["model"]["ridge_alpha"]),
        split_seed=int(config["model"]["split_seed"]),
        null_seed=int(config["model"]["null_seed"]),
        swaps_per_edge=int(config["descriptors"]["swaps_per_edge"]),
        priority_count=int(config["priority"]["priority_count"]),
        true_priority_fraction=float(config["priority"]["true_priority_fraction"]),
    )
    coarse_samples = pd.read_csv(coarse_dir / "sample-metrics.csv")
    comparison_samples = pd.concat([sample_metrics, coarse_samples], ignore_index=True)
    references = [
        "population_mean",
        "degree_preserving_random_structural_ridge",
        "family_score_ridge",
        "all_visible_ridge",
    ]
    comparisons = [
        bootstrap_subject_model_difference(
            comparison_samples,
            reference=reference,
            challenger="structural_descriptor_ridge",
            draws=int(config["inference"]["bootstrap_draws"]),
            seed=int(config["inference"]["seed"]) + index,
        )
        for index, reference in enumerate(references)
    ]
    by_reference = {item["reference"]: item for item in comparisons}
    threshold = float(config["gates"]["bootstrap_ci_lower_must_exceed"])
    null_comparison = by_reference["degree_preserving_random_structural_ridge"]
    gate_results = {
        "null_degrees_preserved": bool(
            null_audit["row_degrees_preserved"].all()
            and null_audit["column_degrees_preserved"].all()
        ),
        "rmse_beats_degree_preserving_null": null_comparison["rmse_ci_lower"] > threshold,
        "precision_beats_degree_preserving_null": null_comparison["precision_ci_lower"]
        > threshold,
    }
    if all(gate_results.values()):
        decision = "ADAPTIVE_STRUCTURE_RESOLUTION_SIGNAL"
    elif gate_results["rmse_beats_degree_preserving_null"]:
        decision = "ADAPTIVE_STRUCTURAL_RECONSTRUCTION_ONLY"
    else:
        decision = "ADAPTIVE_STRUCTURE_GATE_FAILED"

    output_dir.mkdir(parents=True, exist_ok=True)
    incidence_path = output_dir / "structural-descriptor-incidence.csv"
    metrics_path = output_dir / "model-metrics.csv"
    sample_path = output_dir / "sample-metrics.csv"
    null_path = output_dir / "null-audit.csv"
    comparisons_path = output_dir / "bootstrap-comparisons.csv"
    report_path = output_dir / "report.md"
    manifest_path = output_dir / "manifest.json"
    incidence.astype(int).rename_axis("feature").reset_index().to_csv(incidence_path, index=False)
    metrics.to_csv(metrics_path, index=False)
    sample_metrics.to_csv(sample_path, index=False)
    null_audit.to_csv(null_path, index=False)
    pd.DataFrame.from_records(comparisons).to_csv(comparisons_path, index=False)
    manifest: dict[str, Any] = {
        "analysis_id": str(config["analysis_id"]),
        "completed_at": datetime.now(UTC).isoformat(),
        "decision": decision,
        "adaptive_followup": True,
        "config": str(config_path),
        "config_sha256": _sha256(config_path),
        "protocol": str(protocol_path),
        "protocol_sha256": _sha256(protocol_path),
        "implementation_sha256": implementation_sha256(),
        "source_sha256": _sha256(source_path),
        "coarse_artifact_integrity": coarse_integrity,
        "samples": len(values),
        "subjects": int(metadata[subject_column].astype(str).nunique()),
        "eligible_features": len(families),
        "descriptors": incidence.shape[1],
        "model_metrics": metrics.replace({np.nan: None}).to_dict(orient="records"),
        "bootstrap_comparisons": comparisons,
        "gate_results": gate_results,
        "claim_boundary": str(config["claim_boundary"]),
    }
    report_path.write_text(_render_report(manifest, metrics))
    outputs = [incidence_path, metrics_path, sample_path, null_path, comparisons_path, report_path]
    manifest["output_sha256"] = {path.name: _sha256(path) for path in outputs}
    manifest_path.write_text(json.dumps(manifest, indent=2, allow_nan=False) + "\n")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config", type=Path, default=Path("config/pathway-score-st002081-structural.yaml")
    )
    parser.add_argument(
        "--output", type=Path, default=Path("artifacts/pathway-score-st002081-structural")
    )
    args = parser.parse_args()
    print(json.dumps(run(args.config, args.output), indent=2))


if __name__ == "__main__":
    main()
