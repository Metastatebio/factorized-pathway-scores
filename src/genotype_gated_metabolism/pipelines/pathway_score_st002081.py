"""Benchmark factorized lipid-family scores on hidden ST002081 panels."""

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
    masked_panel_benchmark,
)
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
    for implementation in (run, masked_panel_benchmark, disjoint_family_masks):
        path = Path(inspect.getfile(implementation))
        digest.update(path.name.encode())
        digest.update(_sha256(path).encode())
    return digest.hexdigest()


def _render_report(manifest: dict[str, Any], metrics: pd.DataFrame) -> str:
    display = metrics[
        [
            "model",
            "row_weighted_rmse_sd",
            "equal_subject_weighted_rmse_sd",
            "pooled_r2",
            "precision_at_k",
            "ndcg_at_k",
        ]
    ].copy()
    return "\n".join(
        [
            "# ST002081 factorized pathway missing-panel benchmark",
            "",
            f"**Decision:** `{manifest['decision']}`  ",
            "**Directional forecasting tested:** false  ",
            "**Participant genomics available:** false",
            "",
            "## Result",
            "",
            (
                f"The locked benchmark evaluated {manifest['eligible_features']} complete lipid "
                f"features from {manifest['eligible_families']} biochemical families in "
                f"{manifest['samples']} samples from {manifest['subjects']} people. Every feature "
                f"was hidden exactly once across {manifest['masks']} disjoint masks."
            ),
            "",
            dataframe_to_markdown(display.round(4)),
            "",
            "## Paired participant bootstrap",
            "",
            *[
                (
                    f"- `{item['challenger']}` versus `{item['reference']}`: RMSE improvement "
                    f"{item['rmse_improvement_sd']:.4f} SD (95% interval "
                    f"{item['rmse_ci_lower']:.4f} to {item['rmse_ci_upper']:.4f}); precision@k "
                    f"improvement {item['precision_improvement']:.4f} "
                    f"({item['precision_ci_lower']:.4f} to "
                    f"{item['precision_ci_upper']:.4f})."
                )
                for item in manifest["bootstrap_comparisons"]
            ],
            "",
            "## Interpretation boundary",
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
    dataset = load_mwtab(source_path)
    subject_column = str(config["columns"]["subject"])
    excluded = {str(value) for value in config["eligibility"]["exclude_subject_values"]}
    metadata = dataset.sample_metadata.copy()
    keep = ~metadata[subject_column].astype(str).isin(excluded)
    metadata = metadata.loc[keep]
    metabolomics = dataset.blocks["metabolomics"].loc[metadata.index]
    complete = metabolomics.columns[metabolomics.notna().all()]
    family_by_feature = eligible_lipid_families(
        complete,
        minimum_family_features=int(config["eligibility"]["minimum_family_features"]),
    )
    values = metabolomics.loc[:, family_by_feature.index]
    masks = disjoint_family_masks(
        family_by_feature,
        masks=int(config["masking"]["masks"]),
        seed=int(config["masking"]["seed"]),
    )
    model_metrics, sample_metrics, target_metrics = masked_panel_benchmark(
        values,
        metadata[subject_column],
        family_by_feature,
        masks,
        folds=int(config["model"]["folds"]),
        ridge_alpha=float(config["model"]["ridge_alpha"]),
        random_group_seed=int(config["model"]["random_group_seed"]),
        priority_count=int(config["priority"]["priority_count"]),
        true_priority_fraction=float(config["priority"]["true_priority_fraction"]),
    )
    comparisons = [
        bootstrap_subject_model_difference(
            sample_metrics,
            reference=reference,
            challenger="family_score_ridge",
            draws=int(config["inference"]["bootstrap_draws"]),
            seed=int(config["inference"]["seed"]) + index,
        )
        for index, reference in enumerate(
            ["population_mean", "random_group_score_ridge", "all_visible_ridge"]
        )
    ]
    comparison_by_reference = {item["reference"]: item for item in comparisons}
    threshold = float(config["gates"]["bootstrap_ci_lower_must_exceed"])
    gate_results = {
        "data_scale": len(family_by_feature) >= int(
            config["eligibility"]["minimum_total_features"]
        )
        and family_by_feature.nunique() >= int(config["eligibility"]["minimum_families"]),
        "rmse_beats_population": comparison_by_reference["population_mean"][
            "rmse_ci_lower"
        ]
        > threshold,
        "rmse_beats_random_grouping": comparison_by_reference[
            "random_group_score_ridge"
        ]["rmse_ci_lower"]
        > threshold,
        "precision_beats_random_grouping": comparison_by_reference[
            "random_group_score_ridge"
        ]["precision_ci_lower"]
        > threshold,
    }
    decision = (
        "PATHWAY_MISSING_PANEL_SIGNAL"
        if all(gate_results.values())
        else "RESOURCE_ONLY"
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    family_path = output_dir / "feature-families.csv"
    masks_path = output_dir / "disjoint-masks.csv"
    model_path = output_dir / "model-metrics.csv"
    sample_path = output_dir / "sample-metrics.csv"
    target_path = output_dir / "target-metrics.csv"
    comparisons_path = output_dir / "bootstrap-comparisons.csv"
    report_path = output_dir / "report.md"
    manifest_path = output_dir / "manifest.json"
    family_by_feature.rename_axis("feature").reset_index().to_csv(family_path, index=False)
    masks.to_csv(masks_path, index=False)
    model_metrics.to_csv(model_path, index=False)
    sample_metrics.to_csv(sample_path, index=False)
    target_metrics.to_csv(target_path, index=False)
    pd.DataFrame.from_records(comparisons).to_csv(comparisons_path, index=False)
    manifest: dict[str, Any] = {
        "analysis_id": str(config["analysis_id"]),
        "completed_at": datetime.now(UTC).isoformat(),
        "decision": decision,
        "config": str(config_path),
        "config_sha256": _sha256(config_path),
        "protocol": str(protocol_path),
        "protocol_sha256": _sha256(protocol_path),
        "implementation_sha256": implementation_sha256(),
        "source_sha256": _sha256(source_path),
        "samples": len(values),
        "subjects": int(metadata[subject_column].astype(str).nunique()),
        "source_features": metabolomics.shape[1],
        "complete_features": len(complete),
        "eligible_features": len(family_by_feature),
        "eligible_families": int(family_by_feature.nunique()),
        "masks": int(config["masking"]["masks"]),
        "model_metrics": model_metrics.replace({np.nan: None}).to_dict(orient="records"),
        "bootstrap_comparisons": comparisons,
        "gate_results": gate_results,
        "claim_boundary": str(config["claim_boundary"]),
    }
    report_path.write_text(_render_report(manifest, model_metrics))
    output_paths = [
        family_path,
        masks_path,
        model_path,
        sample_path,
        target_path,
        comparisons_path,
        report_path,
    ]
    manifest["output_sha256"] = {path.name: _sha256(path) for path in output_paths}
    manifest_path.write_text(json.dumps(manifest, indent=2, allow_nan=False) + "\n")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=Path("config/pathway-score-st002081.yaml"))
    parser.add_argument(
        "--output", type=Path, default=Path("artifacts/pathway-score-st002081")
    )
    args = parser.parse_args()
    print(json.dumps(run(args.config, args.output), indent=2))


if __name__ == "__main__":
    main()
