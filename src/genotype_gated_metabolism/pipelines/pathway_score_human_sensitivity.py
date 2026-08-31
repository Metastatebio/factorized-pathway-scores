"""Run post-result human structural pathway-score sensitivity grids."""

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

from ..analysis.intervention_registry import parse_factor_string
from ..analysis.pathway_scores import (
    bootstrap_subject_model_difference,
    disjoint_family_masks,
    eligible_lipid_families,
    masked_structural_descriptor_benchmark,
    structural_descriptor_incidence,
)
from ..analysis.publication_readiness import load_verified_manifest
from ..datasets.metabolomics_workbench import MetabolomicsWorkbenchClient
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


def _st002081_inputs(base_path: Path) -> tuple[pd.DataFrame, pd.Series, pd.DataFrame, dict[str, Any]]:
    config = yaml.safe_load(base_path.read_text())
    source_path = _resolve(base_path, str(config["source"]["mwtab"]))
    dataset = load_mwtab(source_path)
    subject_column = str(config["columns"]["subject"])
    excluded = {str(value) for value in config["eligibility"]["exclude_subject_values"]}
    metadata = dataset.sample_metadata.copy()
    metadata = metadata.loc[~metadata[subject_column].astype(str).isin(excluded)]
    metabolomics = dataset.blocks["metabolomics"].loc[metadata.index]
    complete = metabolomics.columns[metabolomics.notna().all()]
    families = eligible_lipid_families(
        complete,
        minimum_family_features=int(config["eligibility"]["minimum_family_features"]),
    )
    masks = disjoint_family_masks(
        families,
        masks=int(config["masking"]["masks"]),
        seed=int(config["masking"]["seed"]),
    )
    return metabolomics.loc[:, families.index], metadata[subject_column], masks, config


def _st000818_inputs(base_path: Path) -> tuple[pd.DataFrame, pd.Series, pd.DataFrame, dict[str, Any]]:
    config = yaml.safe_load(base_path.read_text())
    cache_dir = _resolve(base_path, str(config["source"]["cache_dir"]))
    study_id = str(config["source"]["study_id"])
    analysis_id = str(config["source"]["analysis_id"])
    client = MetabolomicsWorkbenchClient(cache_dir=cache_dir)
    factors = client.factors(study_id)
    factor_name = str(config["columns"]["validation_group_factor"])
    parsed = factors["factors"].map(parse_factor_string)
    groups = pd.Series(
        parsed.map(lambda values: values.get(factor_name)).to_numpy(),
        index=factors["local_sample_id"].astype(str),
        dtype="object",
    ).dropna()
    selected = client.measurements(study_id)
    selected = selected.loc[selected["analysis_id"].eq(analysis_id)]
    matrix = selected.pivot_table(
        index="local_sample_id", columns="metabolite_name", values="value", aggfunc="mean"
    )
    shared = matrix.index.astype(str).intersection(groups.index)
    matrix = matrix.loc[shared]
    groups = groups.loc[shared]
    complete = matrix.columns[matrix.notna().all()]
    families = eligible_lipid_families(
        complete,
        minimum_family_features=int(config["eligibility"]["minimum_family_features"]),
    )
    masks = disjoint_family_masks(
        families,
        masks=int(config["masking"]["masks"]),
        seed=int(config["masking"]["seed"]),
    )
    return matrix.loc[:, families.index], groups, masks, config


def _run_setting(
    dataset_id: str,
    setting_index: int,
    setting: dict[str, Any],
    values: pd.DataFrame,
    groups: pd.Series,
    masks: pd.DataFrame,
    base: dict[str, Any],
    bootstrap_draws: int,
    bootstrap_seed: int,
) -> dict[str, Any]:
    incidence = structural_descriptor_incidence(
        values.columns,
        minimum_features=int(setting["minimum_features"]),
        maximum_feature_fraction=float(setting["maximum_feature_fraction"]),
    )
    metrics, samples, audit = masked_structural_descriptor_benchmark(
        values,
        groups,
        masks,
        incidence,
        folds=int(base["model"]["folds"]),
        ridge_alpha=float(setting["ridge_alpha"]),
        split_seed=int(base["model"]["split_seed"]),
        null_seed=bootstrap_seed + 1000 + setting_index * 10,
        swaps_per_edge=int(setting["swaps_per_edge"]),
        priority_count=int(base["priority"]["priority_count"]),
        true_priority_fraction=float(base["priority"]["true_priority_fraction"]),
    )
    comparison = bootstrap_subject_model_difference(
        samples,
        reference="degree_preserving_random_structural_ridge",
        challenger="structural_descriptor_ridge",
        draws=bootstrap_draws,
        seed=bootstrap_seed + setting_index,
    )
    by_model = metrics.set_index("model")
    structural = by_model.loc["structural_descriptor_ridge"]
    random = by_model.loc["degree_preserving_random_structural_ridge"]
    return {
        "dataset": dataset_id,
        "setting": str(setting["id"]),
        "ridge_alpha": float(setting["ridge_alpha"]),
        "minimum_features": int(setting["minimum_features"]),
        "maximum_feature_fraction": float(setting["maximum_feature_fraction"]),
        "swaps_per_edge": int(setting["swaps_per_edge"]),
        "descriptors": incidence.shape[1],
        "structural_rmse_sd": float(structural["row_weighted_rmse_sd"]),
        "random_rmse_sd": float(random["row_weighted_rmse_sd"]),
        "structural_precision_at_k": float(structural["precision_at_k"]),
        "random_precision_at_k": float(random["precision_at_k"]),
        **comparison,
        "row_degrees_preserved": bool(audit["row_degrees_preserved"].all()),
        "column_degrees_preserved": bool(audit["column_degrees_preserved"].all()),
    }


def adjudicate_human_sensitivity(
    results: pd.DataFrame, *, point_threshold: float, minimum_ci_pass_rate: float
) -> tuple[str, pd.DataFrame]:
    """Summarize predeclared robustness requirements by cohort."""
    records = []
    for dataset, frame in results.groupby("dataset", sort=True):
        records.append(
            {
                "dataset": dataset,
                "settings": len(frame),
                "all_degrees_preserved": bool(
                    frame["row_degrees_preserved"].all()
                    and frame["column_degrees_preserved"].all()
                ),
                "all_rmse_points_positive": bool(
                    frame["rmse_improvement_sd"].gt(point_threshold).all()
                ),
                "all_precision_points_positive": bool(
                    frame["precision_improvement"].gt(point_threshold).all()
                ),
                "rmse_ci_pass_rate": float(frame["rmse_ci_lower"].gt(0).mean()),
                "precision_ci_pass_rate": float(
                    frame["precision_ci_lower"].gt(0).mean()
                ),
            }
        )
    summary = pd.DataFrame.from_records(records)
    passed = bool(
        summary["all_degrees_preserved"].all()
        and summary["all_rmse_points_positive"].all()
        and summary["all_precision_points_positive"].all()
        and summary["rmse_ci_pass_rate"].ge(minimum_ci_pass_rate).all()
        and summary["precision_ci_pass_rate"].ge(minimum_ci_pass_rate).all()
    )
    return ("HUMAN_STRUCTURAL_SENSITIVITY_ROBUST" if passed else "SENSITIVITY_WARNING"), summary


def _render_report(
    manifest: dict[str, Any], results: pd.DataFrame, summary: pd.DataFrame
) -> str:
    return "\n".join(
        [
            "# Human structural pathway-score sensitivity",
            "",
            f"**Decision:** `{manifest['decision']}`",
            "",
            "## Adjudication",
            "",
            dataframe_to_markdown(summary.round(4)),
            "",
            "## Complete grid",
            "",
            dataframe_to_markdown(results.round(4)),
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
    st_config = _resolve(config_path, str(config["sources"]["st002081_config"]))
    external_config = _resolve(config_path, str(config["sources"]["st000818_config"]))
    artifact_paths = {
        "st002081": _resolve(config_path, str(config["sources"]["st002081_artifact"])),
        "st000818": _resolve(config_path, str(config["sources"]["st000818_artifact"])),
    }
    integrity = {}
    for key, path in artifact_paths.items():
        _, status = load_verified_manifest(path)
        integrity[key] = status
    if any(value != "verified" for value in integrity.values()):
        raise ValueError(f"Source integrity failed: {integrity}")
    inputs = {
        "ST002081": _st002081_inputs(st_config),
        "ST000818": _st000818_inputs(external_config),
    }
    settings = list(config["settings"])
    tasks = []
    for dataset_index, (dataset_id, (values, groups, masks, base)) in enumerate(inputs.items()):
        for setting_index, setting in enumerate(settings):
            tasks.append(
                (
                    dataset_id,
                    dataset_index * len(settings) + setting_index,
                    setting,
                    values,
                    groups,
                    masks,
                    base,
                )
            )
    jobs = min(int(config["execution"]["parallel_jobs"]), len(tasks))
    with parallel_backend("loky", inner_max_num_threads=1):
        records = Parallel(n_jobs=jobs)(
            delayed(_run_setting)(
                dataset_id,
                setting_index,
                setting,
                values,
                groups,
                masks,
                base,
                int(config["inference"]["bootstrap_draws"]),
                int(config["inference"]["seed"]),
            )
            for dataset_id, setting_index, setting, values, groups, masks, base in tasks
        )
    results = pd.DataFrame.from_records(records).sort_values(["dataset", "setting"])
    decision, summary = adjudicate_human_sensitivity(
        results,
        point_threshold=float(config["gates"]["point_improvement_threshold"]),
        minimum_ci_pass_rate=float(config["gates"]["minimum_ci_pass_rate"]),
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    results_path = output_dir / "sensitivity-grid.csv"
    summary_path = output_dir / "sensitivity-summary.csv"
    report_path = output_dir / "report.md"
    manifest_path = output_dir / "manifest.json"
    results.to_csv(results_path, index=False)
    summary.to_csv(summary_path, index=False)
    manifest: dict[str, Any] = {
        "analysis_id": str(config["analysis_id"]),
        "completed_at": datetime.now(UTC).isoformat(),
        "decision": decision,
        "post_result_sensitivity": True,
        "config": str(config_path),
        "config_sha256": _sha256(config_path),
        "protocol": str(protocol_path),
        "protocol_sha256": _sha256(protocol_path),
        "implementation_sha256": implementation_sha256(),
        "source_integrity": integrity,
        "settings_per_dataset": len(settings),
        "parallel_jobs": jobs,
        "adjudication": summary.to_dict(orient="records"),
        "claim_boundary": str(config["claim_boundary"]),
    }
    report_path.write_text(_render_report(manifest, results, summary))
    manifest["output_sha256"] = {
        path.name: _sha256(path) for path in (results_path, summary_path, report_path)
    }
    manifest_path.write_text(json.dumps(manifest, indent=2, allow_nan=False) + "\n")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config", type=Path, default=Path("config/pathway-score-human-sensitivity.yaml")
    )
    parser.add_argument(
        "--output", type=Path, default=Path("artifacts/pathway-score-human-sensitivity")
    )
    args = parser.parse_args()
    print(json.dumps(run(args.config, args.output), indent=2))


if __name__ == "__main__":
    main()
