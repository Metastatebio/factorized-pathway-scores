"""Run post-result CCLE reaction-topology and model sensitivity grids."""

from __future__ import annotations

import argparse
import copy
import hashlib
import inspect
import json
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import pandas as pd
import yaml
from joblib import Parallel, delayed, parallel_backend

from ..analysis.ccle_pathway_prediction import (
    bootstrap_target_model_difference,
    build_signature_matrix,
    ccle_pathway_prediction_benchmark,
    target_feature_sets,
)
from ..analysis.publication_readiness import load_verified_manifest
from ..datasets.ccle import load_ccle
from ..features.signatures import parse_gene_signature
from ..reporting import dataframe_to_markdown
from .ccle_proof import _candidate_catalog


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
    for implementation in (
        run,
        ccle_pathway_prediction_benchmark,
        target_feature_sets,
    ):
        path = Path(inspect.getfile(implementation))
        digest.update(path.name.encode())
        digest.update(_sha256(path).encode())
    return digest.hexdigest()


def _catalog_for_setting(
    base_config: dict[str, Any], project_root: Path, setting: dict[str, Any]
) -> tuple[pd.DataFrame, pd.DataFrame]:
    configured = copy.deepcopy(base_config)
    configured["network"]["maximum_metabolites_per_reaction"] = setting[
        "maximum_metabolites_per_reaction"
    ]
    priorities = ["1_direct_reaction"]
    if bool(setting["include_transport"]):
        priorities.append("2_direct_transport")
    configured["candidate_generation"]["included_review_priorities"] = priorities
    return _candidate_catalog(configured, project_root)


def _run_setting(
    setting_index: int,
    setting: dict[str, Any],
    candidates: pd.DataFrame,
    mapped_metabolites: list[str],
    metabolomics: pd.DataFrame,
    expression: pd.DataFrame,
    lineages: pd.Series,
    benchmark: dict[str, Any],
    bootstrap_draws: int,
    bootstrap_seed: int,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    signatures = build_signature_matrix(
        expression,
        candidates["genes"].dropna().astype(str).tolist(),
        aggregation=str(setting["gpr_aggregation"]),
    )
    features = target_feature_sets(
        candidates,
        mapped_metabolites=mapped_metabolites,
        available_signatures=list(signatures.columns),
        seed=int(benchmark["model"]["random_feature_seed"]),
    )
    model_metrics, target_metrics, feature_table = ccle_pathway_prediction_benchmark(
        metabolomics,
        signatures,
        lineages,
        features,
        minimum_target_samples=int(benchmark["eligibility"]["minimum_target_samples"]),
        minimum_lineage_samples=int(benchmark["eligibility"]["minimum_lineage_samples"]),
        folds=int(benchmark["model"]["folds"]),
        repeats=int(benchmark["model"]["repeats"]),
        ridge_alpha=float(setting["ridge_alpha"]),
        seed=int(benchmark["model"]["seed"]),
    )
    comparisons = {}
    for offset, reference in enumerate(
        [
            "random_factorized_ridge",
            "network_metabolites_ridge",
            "network_additive_ridge",
            "all_metabolites_ridge",
        ]
    ):
        comparisons[reference] = bootstrap_target_model_difference(
            target_metrics,
            reference=reference,
            challenger="factorized_interaction_ridge",
            draws=bootstrap_draws,
            seed=bootstrap_seed + setting_index * 10 + offset,
        )
    by_model = model_metrics.set_index("model")
    factorized = by_model.loc["factorized_interaction_ridge"]
    random_comparison = comparisons["random_factorized_ridge"]
    additive_comparison = comparisons["network_additive_ridge"]
    summary = {
        "setting": str(setting["id"]),
        "maximum_metabolites_per_reaction": setting["maximum_metabolites_per_reaction"],
        "include_transport": bool(setting["include_transport"]),
        "ridge_alpha": float(setting["ridge_alpha"]),
        "gpr_aggregation": str(setting["gpr_aggregation"]),
        "candidate_rows": len(candidates),
        "targets": int(feature_table["target"].nunique()),
        "signatures": signatures.shape[1],
        "mean_factorized_rmse_sd": float(factorized["mean_equal_lineage_rmse_sd"]),
        "median_factorized_target_r2": float(factorized["median_target_r2"]),
        "positive_r2_targets": int(factorized["positive_r2_targets"]),
        "rmse_improvement_vs_random_sd": random_comparison["rmse_improvement_sd"],
        "random_ci_lower": random_comparison["ci_lower"],
        "random_ci_upper": random_comparison["ci_upper"],
        "rmse_improvement_vs_additive_sd": additive_comparison["rmse_improvement_sd"],
        "additive_ci_lower": additive_comparison["ci_lower"],
        "additive_ci_upper": additive_comparison["ci_upper"],
    }
    metric_records = model_metrics.assign(setting=str(setting["id"])).to_dict(
        orient="records"
    )
    return summary, metric_records


def adjudicate_ccle_sensitivity(
    summary: pd.DataFrame, *, point_threshold: float, minimum_ci_pass_rate: float
) -> tuple[str, dict[str, Any]]:
    """Apply the frozen post-result reaction-topology robustness rule."""
    adjudication = {
        "settings": len(summary),
        "all_point_improvements_positive": bool(
            summary["rmse_improvement_vs_random_sd"].gt(point_threshold).all()
        ),
        "random_ci_pass_rate": float(summary["random_ci_lower"].gt(0).mean()),
        "positive_median_r2_rate": float(
            summary["median_factorized_target_r2"].gt(0).mean()
        ),
        "interaction_beats_additive_rate": float(
            summary["rmse_improvement_vs_additive_sd"].gt(0).mean()
        ),
    }
    passed = bool(
        adjudication["all_point_improvements_positive"]
        and adjudication["random_ci_pass_rate"] >= minimum_ci_pass_rate
        and adjudication["positive_median_r2_rate"] == 1.0
    )
    return ("CCLE_REACTION_TOPOLOGY_SENSITIVITY_ROBUST" if passed else "SENSITIVITY_WARNING"), adjudication


def _render_report(
    manifest: dict[str, Any], summary: pd.DataFrame, metrics: pd.DataFrame
) -> str:
    return "\n".join(
        [
            "# CCLE reaction-topology sensitivity",
            "",
            f"**Decision:** `{manifest['decision']}`",
            "",
            "## Setting-level results",
            "",
            dataframe_to_markdown(summary.round(4)),
            "",
            "## Complete model metrics",
            "",
            dataframe_to_markdown(metrics.round(4)),
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
    candidate_config_path = _resolve(
        config_path, str(config["sources"]["candidate_config"])
    )
    benchmark_config_path = _resolve(
        config_path, str(config["sources"]["benchmark_config"])
    )
    artifact_path = _resolve(
        config_path, str(config["sources"]["benchmark_artifact"])
    )
    _, integrity = load_verified_manifest(artifact_path)
    if integrity != "verified":
        raise ValueError(f"Primary CCLE artifact integrity failed: {integrity}")
    candidate_base = yaml.safe_load(candidate_config_path.read_text())
    benchmark = yaml.safe_load(benchmark_config_path.read_text())
    project_root = candidate_config_path.parent.parent
    settings = list(config["settings"])
    mapping = None
    catalogs = []
    for setting in settings:
        local_mapping, catalog = _catalog_for_setting(candidate_base, project_root, setting)
        if mapping is None:
            mapping = local_mapping
        catalogs.append(catalog)
    assert mapping is not None
    mapped_metabolites = sorted(
        mapping.loc[mapping["mapping_status"].eq("mapped"), "assay_metabolite_id"].astype(str)
    )
    genes = sorted(
        {
            gene
            for catalog in catalogs
            for signature in catalog["genes"].dropna().astype(str)
            for gene in parse_gene_signature(signature)
        }
    )
    raw_dir = _resolve(benchmark_config_path, str(benchmark["source"]["ccle_raw_dir"]))
    dataset = load_ccle(raw_dir, genes)
    metabolomics = dataset.blocks["metabolomics"].loc[:, mapped_metabolites]
    expression = dataset.blocks["transcriptomics"]
    lineages = dataset.sample_metadata["lineage"]
    jobs = min(int(config["execution"]["parallel_jobs"]), len(settings))
    with parallel_backend("loky", inner_max_num_threads=1):
        outputs = Parallel(n_jobs=jobs)(
            delayed(_run_setting)(
                index,
                setting,
                catalog,
                mapped_metabolites,
                metabolomics,
                expression,
                lineages,
                benchmark,
                int(config["inference"]["bootstrap_draws"]),
                int(config["inference"]["seed"]),
            )
            for index, (setting, catalog) in enumerate(zip(settings, catalogs, strict=True))
        )
    summary = pd.DataFrame.from_records(item[0] for item in outputs)
    metrics = pd.DataFrame.from_records(
        record for _, records in outputs for record in records
    )
    decision, adjudication = adjudicate_ccle_sensitivity(
        summary,
        point_threshold=float(config["gates"]["point_improvement_threshold"]),
        minimum_ci_pass_rate=float(config["gates"]["minimum_ci_pass_rate"]),
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    summary_path = output_dir / "setting-summary.csv"
    metrics_path = output_dir / "model-metrics-grid.csv"
    report_path = output_dir / "report.md"
    manifest_path = output_dir / "manifest.json"
    summary.to_csv(summary_path, index=False)
    metrics.to_csv(metrics_path, index=False)
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
        "primary_artifact_integrity": integrity,
        "candidate_config_sha256": _sha256(candidate_config_path),
        "benchmark_config_sha256": _sha256(benchmark_config_path),
        "settings": len(settings),
        "parallel_jobs": jobs,
        "adjudication": adjudication,
        "claim_boundary": str(config["claim_boundary"]),
    }
    report_path.write_text(_render_report(manifest, summary, metrics))
    manifest["output_sha256"] = {
        path.name: _sha256(path) for path in (summary_path, metrics_path, report_path)
    }
    manifest_path.write_text(json.dumps(manifest, indent=2, allow_nan=False) + "\n")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config", type=Path, default=Path("config/pathway-score-ccle-sensitivity.yaml")
    )
    parser.add_argument(
        "--output", type=Path, default=Path("artifacts/pathway-score-ccle-sensitivity")
    )
    args = parser.parse_args()
    print(json.dumps(run(args.config, args.output), indent=2))


if __name__ == "__main__":
    main()
