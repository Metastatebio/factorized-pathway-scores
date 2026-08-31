"""Benchmark direct HumanGEM/GPR features for held-out CCLE metabolites."""

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

from ..analysis.ccle_pathway_prediction import (
    bootstrap_target_model_difference,
    build_signature_matrix,
    ccle_pathway_prediction_benchmark,
    target_feature_sets,
)
from ..datasets.ccle import load_ccle
from ..features.signatures import parse_gene_signature
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
    for implementation in (run, ccle_pathway_prediction_benchmark, target_feature_sets):
        path = Path(inspect.getfile(implementation))
        digest.update(path.name.encode())
        digest.update(_sha256(path).encode())
    return digest.hexdigest()


def _render_report(manifest: dict[str, Any], metrics: pd.DataFrame) -> str:
    return "\n".join(
        [
            "# CCLE factorized pathway held-out-metabolite benchmark",
            "",
            f"**Decision:** `{manifest['decision']}`  ",
            "**Human personalized prediction supported:** false  ",
            "**Measured flux supported:** false",
            "",
            "## Result",
            "",
            (
                f"The locked benchmark evaluated {manifest['targets']} target metabolites across "
                f"{manifest['lineages']} lineage groups using two repeats of five "
                "lineage-isolated folds."
            ),
            "",
            dataframe_to_markdown(metrics.round(4)),
            "",
            "## Paired target bootstrap",
            "",
            *[
                (
                    f"- `{item['challenger']}` versus `{item['reference']}`: mean equal-lineage "
                    f"RMSE improvement {item['rmse_improvement_sd']:.4f} SD (95% interval "
                    f"{item['ci_lower']:.4f} to {item['ci_upper']:.4f})."
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
    source = config["source"]
    raw_dir = _resolve(config_path, str(source["ccle_raw_dir"]))
    mapping_path = _resolve(config_path, str(source["mapping"]))
    candidates_path = _resolve(config_path, str(source["candidates"]))
    source_manifest_path = _resolve(config_path, str(source["source_manifest"]))
    mapping = pd.read_csv(mapping_path)
    candidates = pd.read_csv(candidates_path)
    mapped_metabolites = sorted(
        mapping.loc[mapping["mapping_status"].eq("mapped"), "assay_metabolite_id"].astype(str)
    )
    requested_genes = sorted(
        {
            gene
            for signature in candidates["genes"].dropna().astype(str)
            for gene in parse_gene_signature(signature)
        }
    )
    dataset = load_ccle(raw_dir, requested_genes)
    metabolomics = dataset.blocks["metabolomics"].loc[:, mapped_metabolites]
    expression = dataset.blocks["transcriptomics"]
    signatures = build_signature_matrix(
        expression,
        candidates["genes"].dropna().astype(str).tolist(),
        aggregation=str(config["model"]["expression_signature_aggregation"]),
    )
    feature_sets = target_feature_sets(
        candidates,
        mapped_metabolites=mapped_metabolites,
        available_signatures=list(signatures.columns),
        seed=int(config["model"]["random_feature_seed"]),
    )
    model_metrics, target_metrics, feature_table = ccle_pathway_prediction_benchmark(
        metabolomics,
        signatures,
        dataset.sample_metadata["lineage"],
        feature_sets,
        minimum_target_samples=int(config["eligibility"]["minimum_target_samples"]),
        minimum_lineage_samples=int(config["eligibility"]["minimum_lineage_samples"]),
        folds=int(config["model"]["folds"]),
        repeats=int(config["model"]["repeats"]),
        ridge_alpha=float(config["model"]["ridge_alpha"]),
        seed=int(config["model"]["seed"]),
    )
    comparisons = [
        bootstrap_target_model_difference(
            target_metrics,
            reference=reference,
            challenger="factorized_interaction_ridge",
            draws=int(config["inference"]["bootstrap_draws"]),
            seed=int(config["inference"]["seed"]) + index,
        )
        for index, reference in enumerate(
            [
                "population_mean",
                "random_factorized_ridge",
                "network_metabolites_ridge",
                "network_additive_ridge",
                "all_metabolites_ridge",
            ]
        )
    ]
    by_reference = {item["reference"]: item for item in comparisons}
    primary = model_metrics.loc[
        model_metrics["model"].eq("factorized_interaction_ridge")
    ].iloc[0]
    threshold = float(config["gates"]["bootstrap_ci_lower_must_exceed"])
    gate_results = {
        "data_scale": int(primary["targets"]) >= int(config["eligibility"]["minimum_targets"])
        and int(feature_table["lineages"].max())
        >= int(config["eligibility"]["minimum_lineages"]),
        "beats_random_factorization": by_reference["random_factorized_ridge"]["ci_lower"]
        > threshold,
        "beats_network_metabolites": by_reference["network_metabolites_ridge"]["ci_lower"]
        > threshold,
        "positive_median_target_r2": float(primary["median_target_r2"]) > 0,
    }
    decision = (
        "CCLE_FACTORIZED_PATHWAY_SIGNAL"
        if all(gate_results.values())
        else "RESOURCE_ONLY"
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    model_path = output_dir / "model-metrics.csv"
    target_path = output_dir / "target-metrics.csv"
    feature_path = output_dir / "target-feature-sets.csv"
    comparisons_path = output_dir / "bootstrap-comparisons.csv"
    report_path = output_dir / "report.md"
    manifest_path = output_dir / "manifest.json"
    model_metrics.to_csv(model_path, index=False)
    target_metrics.to_csv(target_path, index=False)
    feature_table.to_csv(feature_path, index=False)
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
        "source_sha256": {
            "mapping": _sha256(mapping_path),
            "candidates": _sha256(candidates_path),
            "source_manifest": _sha256(source_manifest_path),
        },
        "aligned_samples": len(metabolomics),
        "mapped_metabolites": len(mapped_metabolites),
        "available_signatures": signatures.shape[1],
        "targets": int(feature_table["target"].nunique()),
        "lineages": int(feature_table["lineages"].max()),
        "model_metrics": model_metrics.replace({np.nan: None}).to_dict(orient="records"),
        "bootstrap_comparisons": comparisons,
        "gate_results": gate_results,
        "claim_boundary": str(config["claim_boundary"]),
    }
    report_path.write_text(_render_report(manifest, model_metrics))
    outputs = [model_path, target_path, feature_path, comparisons_path, report_path]
    manifest["output_sha256"] = {path.name: _sha256(path) for path in outputs}
    manifest_path.write_text(json.dumps(manifest, indent=2, allow_nan=False) + "\n")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=Path("config/pathway-score-ccle.yaml"))
    parser.add_argument("--output", type=Path, default=Path("artifacts/pathway-score-ccle"))
    args = parser.parse_args()
    print(json.dumps(run(args.config, args.output), indent=2))


if __name__ == "__main__":
    main()
