"""Run a multi-seed random-feature null ensemble for the CCLE pathway benchmark."""

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

from ..analysis.ccle_pathway_prediction import (
    bootstrap_target_model_difference,
    build_signature_matrix,
    ccle_pathway_prediction_benchmark,
    propensity_matched_target_feature_sets,
    target_feature_sets,
)
from ..analysis.publication_readiness import load_verified_manifest
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
    for implementation in (
        run,
        ccle_pathway_prediction_benchmark,
        target_feature_sets,
        propensity_matched_target_feature_sets,
    ):
        path = Path(inspect.getfile(implementation))
        digest.update(path.name.encode())
        digest.update(_sha256(path).encode())
    return digest.hexdigest()


def _run_seed(
    seed_index: int,
    random_seed: int,
    candidates: pd.DataFrame,
    mapped_metabolites: list[str],
    metabolomics: pd.DataFrame,
    signatures: pd.DataFrame,
    lineages: pd.Series,
    benchmark: dict[str, Any],
    bootstrap_draws: int,
    bootstrap_seed: int,
    matching_mode: str,
) -> tuple[dict[str, Any], pd.DataFrame]:
    if matching_mode == "network_degree_and_coverage":
        features = propensity_matched_target_feature_sets(
            candidates,
            mapped_metabolites=mapped_metabolites,
            available_signatures=list(signatures.columns),
            metabolite_coverage=metabolomics.notna().mean(axis=0),
            seed=random_seed,
        )
    elif matching_mode == "dimension_only":
        features = target_feature_sets(
            candidates,
            mapped_metabolites=mapped_metabolites,
            available_signatures=list(signatures.columns),
            seed=random_seed,
        )
    else:
        raise ValueError(f"Unknown random-feature matching mode: {matching_mode}")
    model_metrics, target_metrics, feature_table = ccle_pathway_prediction_benchmark(
        metabolomics,
        signatures,
        lineages,
        features,
        minimum_target_samples=int(benchmark["eligibility"]["minimum_target_samples"]),
        minimum_lineage_samples=int(benchmark["eligibility"]["minimum_lineage_samples"]),
        folds=int(benchmark["model"]["folds"]),
        repeats=int(benchmark["model"]["repeats"]),
        ridge_alpha=float(benchmark["model"]["ridge_alpha"]),
        seed=int(benchmark["model"]["seed"]),
    )
    comparison = bootstrap_target_model_difference(
        target_metrics,
        reference="random_factorized_ridge",
        challenger="factorized_interaction_ridge",
        draws=bootstrap_draws,
        seed=bootstrap_seed + seed_index,
    )
    by_model = model_metrics.set_index("model")
    targets = sorted(feature_table["target"].astype(str).unique())
    target_digest = hashlib.sha256("\n".join(targets).encode()).hexdigest()
    balance = _feature_balance(
        features,
        candidates=candidates,
        mapped_metabolites=mapped_metabolites,
        available_signatures=list(signatures.columns),
        metabolite_coverage=metabolomics.notna().mean(axis=0),
    )
    summary = {
        "random_feature_seed": random_seed,
        "targets": len(targets),
        "target_digest": target_digest,
        "random_rmse_sd": float(
            by_model.loc["random_factorized_ridge", "mean_equal_lineage_rmse_sd"]
        ),
        "factorized_rmse_sd": float(
            by_model.loc["factorized_interaction_ridge", "mean_equal_lineage_rmse_sd"]
        ),
        "improvement_sd": float(comparison["rmse_improvement_sd"]),
        "ci_lower": float(comparison["ci_lower"]),
        "ci_upper": float(comparison["ci_upper"]),
        **balance,
    }
    paired = target_metrics.loc[
        target_metrics["model"].isin(
            ["random_factorized_ridge", "factorized_interaction_ridge"]
        ),
        ["target", "model", "equal_lineage_weighted_rmse_sd"],
    ].pivot(
        index="target", columns="model", values="equal_lineage_weighted_rmse_sd"
    )
    paired = paired.rename(
        columns={
            "random_factorized_ridge": "random_rmse_sd",
            "factorized_interaction_ridge": "factorized_rmse_sd",
        }
    ).reset_index()
    paired.insert(0, "random_feature_seed", random_seed)
    return summary, paired


def _feature_balance(
    features: list[Any],
    *,
    candidates: pd.DataFrame,
    mapped_metabolites: list[str],
    available_signatures: list[str],
    metabolite_coverage: pd.Series,
) -> dict[str, float]:
    metabolite_set = set(mapped_metabolites)
    signature_set = set(available_signatures)
    neighbors: dict[str, set[str]] = {value: set() for value in metabolite_set}
    signature_degree = {value: 0 for value in signature_set}
    for row in candidates.to_dict(orient="records"):
        left, right, signature = (
            str(row["metabolite_a"]),
            str(row["metabolite_b"]),
            str(row["genes"]),
        )
        if left not in metabolite_set or right not in metabolite_set:
            continue
        if signature not in signature_set:
            continue
        neighbors[left].add(right)
        neighbors[right].add(left)
        signature_degree[signature] += 1
    metabolite_degree = {
        value: float(np.log1p(len(neighbors[value]))) for value in metabolite_set
    }
    signature_log_degree = {
        value: float(np.log1p(signature_degree[value])) for value in signature_set
    }
    metabolite_degree_differences = []
    coverage_differences = []
    signature_degree_differences = []
    for feature in features:
        metabolite_degree_differences.append(
            abs(
                np.mean([metabolite_degree[value] for value in feature.network_metabolites])
                - np.mean([metabolite_degree[value] for value in feature.random_metabolites])
            )
        )
        coverage_differences.append(
            abs(
                np.mean([metabolite_coverage[value] for value in feature.network_metabolites])
                - np.mean([metabolite_coverage[value] for value in feature.random_metabolites])
            )
        )
        signature_degree_differences.append(
            abs(
                np.mean([signature_log_degree[value] for value in feature.network_signatures])
                - np.mean([signature_log_degree[value] for value in feature.random_signatures])
            )
        )
    return {
        "mean_metabolite_log_degree_imbalance": float(
            np.mean(metabolite_degree_differences)
        ),
        "mean_metabolite_coverage_imbalance": float(np.mean(coverage_differences)),
        "mean_signature_log_degree_imbalance": float(
            np.mean(signature_degree_differences)
        ),
    }


def _expected_null_comparison(
    target_ensemble: pd.DataFrame, *, draws: int, seed: int
) -> dict[str, float | int]:
    per_target = (
        target_ensemble.groupby("target", as_index=False)
        .agg(
            expected_random_rmse_sd=("random_rmse_sd", "mean"),
            factorized_rmse_sd=("factorized_rmse_sd", "mean"),
            factorized_range_sd=(
                "factorized_rmse_sd",
                lambda values: float(values.max() - values.min()),
            ),
        )
        .sort_values("target")
    )
    differences = (
        per_target["expected_random_rmse_sd"] - per_target["factorized_rmse_sd"]
    ).to_numpy(dtype=float)
    rng = np.random.default_rng(seed)
    indexes = rng.integers(0, len(differences), size=(draws, len(differences)))
    bootstrap = differences[indexes].mean(axis=1)
    return {
        "targets": len(differences),
        "draws": draws,
        "improvement_sd": float(differences.mean()),
        "ci_lower": float(np.quantile(bootstrap, 0.025)),
        "ci_upper": float(np.quantile(bootstrap, 0.975)),
        "maximum_factorized_target_range_sd": float(
            per_target["factorized_range_sd"].max()
        ),
    }


def adjudicate_null_ensemble(
    seed_summary: pd.DataFrame,
    *,
    ensemble_ci_lower: float,
    minimum_seed_ci_pass_rate: float,
    improvement_threshold: float,
    invariant_tolerance: float,
) -> tuple[str, dict[str, Any]]:
    """Apply the frozen stochastic-null robustness gate."""
    audit = {
        "seeds": len(seed_summary),
        "target_sets_invariant": bool(seed_summary["target_digest"].nunique() == 1),
        "target_counts_invariant": bool(seed_summary["targets"].nunique() == 1),
        "factorized_rmse_range_sd": float(
            seed_summary["factorized_rmse_sd"].max()
            - seed_summary["factorized_rmse_sd"].min()
        ),
        "all_seed_point_improvements_positive": bool(
            seed_summary["improvement_sd"].gt(improvement_threshold).all()
        ),
        "seed_ci_pass_rate": float(seed_summary["ci_lower"].gt(0).mean()),
        "ensemble_ci_lower": float(ensemble_ci_lower),
    }
    passed = bool(
        audit["target_sets_invariant"]
        and audit["target_counts_invariant"]
        and audit["factorized_rmse_range_sd"] <= invariant_tolerance
        and audit["all_seed_point_improvements_positive"]
        and audit["seed_ci_pass_rate"] >= minimum_seed_ci_pass_rate
        and ensemble_ci_lower > improvement_threshold
    )
    return (
        "CCLE_RANDOM_FEATURE_NULL_ENSEMBLE_ROBUST"
        if passed
        else "NULL_ENSEMBLE_WARNING",
        audit,
    )


def _render_report(
    manifest: dict[str, Any], seed_summary: pd.DataFrame, ensemble: dict[str, Any]
) -> str:
    return "\n".join(
        [
            "# CCLE random-feature null ensemble",
            "",
            f"**Decision:** `{manifest['decision']}`",
            "",
            "## Expected-null contrast",
            "",
            (
                f"Across {manifest['random_feature_seeds']} random feature realizations, the "
                f"expected random-minus-topology RMSE improvement was "
                f"{ensemble['improvement_sd']:.4f} SD (95% target-bootstrap interval "
                f"{ensemble['ci_lower']:.4f} to {ensemble['ci_upper']:.4f})."
            ),
            "",
            "## Complete seed results",
            "",
            dataframe_to_markdown(seed_summary.round(6)),
            "",
            "## Adjudication",
            "",
            dataframe_to_markdown(pd.DataFrame([manifest["adjudication"]]).round(6)),
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
    benchmark_config_path = _resolve(
        config_path, str(config["sources"]["benchmark_config"])
    )
    benchmark_artifact = _resolve(
        config_path, str(config["sources"]["benchmark_artifact"])
    )
    _, integrity = load_verified_manifest(benchmark_artifact)
    if integrity != "verified":
        raise ValueError(f"Primary CCLE artifact integrity failed: {integrity}")
    benchmark = yaml.safe_load(benchmark_config_path.read_text())
    mapping_path = _resolve(benchmark_config_path, str(benchmark["source"]["mapping"]))
    candidates_path = _resolve(
        benchmark_config_path, str(benchmark["source"]["candidates"])
    )
    mapping = pd.read_csv(mapping_path)
    candidates = pd.read_csv(candidates_path)
    mapped_metabolites = sorted(
        mapping.loc[
            mapping["mapping_status"].eq("mapped"), "assay_metabolite_id"
        ].astype(str)
    )
    genes = sorted(
        {
            gene
            for signature in candidates["genes"].dropna().astype(str)
            for gene in parse_gene_signature(signature)
        }
    )
    raw_dir = _resolve(
        benchmark_config_path, str(benchmark["source"]["ccle_raw_dir"])
    )
    dataset = load_ccle(raw_dir, genes)
    metabolomics = dataset.blocks["metabolomics"].loc[:, mapped_metabolites]
    expression = dataset.blocks["transcriptomics"]
    signatures = build_signature_matrix(
        expression,
        candidates["genes"].dropna().astype(str).tolist(),
        aggregation=str(benchmark["model"]["expression_signature_aggregation"]),
    )
    seeds = [int(seed) for seed in config["random_feature_seeds"]]
    matching_mode = str(config.get("matching", {}).get("mode", "dimension_only"))
    jobs = min(int(config["execution"]["parallel_jobs"]), len(seeds))
    with parallel_backend("loky", inner_max_num_threads=1):
        outputs = Parallel(n_jobs=jobs)(
            delayed(_run_seed)(
                index,
                random_seed,
                candidates,
                mapped_metabolites,
                metabolomics,
                signatures,
                dataset.sample_metadata["lineage"],
                benchmark,
                int(config["inference"]["per_seed_bootstrap_draws"]),
                int(config["inference"]["seed"]),
                matching_mode,
            )
            for index, random_seed in enumerate(seeds)
        )
    seed_summary = pd.DataFrame.from_records(item[0] for item in outputs).sort_values(
        "random_feature_seed"
    )
    target_ensemble = pd.concat((item[1] for item in outputs), ignore_index=True)
    target_ensemble = target_ensemble.sort_values(
        ["random_feature_seed", "target"]
    ).reset_index(drop=True)
    expected_null = _expected_null_comparison(
        target_ensemble,
        draws=int(config["inference"]["ensemble_bootstrap_draws"]),
        seed=int(config["inference"]["seed"]) + len(seeds),
    )
    decision, adjudication = adjudicate_null_ensemble(
        seed_summary,
        ensemble_ci_lower=float(expected_null["ci_lower"]),
        minimum_seed_ci_pass_rate=float(
            config["gates"]["minimum_seed_ci_pass_rate"]
        ),
        improvement_threshold=float(config["gates"]["improvement_threshold"]),
        invariant_tolerance=float(config["gates"]["invariant_tolerance"]),
    )
    balance_gates = {
        "metabolite_log_degree_balance": bool(
            seed_summary["mean_metabolite_log_degree_imbalance"].max()
            <= float(
                config["gates"].get(
                    "maximum_metabolite_log_degree_imbalance", float("inf")
                )
            )
        ),
        "metabolite_coverage_balance": bool(
            seed_summary["mean_metabolite_coverage_imbalance"].max()
            <= float(
                config["gates"].get(
                    "maximum_metabolite_coverage_imbalance", float("inf")
                )
            )
        ),
        "signature_log_degree_balance": bool(
            seed_summary["mean_signature_log_degree_imbalance"].max()
            <= float(
                config["gates"].get(
                    "maximum_signature_log_degree_imbalance", float("inf")
                )
            )
        ),
    }
    adjudication["matching_mode"] = matching_mode
    adjudication.update(balance_gates)
    if not all(balance_gates.values()):
        decision = "NULL_ENSEMBLE_WARNING"
    output_dir.mkdir(parents=True, exist_ok=True)
    seed_path = output_dir / "seed-summary.csv"
    target_path = output_dir / "target-null-ensemble.csv"
    report_path = output_dir / "report.md"
    manifest_path = output_dir / "manifest.json"
    seed_summary.to_csv(seed_path, index=False)
    target_ensemble.to_csv(target_path, index=False)
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
        "benchmark_config_sha256": _sha256(benchmark_config_path),
        "mapping_sha256": _sha256(mapping_path),
        "candidates_sha256": _sha256(candidates_path),
        "random_feature_seeds": len(seeds),
        "matching_mode": matching_mode,
        "parallel_jobs": jobs,
        "expected_null_comparison": expected_null,
        "adjudication": adjudication,
        "claim_boundary": str(config["claim_boundary"]),
    }
    report_path.write_text(_render_report(manifest, seed_summary, expected_null))
    manifest["output_sha256"] = {
        path.name: _sha256(path) for path in (seed_path, target_path, report_path)
    }
    manifest_path.write_text(json.dumps(manifest, indent=2, allow_nan=False) + "\n")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        type=Path,
        default=Path("config/pathway-score-ccle-null-ensemble.yaml"),
    )
    parser.add_argument(
        "--output", type=Path, default=Path("artifacts/pathway-score-ccle-null-ensemble")
    )
    args = parser.parse_args()
    print(json.dumps(run(args.config, args.output), indent=2))


if __name__ == "__main__":
    main()
