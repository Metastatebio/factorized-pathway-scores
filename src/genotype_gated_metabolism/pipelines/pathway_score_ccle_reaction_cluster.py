"""Audit CCLE topology effects under reaction-subsystem clustered uncertainty."""

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

from ..analysis.publication_readiness import load_verified_manifest
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


def _dominant_subsystem(candidates: pd.DataFrame, target: str) -> str:
    block = candidates.loc[
        candidates["metabolite_a"].astype(str).eq(target)
        | candidates["metabolite_b"].astype(str).eq(target)
    ]
    counts: dict[str, int] = {}
    for value in block["subsystems"].dropna().astype(str):
        for subsystem in value.split(";"):
            subsystem = subsystem.strip()
            if subsystem:
                counts[subsystem] = counts.get(subsystem, 0) + 1
    if not counts:
        return "Unassigned"
    return min(counts, key=lambda value: (-counts[value], value))


def _cluster_bootstrap(
    effects: pd.DataFrame, *, draws: int, seed: int
) -> dict[str, float | int]:
    groups = [
        block["effect_sd"].to_numpy(dtype=float)
        for _, block in effects.groupby("subsystem", sort=True)
    ]
    rng = np.random.default_rng(seed)
    values = np.empty(draws, dtype=float)
    for draw in range(draws):
        sampled = rng.integers(0, len(groups), size=len(groups))
        values[draw] = float(
            np.concatenate([groups[int(index)] for index in sampled]).mean()
        )
    return {
        "clusters": len(groups),
        "targets": len(effects),
        "draws": draws,
        "effect_sd": float(effects["effect_sd"].mean()),
        "ci_lower": float(np.quantile(values, 0.025)),
        "ci_upper": float(np.quantile(values, 0.975)),
    }


def adjudicate_cluster_robustness(
    *,
    clusters: int,
    ci_lower: float,
    minimum_leave_one_cluster_effect: float,
    minimum_clusters: int,
    threshold: float,
) -> tuple[str, dict[str, bool | float | int]]:
    """Apply the frozen cluster-robustness gate."""
    audit: dict[str, bool | float | int] = {
        "clusters": clusters,
        "minimum_clusters_pass": clusters >= minimum_clusters,
        "cluster_ci_lower": ci_lower,
        "cluster_ci_pass": ci_lower > threshold,
        "minimum_leave_one_cluster_effect_sd": minimum_leave_one_cluster_effect,
        "leave_one_cluster_pass": minimum_leave_one_cluster_effect > threshold,
    }
    passed = bool(
        audit["minimum_clusters_pass"]
        and audit["cluster_ci_pass"]
        and audit["leave_one_cluster_pass"]
    )
    return (
        "CCLE_REACTION_CLUSTER_ROBUST" if passed else "REACTION_CLUSTER_WARNING",
        audit,
    )


def _render_report(
    manifest: dict[str, Any],
    inference: dict[str, Any],
    cluster_summary: pd.DataFrame,
) -> str:
    return "\n".join(
        [
            "# CCLE reaction-cluster robustness",
            "",
            f"**Decision:** `{manifest['decision']}`",
            "",
            "## Cluster-resampled inference",
            "",
            (
                f"The property-matched expected-null effect was {inference['effect_sd']:.4f} SD "
                f"across {inference['targets']} targets assigned to {inference['clusters']} "
                f"subsystems. The 95% cluster-bootstrap interval was "
                f"{inference['ci_lower']:.4f} to {inference['ci_upper']:.4f}."
            ),
            "",
            "## Subsystem summary",
            "",
            dataframe_to_markdown(cluster_summary.round(6)),
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
    null_artifact = _resolve(config_path, str(config["sources"]["null_artifact"]))
    candidates_path = _resolve(config_path, str(config["sources"]["candidates"]))
    null_manifest, integrity = load_verified_manifest(null_artifact)
    if integrity != "verified" or null_manifest is None:
        raise ValueError(f"Property-matched null integrity failed: {integrity}")
    target_ensemble = pd.read_csv(null_artifact / "target-null-ensemble.csv")
    candidates = pd.read_csv(candidates_path)
    effects = (
        target_ensemble.groupby("target", as_index=False)
        .agg(
            expected_random_rmse_sd=("random_rmse_sd", "mean"),
            factorized_rmse_sd=("factorized_rmse_sd", "mean"),
        )
        .sort_values("target")
    )
    effects["effect_sd"] = (
        effects["expected_random_rmse_sd"] - effects["factorized_rmse_sd"]
    )
    effects["subsystem"] = [
        _dominant_subsystem(candidates, str(target)) for target in effects["target"]
    ]
    inference = _cluster_bootstrap(
        effects,
        draws=int(config["inference"]["bootstrap_draws"]),
        seed=int(config["inference"]["seed"]),
    )
    cluster_summary = (
        effects.groupby("subsystem", as_index=False)
        .agg(targets=("target", "size"), mean_effect_sd=("effect_sd", "mean"))
        .sort_values("mean_effect_sd")
    )
    leave_one_out = []
    for subsystem in sorted(effects["subsystem"].unique()):
        retained = effects.loc[~effects["subsystem"].eq(subsystem), "effect_sd"]
        leave_one_out.append(
            {"excluded_subsystem": subsystem, "effect_sd": float(retained.mean())}
        )
    leave_one_out_table = pd.DataFrame.from_records(leave_one_out).sort_values(
        "effect_sd"
    )
    minimum_leave_one_out = float(leave_one_out_table["effect_sd"].min())
    decision, adjudication = adjudicate_cluster_robustness(
        clusters=int(inference["clusters"]),
        ci_lower=float(inference["ci_lower"]),
        minimum_leave_one_cluster_effect=minimum_leave_one_out,
        minimum_clusters=int(config["gates"]["minimum_clusters"]),
        threshold=float(config["gates"]["lower_bound_must_exceed"]),
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    effects_path = output_dir / "target-effects-with-subsystem.csv"
    cluster_path = output_dir / "subsystem-summary.csv"
    leave_out_path = output_dir / "leave-one-subsystem-out.csv"
    report_path = output_dir / "report.md"
    manifest_path = output_dir / "manifest.json"
    effects.to_csv(effects_path, index=False)
    cluster_summary.to_csv(cluster_path, index=False)
    leave_one_out_table.to_csv(leave_out_path, index=False)
    manifest: dict[str, Any] = {
        "analysis_id": str(config["analysis_id"]),
        "completed_at": datetime.now(UTC).isoformat(),
        "decision": decision,
        "post_result_sensitivity": True,
        "config": str(config_path),
        "config_sha256": _sha256(config_path),
        "protocol": str(protocol_path),
        "protocol_sha256": _sha256(protocol_path),
        "implementation_sha256": _sha256(Path(inspect.getfile(run))),
        "null_artifact_integrity": integrity,
        "null_manifest_sha256": _sha256(null_artifact / "manifest.json"),
        "candidates_sha256": _sha256(candidates_path),
        "inference": inference,
        "adjudication": adjudication,
        "claim_boundary": str(config["claim_boundary"]),
    }
    report_path.write_text(_render_report(manifest, inference, cluster_summary))
    manifest["output_sha256"] = {
        path.name: _sha256(path)
        for path in (effects_path, cluster_path, leave_out_path, report_path)
    }
    manifest_path.write_text(json.dumps(manifest, indent=2, allow_nan=False) + "\n")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        type=Path,
        default=Path("config/pathway-score-ccle-reaction-cluster.yaml"),
    )
    parser.add_argument(
        "--output", type=Path, default=Path("artifacts/pathway-score-ccle-reaction-cluster")
    )
    args = parser.parse_args()
    print(json.dumps(run(args.config, args.output), indent=2))


if __name__ == "__main__":
    main()
