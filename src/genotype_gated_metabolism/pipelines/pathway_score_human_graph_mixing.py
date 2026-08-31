"""Measure graph movement and diversity in the human fixed-degree null ensembles."""

from __future__ import annotations

import argparse
import hashlib
import inspect
import itertools
import json
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import yaml

from ..analysis.pathway_scores import (
    degree_preserving_descriptor_null,
    disjoint_family_masks,
    structural_descriptor_incidence,
)
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


def _edge_jaccard(left: np.ndarray, right: np.ndarray) -> float:
    intersection = int(np.logical_and(left, right).sum())
    union = int(np.logical_or(left, right).sum())
    return float(intersection / union) if union else 1.0


def _audit_dataset(
    dataset: str,
    incidence: pd.DataFrame,
    masks: pd.DataFrame,
    *,
    first_seed: int,
    repeats: int,
    swaps_per_edge: int,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    null_records: list[dict[str, Any]] = []
    mask_records: list[dict[str, Any]] = []
    feature_order = set(incidence.index.astype(str))
    for mask_id in sorted(masks["mask"].unique()):
        hidden = set(
            masks.loc[masks["mask"].eq(mask_id), "feature"].astype(str)
        )
        visible = sorted(feature_order - hidden)
        visible_incidence = incidence.loc[visible]
        usable = visible_incidence.columns[visible_incidence.sum(axis=0).ge(2)]
        observed = visible_incidence.loc[:, usable].astype(bool)
        observed_array = observed.to_numpy(dtype=bool)
        edges = int(observed_array.sum())
        null_arrays: list[np.ndarray] = []
        digests: list[str] = []
        for index in range(repeats):
            seed = first_seed + index + int(mask_id)
            null = degree_preserving_descriptor_null(
                observed, swaps_per_edge=swaps_per_edge, seed=seed
            )
            null_array = null.to_numpy(dtype=bool)
            shared = int(np.logical_and(observed_array, null_array).sum())
            digest = hashlib.sha256(np.packbits(null_array).tobytes()).hexdigest()
            null_arrays.append(null_array)
            digests.append(digest)
            null_records.append(
                {
                    "dataset": dataset,
                    "mask": int(mask_id),
                    "null_index": index,
                    "null_seed": seed,
                    "features": observed.shape[0],
                    "descriptors": observed.shape[1],
                    "edges": edges,
                    "shared_observed_edges": shared,
                    "edge_replacement_fraction": float(1.0 - shared / edges),
                    "observed_null_jaccard": _edge_jaccard(
                        observed_array, null_array
                    ),
                    "row_degrees_preserved": bool(
                        np.array_equal(
                            observed_array.sum(axis=1), null_array.sum(axis=1)
                        )
                    ),
                    "column_degrees_preserved": bool(
                        np.array_equal(
                            observed_array.sum(axis=0), null_array.sum(axis=0)
                        )
                    ),
                    "edge_digest": digest,
                }
            )
        pairwise = [
            _edge_jaccard(null_arrays[left], null_arrays[right])
            for left, right in itertools.combinations(range(repeats), 2)
        ]
        mask_records.append(
            {
                "dataset": dataset,
                "mask": int(mask_id),
                "expected_nulls": repeats,
                "unique_nulls": len(set(digests)),
                "minimum_pairwise_jaccard": float(min(pairwise)),
                "mean_pairwise_jaccard": float(np.mean(pairwise)),
                "maximum_pairwise_jaccard": float(max(pairwise)),
            }
        )
    return null_records, mask_records


def adjudicate_graph_mixing(
    nulls: pd.DataFrame,
    masks: pd.DataFrame,
    *,
    minimum_edge_replacement_fraction: float,
    maximum_pairwise_null_jaccard: float,
) -> tuple[str, dict[str, Any]]:
    """Apply the frozen graph-mixing quality gate."""
    audit = {
        "all_degrees_preserved": bool(
            nulls["row_degrees_preserved"].all()
            and nulls["column_degrees_preserved"].all()
        ),
        "minimum_edge_replacement_fraction": float(
            nulls["edge_replacement_fraction"].min()
        ),
        "edge_replacement_pass": bool(
            nulls["edge_replacement_fraction"]
            .gt(minimum_edge_replacement_fraction)
            .all()
        ),
        "all_nulls_unique": bool(
            masks["unique_nulls"].eq(masks["expected_nulls"]).all()
        ),
        "maximum_pairwise_null_jaccard": float(
            masks["maximum_pairwise_jaccard"].max()
        ),
        "pairwise_diversity_pass": bool(
            masks["maximum_pairwise_jaccard"]
            .lt(maximum_pairwise_null_jaccard)
            .all()
        ),
    }
    return (
        "HUMAN_GRAPH_NULL_MIXING_ADEQUATE"
        if all(
            [
                audit["all_degrees_preserved"],
                audit["edge_replacement_pass"],
                audit["all_nulls_unique"],
                audit["pairwise_diversity_pass"],
            ]
        )
        else "GRAPH_NULL_MIXING_WARNING",
        audit,
    )


def _render_report(
    manifest: dict[str, Any], nulls: pd.DataFrame, masks: pd.DataFrame
) -> str:
    dataset_summary = (
        nulls.groupby("dataset", as_index=False)
        .agg(
            null_panel_rows=("null_seed", "size"),
            minimum_edge_replacement=("edge_replacement_fraction", "min"),
            maximum_observed_null_jaccard=("observed_null_jaccard", "max"),
        )
        .merge(
            masks.groupby("dataset", as_index=False).agg(
                minimum_unique_nulls=("unique_nulls", "min"),
                maximum_pairwise_null_jaccard=("maximum_pairwise_jaccard", "max"),
            ),
            on="dataset",
            validate="one_to_one",
        )
    )
    return "\n".join(
        [
            "# Human graph-null mixing diagnostic",
            "",
            f"**Decision:** `{manifest['decision']}`",
            "",
            "## Dataset summary",
            "",
            dataframe_to_markdown(dataset_summary.round(6)),
            "",
            "## Mask-level null diversity",
            "",
            dataframe_to_markdown(masks.round(6)),
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
    sources = {
        key: _resolve(config_path, str(value))
        for key, value in config["sources"].items()
    }
    _, discovery_integrity = load_verified_manifest(sources["discovery_manifest"])
    _, replication_integrity = load_verified_manifest(sources["replication_manifest"])
    if discovery_integrity != "verified" or replication_integrity != "verified":
        raise ValueError("A contributing human artifact failed integrity verification.")

    discovery_incidence = pd.read_csv(
        sources["discovery_incidence"], index_col="feature"
    ).astype(bool)
    discovery_masks = pd.read_csv(sources["discovery_masks"])

    replication_config = yaml.safe_load(sources["replication_config"].read_text())
    replication_features = pd.read_csv(sources["replication_features"])
    replication_incidence = structural_descriptor_incidence(
        pd.Index(replication_features["feature"].astype(str)),
        minimum_features=int(replication_config["descriptors"]["minimum_features"]),
        maximum_feature_fraction=float(
            replication_config["descriptors"]["maximum_feature_fraction"]
        ),
    )
    replication_families = replication_features.set_index("feature")["family"]
    replication_masks = disjoint_family_masks(
        replication_families,
        masks=int(replication_config["masking"]["masks"]),
        seed=int(replication_config["masking"]["seed"]),
    )

    randomization = config["randomization"]
    discovery_rows, discovery_mask_rows = _audit_dataset(
        "ST002081",
        discovery_incidence,
        discovery_masks,
        first_seed=int(randomization["discovery_first_seed"]),
        repeats=int(randomization["repeats"]),
        swaps_per_edge=int(randomization["swaps_per_edge"]),
    )
    replication_rows, replication_mask_rows = _audit_dataset(
        "ST000818",
        replication_incidence,
        replication_masks,
        first_seed=int(randomization["replication_first_seed"]),
        repeats=int(randomization["repeats"]),
        swaps_per_edge=int(randomization["swaps_per_edge"]),
    )
    nulls = pd.DataFrame.from_records(discovery_rows + replication_rows)
    masks = pd.DataFrame.from_records(discovery_mask_rows + replication_mask_rows)
    decision, adjudication = adjudicate_graph_mixing(
        nulls,
        masks,
        minimum_edge_replacement_fraction=float(
            config["gates"]["minimum_edge_replacement_fraction"]
        ),
        maximum_pairwise_null_jaccard=float(
            config["gates"]["maximum_pairwise_null_jaccard"]
        ),
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    null_path = output_dir / "null-distance-audit.csv"
    mask_path = output_dir / "mask-diversity-summary.csv"
    report_path = output_dir / "report.md"
    manifest_path = output_dir / "manifest.json"
    nulls.to_csv(null_path, index=False)
    masks.to_csv(mask_path, index=False)
    manifest: dict[str, Any] = {
        "analysis_id": str(config["analysis_id"]),
        "completed_at": datetime.now(UTC).isoformat(),
        "decision": decision,
        "post_result_diagnostic": True,
        "config": str(config_path),
        "config_sha256": _sha256(config_path),
        "protocol": str(protocol_path),
        "protocol_sha256": _sha256(protocol_path),
        "implementation_sha256": _sha256(Path(inspect.getfile(run))),
        "source_integrity": {
            "discovery": discovery_integrity,
            "replication": replication_integrity,
        },
        "datasets": int(nulls["dataset"].nunique()),
        "null_panel_rows": len(nulls),
        "adjudication": adjudication,
        "claim_boundary": str(config["claim_boundary"]),
    }
    report_path.write_text(_render_report(manifest, nulls, masks))
    manifest["output_sha256"] = {
        path.name: _sha256(path) for path in (null_path, mask_path, report_path)
    }
    manifest_path.write_text(json.dumps(manifest, indent=2, allow_nan=False) + "\n")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        type=Path,
        default=Path("config/pathway-score-human-graph-mixing.yaml"),
    )
    parser.add_argument(
        "--output", type=Path, default=Path("artifacts/pathway-score-human-graph-mixing")
    )
    args = parser.parse_args()
    print(json.dumps(run(args.config, args.output), indent=2))


if __name__ == "__main__":
    main()
