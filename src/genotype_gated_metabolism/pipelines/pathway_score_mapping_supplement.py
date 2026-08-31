"""Generate accepted/rejected mapping and QC supplements for pathway-score papers."""

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

from ..analysis.ccle_pathway_prediction import build_signature_matrix, target_feature_sets
from ..analysis.intervention_registry import parse_factor_string
from ..analysis.pathway_scores import (
    disjoint_family_masks,
    lipid_family,
    lipid_structure_descriptors,
    structural_descriptor_incidence,
)
from ..analysis.publication_readiness import load_verified_manifest
from ..datasets.ccle import load_ccle
from ..datasets.metabolomics_workbench import MetabolomicsWorkbenchClient
from ..datasets.mwtab import load_mwtab
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
    return _sha256(Path(inspect.getfile(run)))


def feature_exclusion_reason(
    *, complete: bool, accepted: bool, recognized_family: bool
) -> str:
    """Return one mutually exclusive assay-feature eligibility reason."""
    if accepted:
        return "accepted"
    if not complete:
        return "incomplete_across_eligible_samples"
    if not recognized_family:
        return "unrecognized_or_undersized_lipid_family"
    return "excluded_by_frozen_eligibility"


def _human_feature_table(
    dataset_id: str,
    matrix: pd.DataFrame,
    accepted_features: pd.Index,
    accepted_families: pd.Series,
    masks: pd.DataFrame,
    incidence: pd.DataFrame,
) -> pd.DataFrame:
    accepted = set(accepted_features.astype(str))
    family_by_feature = accepted_families.astype(str).to_dict()
    mask_by_feature = masks.set_index("feature")["mask"].to_dict()
    accepted_family_names = set(accepted_families.astype(str))
    records = []
    for feature in matrix.columns.astype(str):
        complete = bool(matrix[feature].notna().all())
        family = lipid_family(feature)
        is_accepted = feature in accepted
        raw_descriptors = lipid_structure_descriptors(feature)
        filtered = (
            tuple(incidence.columns[incidence.loc[feature]].astype(str))
            if is_accepted and feature in incidence.index
            else ()
        )
        records.append(
            {
                "dataset": dataset_id,
                "feature": feature,
                "family": family_by_feature.get(feature, family),
                "complete": complete,
                "accepted": is_accepted,
                "exclusion_reason": feature_exclusion_reason(
                    complete=complete,
                    accepted=is_accepted,
                    recognized_family=family in accepted_family_names,
                ),
                "mask": mask_by_feature.get(feature),
                "raw_descriptors": ";".join(raw_descriptors),
                "filtered_descriptors": ";".join(filtered),
                "filtered_descriptor_count": len(filtered),
            }
        )
    return pd.DataFrame.from_records(records)


def _st002081_table(base_path: Path, coarse_dir: Path, structural_dir: Path) -> pd.DataFrame:
    config = yaml.safe_load(base_path.read_text())
    dataset = load_mwtab(_resolve(base_path, str(config["source"]["mwtab"])))
    subject_column = str(config["columns"]["subject"])
    excluded = {str(value) for value in config["eligibility"]["exclude_subject_values"]}
    metadata = dataset.sample_metadata
    keep = ~metadata[subject_column].astype(str).isin(excluded)
    matrix = dataset.blocks["metabolomics"].loc[metadata.index[keep]]
    families = pd.read_csv(coarse_dir / "feature-families.csv").set_index("feature")[
        "family"
    ]
    masks = pd.read_csv(coarse_dir / "disjoint-masks.csv")
    incidence = pd.read_csv(structural_dir / "structural-descriptor-incidence.csv").set_index(
        "feature"
    ).astype(bool)
    return _human_feature_table(
        "ST002081", matrix, families.index, families, masks, incidence
    )


def _st000818_table(base_path: Path, artifact_dir: Path) -> pd.DataFrame:
    config = yaml.safe_load(base_path.read_text())
    cache_dir = _resolve(base_path, str(config["source"]["cache_dir"]))
    client = MetabolomicsWorkbenchClient(cache_dir=cache_dir)
    study_id = str(config["source"]["study_id"])
    analysis_id = str(config["source"]["analysis_id"])
    factors = client.factors(study_id)
    factor_name = str(config["columns"]["validation_group_factor"])
    parsed = factors["factors"].map(parse_factor_string)
    groups = pd.Series(
        parsed.map(lambda values: values.get(factor_name)).to_numpy(),
        index=factors["local_sample_id"].astype(str),
    ).dropna()
    measurements = client.measurements(study_id)
    selected = measurements.loc[measurements["analysis_id"].eq(analysis_id)]
    matrix = selected.pivot_table(
        index="local_sample_id", columns="metabolite_name", values="value", aggfunc="mean"
    )
    matrix = matrix.loc[matrix.index.astype(str).intersection(groups.index)]
    accepted_table = pd.read_csv(artifact_dir / "eligible-features.csv")
    families = accepted_table.set_index("feature")["family"]
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
    return _human_feature_table(
        "ST000818", matrix, families.index, families, masks, incidence
    )


def _ccle_tables(
    base_path: Path, artifact_dir: Path
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    config = yaml.safe_load(base_path.read_text())
    mapping = pd.read_csv(_resolve(base_path, str(config["source"]["mapping"])))
    candidates = pd.read_csv(_resolve(base_path, str(config["source"]["candidates"])))
    genes = sorted(
        {
            gene
            for signature in candidates["genes"].dropna().astype(str)
            for gene in parse_gene_signature(signature)
        }
    )
    raw_dir = _resolve(base_path, str(config["source"]["ccle_raw_dir"]))
    dataset = load_ccle(raw_dir, genes)
    mapped = set(mapping["assay_metabolite_id"].astype(str))
    metabolite_records = []
    mapping_by_assay = mapping.set_index("assay_metabolite_id").to_dict(orient="index")
    for assay_id in dataset.blocks["metabolomics"].columns.astype(str):
        source = mapping_by_assay.get(assay_id, {})
        metabolite_records.append(
            {
                "assay_metabolite_id": assay_id,
                "mapping_status": "mapped" if assay_id in mapped else "not_nominated",
                "mapping_reason": (
                    "unique_stable_or_exact_model_mapping"
                    if assay_id in mapped
                    else "absent_from_frozen_unique_stable_mapping_panel"
                ),
                **source,
            }
        )
    metabolite_table = pd.DataFrame.from_records(metabolite_records)
    mapped_metabolites = sorted(mapped)
    signatures = build_signature_matrix(
        dataset.blocks["transcriptomics"],
        candidates["genes"].dropna().astype(str).tolist(),
        aggregation=str(config["model"]["expression_signature_aggregation"]),
    )
    features = target_feature_sets(
        candidates,
        mapped_metabolites=mapped_metabolites,
        available_signatures=list(signatures.columns),
        seed=int(config["model"]["random_feature_seed"]),
    )
    used_signatures = {
        signature for feature in features for signature in feature.network_signatures
    }
    gpr_records = []
    for signature, frame in candidates.groupby("genes", sort=True):
        genes_in_signature = parse_gene_signature(str(signature))
        gpr_records.append(
            {
                "signature": str(signature),
                "genes": ";".join(genes_in_signature),
                "gene_count": len(genes_in_signature),
                "candidate_rows": len(frame),
                "reaction_ids": ";".join(
                    sorted(
                        {
                            value
                            for encoded in frame["reaction_ids"].dropna().astype(str)
                            for value in encoded.split(";")
                            if value
                        }
                    )
                ),
                "subsystems": ";".join(
                    sorted(
                        {
                            value
                            for encoded in frame["subsystems"].dropna().astype(str)
                            for value in encoded.split(";")
                            if value
                        }
                    )
                ),
                "complete_expression_signature": str(signature) in signatures.columns,
                "used_by_prediction_target": str(signature) in used_signatures,
            }
        )
    gpr_table = pd.DataFrame.from_records(gpr_records)
    target_table = pd.read_csv(artifact_dir / "target-feature-sets.csv")
    return metabolite_table, gpr_table, target_table


def _qc_summary(
    human: pd.DataFrame,
    ccle_metabolites: pd.DataFrame,
    gpr: pd.DataFrame,
    targets: pd.DataFrame,
) -> pd.DataFrame:
    records = []
    for dataset, frame in human.groupby("dataset", sort=True):
        records.append(
            {
                "domain": f"{dataset} lipid features",
                "total": len(frame),
                "accepted": int(frame["accepted"].sum()),
                "rejected": int((~frame["accepted"]).sum()),
            }
        )
    records.extend(
        [
            {
                "domain": "CCLE assay metabolites",
                "total": len(ccle_metabolites),
                "accepted": int(ccle_metabolites["mapping_status"].eq("mapped").sum()),
                "rejected": int(
                    (~ccle_metabolites["mapping_status"].eq("mapped")).sum()
                ),
            },
            {
                "domain": "CCLE GPR signatures",
                "total": len(gpr),
                "accepted": int(gpr["used_by_prediction_target"].sum()),
                "rejected": int((~gpr["used_by_prediction_target"]).sum()),
            },
            {
                "domain": "CCLE prediction targets",
                "total": int(targets["target"].nunique()),
                "accepted": int(targets["target"].nunique()),
                "rejected": 0,
            },
        ]
    )
    return pd.DataFrame.from_records(records)


def _render_report(manifest: dict[str, Any], summary: pd.DataFrame) -> str:
    return "\n".join(
        [
            "# Pathway-score mapping and QC supplement",
            "",
            dataframe_to_markdown(summary),
            "",
            (
                "The accompanying tables retain accepted and rejected assay features, explicit "
                "exclusion reasons, structural descriptors, hidden-panel assignments, CCLE "
                "mapping provenance, GPR availability and target-level feature sets."
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
    source = config["sources"]
    st_config = _resolve(config_path, str(source["st002081_config"]))
    st_coarse = _resolve(config_path, str(source["st002081_coarse_artifact"]))
    st_structural = _resolve(config_path, str(source["st002081_structural_artifact"]))
    external_config = _resolve(config_path, str(source["st000818_config"]))
    external_artifact = _resolve(config_path, str(source["st000818_artifact"]))
    ccle_config = _resolve(config_path, str(source["ccle_config"]))
    ccle_artifact = _resolve(config_path, str(source["ccle_artifact"]))
    artifact_paths = {
        "st002081_coarse": st_coarse,
        "st002081_structural": st_structural,
        "st000818": external_artifact,
        "ccle": ccle_artifact,
    }
    integrity = {}
    for key, path in artifact_paths.items():
        _, status = load_verified_manifest(path)
        integrity[key] = status
    if any(status != "verified" for status in integrity.values()):
        raise ValueError(f"Source artifact integrity failed: {integrity}")

    human = pd.concat(
        [
            _st002081_table(st_config, st_coarse, st_structural),
            _st000818_table(external_config, external_artifact),
        ],
        ignore_index=True,
    )
    ccle_metabolites, gpr, targets = _ccle_tables(ccle_config, ccle_artifact)
    summary = _qc_summary(human, ccle_metabolites, gpr, targets)
    output_dir.mkdir(parents=True, exist_ok=True)
    human_path = output_dir / "human-lipid-feature-map.csv"
    ccle_map_path = output_dir / "ccle-metabolite-map.csv"
    gpr_path = output_dir / "ccle-gpr-map.csv"
    target_path = output_dir / "ccle-target-feature-sets.csv"
    summary_path = output_dir / "mapping-qc-summary.csv"
    report_path = output_dir / "report.md"
    manifest_path = output_dir / "manifest.json"
    human.to_csv(human_path, index=False)
    ccle_metabolites.to_csv(ccle_map_path, index=False)
    gpr.to_csv(gpr_path, index=False)
    targets.to_csv(target_path, index=False)
    summary.to_csv(summary_path, index=False)
    manifest: dict[str, Any] = {
        "analysis_id": str(config["analysis_id"]),
        "completed_at": datetime.now(UTC).isoformat(),
        "decision": "MAPPING_QC_SUPPLEMENT_COMPLETE",
        "config": str(config_path),
        "config_sha256": _sha256(config_path),
        "implementation_sha256": implementation_sha256(),
        "source_integrity": integrity,
        "summary": summary.to_dict(orient="records"),
        "claim_boundary": str(config["claim_boundary"]),
    }
    report_path.write_text(_render_report(manifest, summary))
    outputs = [
        human_path,
        ccle_map_path,
        gpr_path,
        target_path,
        summary_path,
        report_path,
    ]
    manifest["output_sha256"] = {path.name: _sha256(path) for path in outputs}
    manifest_path.write_text(json.dumps(manifest, indent=2, allow_nan=False) + "\n")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config", type=Path, default=Path("config/pathway-score-mapping-supplement.yaml")
    )
    parser.add_argument(
        "--output", type=Path, default=Path("artifacts/pathway-score-mapping-supplement")
    )
    args = parser.parse_args()
    print(json.dumps(run(args.config, args.output), indent=2))


if __name__ == "__main__":
    main()
