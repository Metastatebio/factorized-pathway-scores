"""Map assay metabolites to compartment-independent Human-GEM identifiers."""

from __future__ import annotations

import re
import unicodedata
from collections import defaultdict
from collections.abc import Iterable, Mapping
from pathlib import Path

import pandas as pd


def _identifiers(value: object) -> set[str]:
    if pd.isna(value):
        return set()
    return {item.strip().upper() for item in str(value).split(";") if item.strip()}


def normalize_metabolite_name(value: object) -> str:
    """Normalize punctuation without erasing biochemical qualifiers or chirality."""
    if pd.isna(value):
        return ""
    text = unicodedata.normalize("NFKD", str(value)).lower()
    greek = {
        "α": "alpha",
        "β": "beta",
        "γ": "gamma",
        "δ": "delta",
    }
    for symbol, name in greek.items():
        text = text.replace(symbol, name)
    return re.sub(r"[^a-z0-9]+", "", text)


def expand_panel_by_exact_model_names(
    curated_panel: pd.DataFrame,
    assay_metabolite_ids: Iterable[str],
    model_metabolite_names: Mapping[str, str],
    annotations: pd.DataFrame,
) -> pd.DataFrame:
    """Add unique exact normalized-name matches with stable model identifiers.

    Curated rows always win. Automatic rows are admitted only when the normalized
    assay name and normalized Human-GEM name are each unique and the matched model
    metabolite has an HMDB or ChEBI annotation. Stable-ID mapping remains a separate
    validation step in :func:`map_metabolite_panel`.
    """
    required = {"assay_metabolite_id", "display_name", "hmdb_id", "chebi_id"}
    missing = required - set(curated_panel.columns)
    if missing:
        raise ValueError(f"Assay panel is missing: {', '.join(sorted(missing))}")
    if curated_panel["assay_metabolite_id"].duplicated().any():
        raise ValueError("assay_metabolite_id values must be unique.")

    annotation_required = {"metsNoComp", "metHMDBID", "metChEBIID"}
    annotation_missing = annotation_required - set(annotations.columns)
    if annotation_missing:
        raise ValueError(
            "Human-GEM annotations are missing: " + ", ".join(sorted(annotation_missing))
        )

    assay_by_name: dict[str, set[str]] = defaultdict(set)
    for assay_id in assay_metabolite_ids:
        normalized = normalize_metabolite_name(assay_id)
        if normalized:
            assay_by_name[normalized].add(str(assay_id))

    model_by_name: dict[str, set[str]] = defaultdict(set)
    for base_id, name in model_metabolite_names.items():
        normalized = normalize_metabolite_name(name)
        if normalized:
            model_by_name[normalized].add(str(base_id))

    stable_ids: dict[str, dict[str, set[str]]] = defaultdict(
        lambda: {"hmdb": set(), "chebi": set()}
    )
    for row in annotations.to_dict(orient="records"):
        base_id = str(row["metsNoComp"])
        stable_ids[base_id]["hmdb"].update(_identifiers(row["metHMDBID"]))
        stable_ids[base_id]["chebi"].update(_identifiers(row["metChEBIID"]))

    curated = curated_panel.copy()
    curated["mapping_origin"] = "curated_stable_id"
    curated_ids = set(curated["assay_metabolite_id"].astype(str))
    records: list[dict[str, str]] = []
    for normalized in sorted(set(assay_by_name) & set(model_by_name)):
        assay_matches = assay_by_name[normalized]
        model_matches = model_by_name[normalized]
        if len(assay_matches) != 1 or len(model_matches) != 1:
            continue
        assay_id = next(iter(assay_matches))
        base_id = next(iter(model_matches))
        if assay_id in curated_ids:
            continue
        identifiers = stable_ids[base_id]
        if not identifiers["hmdb"] and not identifiers["chebi"]:
            continue
        records.append(
            {
                "assay_metabolite_id": assay_id,
                "display_name": str(model_metabolite_names[base_id]),
                "hmdb_id": ";".join(sorted(identifiers["hmdb"])),
                "chebi_id": ";".join(sorted(identifiers["chebi"])),
                "mapping_origin": "unique_exact_normalized_name",
            }
        )
    automatic = pd.DataFrame.from_records(
        records,
        columns=[
            "assay_metabolite_id",
            "display_name",
            "hmdb_id",
            "chebi_id",
            "mapping_origin",
        ],
    )
    return pd.concat([curated, automatic], ignore_index=True).sort_values(
        ["mapping_origin", "assay_metabolite_id"], kind="stable"
    )


def build_identifier_indexes(
    annotations: pd.DataFrame,
) -> tuple[dict[str, set[str]], dict[str, set[str]]]:
    """Return external-ID to compartment-independent Human-GEM ID indexes."""
    required = {"metsNoComp", "metHMDBID", "metChEBIID"}
    missing = required - set(annotations.columns)
    if missing:
        raise ValueError(f"Human-GEM annotations are missing: {', '.join(sorted(missing))}")

    hmdb_index: dict[str, set[str]] = defaultdict(set)
    chebi_index: dict[str, set[str]] = defaultdict(set)
    for row in annotations.to_dict(orient="records"):
        base_id = str(row["metsNoComp"])
        for identifier in _identifiers(row["metHMDBID"]):
            hmdb_index[identifier].add(base_id)
        for identifier in _identifiers(row["metChEBIID"]):
            chebi_index[identifier].add(base_id)
    return dict(hmdb_index), dict(chebi_index)


def map_metabolite_panel(panel: pd.DataFrame, annotations: pd.DataFrame) -> pd.DataFrame:
    """Map by stable IDs only; report ambiguity or cross-database conflicts."""
    required = {"assay_metabolite_id", "display_name", "hmdb_id", "chebi_id"}
    missing = required - set(panel.columns)
    if missing:
        raise ValueError(f"Assay panel is missing: {', '.join(sorted(missing))}")
    if panel["assay_metabolite_id"].duplicated().any():
        raise ValueError("assay_metabolite_id values must be unique.")

    hmdb_index, chebi_index = build_identifier_indexes(annotations)
    records: list[dict[str, object]] = []
    for row in panel.to_dict(orient="records"):
        hmdb_ids = _identifiers(row["hmdb_id"])
        chebi_ids = _identifiers(row["chebi_id"])
        hmdb_matches = set().union(*(hmdb_index.get(item, set()) for item in hmdb_ids))
        chebi_matches = set().union(*(chebi_index.get(item, set()) for item in chebi_ids))

        if hmdb_matches and chebi_matches and not hmdb_matches & chebi_matches:
            candidates = hmdb_matches | chebi_matches
            status = "conflicting_identifiers"
            matched_by = "HMDB;ChEBI"
        else:
            candidates = (hmdb_matches & chebi_matches) or hmdb_matches or chebi_matches
            matched_by_parts = []
            if hmdb_matches:
                matched_by_parts.append("HMDB")
            if chebi_matches:
                matched_by_parts.append("ChEBI")
            matched_by = ";".join(matched_by_parts)
            if len(candidates) == 1:
                status = "mapped"
            elif len(candidates) > 1:
                status = "ambiguous"
            else:
                status = "unmapped"

        base_id = next(iter(candidates)) if status == "mapped" else ""
        records.append(
            {
                "assay_metabolite_id": row["assay_metabolite_id"],
                "human_gem_base_id": base_id,
                "display_name": row["display_name"],
                "mapping_status": status,
                "matched_by": matched_by,
                "candidate_human_gem_base_ids": ";".join(sorted(candidates)),
                "hmdb_id": row["hmdb_id"],
                "chebi_id": row["chebi_id"],
                "mapping_origin": row.get("mapping_origin", "unspecified"),
            }
        )
    return pd.DataFrame.from_records(records)


def read_human_gem_metabolites(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Missing Human-GEM metabolite annotations: {path}")
    return pd.read_csv(path, sep="\t", dtype=str)
