"""Reusable reader for Metabolomics Workbench mwTab study exports."""

from __future__ import annotations

import hashlib
from pathlib import Path

import pandas as pd

from .base import CrossOmicsDataset


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _pairs(value: str, separator: str) -> dict[str, str]:
    parsed: dict[str, str] = {}
    for item in value.split(separator):
        if ":" in item:
            key, content = item.split(":", 1)
        elif "=" in item:
            key, content = item.split("=", 1)
        else:
            continue
        parsed[key.strip()] = content.strip()
    return parsed


def load_mwtab(path: Path) -> CrossOmicsDataset:
    """Load sample metadata and the named-metabolite matrix from one mwTab file."""
    sample_records: list[dict[str, object]] = []
    matrix_rows: list[list[str]] = []
    in_matrix = False
    with path.open(errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\r\n")
            fields = line.split("\t")
            record_type = fields[0].strip() if fields else ""
            if record_type == "SUBJECT_SAMPLE_FACTORS" and len(fields) >= 5:
                factors = _pairs(fields[3], "|")
                additional = _pairs(fields[4], ";")
                sample_records.append(
                    {
                        "subject_event_id": fields[1].strip(),
                        "sample_id": fields[2].strip(),
                        **{f"factor_{key}": value for key, value in factors.items()},
                        **{f"additional_{key}": value for key, value in additional.items()},
                    }
                )
            elif record_type == "MS_METABOLITE_DATA_START":
                in_matrix = True
            elif record_type == "MS_METABOLITE_DATA_END":
                in_matrix = False
            elif in_matrix and line:
                matrix_rows.append(fields)

    if not sample_records:
        raise ValueError("mwTab file contains no SUBJECT_SAMPLE_FACTORS records.")
    if len(matrix_rows) < 3 or matrix_rows[0][0] != "Samples":
        raise ValueError("mwTab named-metabolite matrix is absent or malformed.")
    sample_ids = matrix_rows[0][1:]
    if len(sample_ids) != len(set(sample_ids)):
        raise ValueError("mwTab matrix contains duplicate sample identifiers.")

    feature_names: list[str] = []
    values: list[list[str | None]] = []
    expected = len(sample_ids)
    for row in matrix_rows[2:]:  # the second matrix line contains redundant factor labels
        if not row or not row[0]:
            continue
        feature_names.append(row[0])
        feature_values = row[1 : expected + 1]
        values.append([*feature_values, *([None] * (expected - len(feature_values)))])
    if len(feature_names) != len(set(feature_names)):
        raise ValueError("mwTab matrix contains duplicate metabolite identifiers.")

    metabolomics = pd.DataFrame(values, index=feature_names, columns=sample_ids).T
    metabolomics = metabolomics.apply(pd.to_numeric, errors="coerce")
    metadata = pd.DataFrame.from_records(sample_records).set_index("sample_id")
    for column in metadata.columns:
        if column.lower().endswith("ageatcollection(month)"):
            metadata[column] = pd.to_numeric(metadata[column], errors="coerce")
    dataset = CrossOmicsDataset(
        sample_metadata=metadata,
        blocks={"metabolomics": metabolomics},
        provenance={
            "format": "Metabolomics Workbench mwTab",
            "path": str(path.resolve()),
            "sha256": _sha256(path),
        },
    )
    return dataset.aligned(["metabolomics"])
