"""Publication-gate status helpers with checksummed artifact integrity."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_verified_manifest(artifact_dir: Path) -> tuple[dict[str, Any] | None, str]:
    """Load a strict-JSON manifest and fail closed on its declared output checksums."""
    manifest_path = artifact_dir / "manifest.json"
    if not manifest_path.exists():
        return None, "manifest_missing"
    try:
        manifest = json.loads(
            manifest_path.read_text(),
            parse_constant=lambda value: (_ for _ in ()).throw(
                ValueError(f"Non-finite JSON constant: {value}")
            ),
        )
    except (json.JSONDecodeError, ValueError):
        return None, "manifest_invalid_json"
    outputs = manifest.get("output_sha256")
    if not isinstance(outputs, dict) or not outputs:
        return manifest, "output_checksums_missing"
    root = artifact_dir.resolve()
    for filename, expected in outputs.items():
        path = (artifact_dir / filename).resolve()
        if root not in path.parents:
            return manifest, f"output_path_invalid:{filename}"
        if not path.exists():
            return manifest, f"output_missing:{filename}"
        if sha256(path) != expected:
            return manifest, f"checksum_mismatch:{filename}"
    return manifest, "verified"


def gate_record(
    gate_id: str,
    requirement: str,
    status: str,
    result: str,
    evidence_path: Path,
    integrity: str,
) -> dict[str, str]:
    """Build one normalized publication-gate record."""
    if status not in {"PASS", "PARTIAL", "OPEN", "FAIL"}:
        raise ValueError(f"Unknown publication-gate status: {status}")
    return {
        "gate_id": gate_id,
        "requirement": requirement,
        "status": status,
        "result": result,
        "evidence_path": str(evidence_path),
        "artifact_integrity": integrity,
    }
