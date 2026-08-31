from __future__ import annotations

import hashlib
import json

from genotype_gated_metabolism.analysis.publication_readiness import (
    gate_record,
    load_verified_manifest,
)


def test_manifest_output_integrity_is_verified(tmp_path) -> None:
    output = tmp_path / "result.csv"
    output.write_text("a\n1\n")
    checksum = hashlib.sha256(output.read_bytes()).hexdigest()
    (tmp_path / "manifest.json").write_text(
        json.dumps({"output_sha256": {"result.csv": checksum}})
    )

    manifest, status = load_verified_manifest(tmp_path)

    assert manifest is not None
    assert status == "verified"


def test_gate_record_rejects_unknown_status(tmp_path) -> None:
    record = gate_record("g1", "test", "PASS", "ok", tmp_path, "verified")
    assert record["status"] == "PASS"


def test_manifest_without_output_checksums_fails_closed(tmp_path) -> None:
    (tmp_path / "manifest.json").write_text(json.dumps({"analysis": "unchecked"}))
    _, status = load_verified_manifest(tmp_path)
    assert status == "output_checksums_missing"


def test_manifest_with_nan_is_invalid_json(tmp_path) -> None:
    (tmp_path / "manifest.json").write_text('{"value": NaN, "output_sha256": {"x": "y"}}')
    manifest, status = load_verified_manifest(tmp_path)
    assert manifest is None
    assert status == "manifest_invalid_json"


def test_manifest_cannot_hash_path_outside_artifact_directory(tmp_path) -> None:
    (tmp_path / "manifest.json").write_text(
        json.dumps({"output_sha256": {"../outside.csv": "not-used"}})
    )
    _, status = load_verified_manifest(tmp_path)
    assert status == "output_path_invalid:../outside.csv"
