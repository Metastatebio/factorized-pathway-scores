"""Verify released artifacts and principal manuscript-facing invariants."""

from __future__ import annotations

import json
import subprocess
from pathlib import Path

from genotype_gated_metabolism.analysis.publication_readiness import load_verified_manifest

ROOT = Path(__file__).resolve().parents[1]
CORE_ARTIFACTS = [
    "pathway-score-st002081",
    "pathway-score-st002081-structural",
    "pathway-score-st002081-structural-null-sensitivity",
    "pathway-score-st000818-replication",
    "pathway-score-human-sensitivity",
    "pathway-score-human-graph-mixing",
    "pathway-score-ccle",
    "pathway-score-ccle-sensitivity",
    "pathway-score-ccle-null-ensemble",
    "pathway-score-ccle-property-matched-null",
    "pathway-score-ccle-reaction-cluster",
    "pathway-score-mapping-supplement",
    "pathway-score-publication-figures",
    "factorized-pathway-publication",
]


def _manifest(name: str) -> dict[str, object]:
    manifest, integrity = load_verified_manifest(ROOT / "artifacts" / name)
    if integrity != "verified" or manifest is None:
        raise RuntimeError(f"{name}: {integrity}")
    return manifest


def main() -> None:
    manifests = {name: _manifest(name) for name in CORE_ARTIFACTS}
    coarse = manifests["pathway-score-st002081"]
    structural = manifests["pathway-score-st002081-structural"]
    replication = manifests["pathway-score-st000818-replication"]
    ccle = manifests["pathway-score-ccle"]
    human_sensitivity = manifests["pathway-score-human-sensitivity"]
    ccle_sensitivity = manifests["pathway-score-ccle-sensitivity"]
    hard_null = manifests["pathway-score-ccle-property-matched-null"]
    synthesis = manifests["factorized-pathway-publication"]

    assert coarse["samples"] == 1539 and coarse["subjects"] == 112
    assert coarse["gate_results"]["rmse_beats_random_grouping"] is False
    assert structural["decision"] == "ADAPTIVE_STRUCTURE_RESOLUTION_SIGNAL"
    assert replication["samples"] == 450 and replication["validation_groups"] == 15
    assert replication["decision"] == "EXTERNAL_STRUCTURAL_REPLICATION"
    assert ccle["aligned_samples"] == 913 and ccle["targets"] == 60
    assert human_sensitivity["settings_per_dataset"] == 9
    assert ccle_sensitivity["settings"] == 10
    assert hard_null["random_feature_seeds"] == 20
    assert synthesis["decision"] == "ROBUST_HUMAN_STRUCTURAL_AND_REACTION_TOPOLOGY_SIGNAL"

    tracked = subprocess.run(
        ["git", "-C", str(ROOT), "ls-files", "data/raw"],
        check=False,
        capture_output=True,
        text=True,
    )
    tracked_raw_files = [line for line in tracked.stdout.splitlines() if line]
    if tracked_raw_files:
        raise RuntimeError(f"Raw source data are tracked: {tracked_raw_files}")

    summary = {
        "artifact_manifests_verified": len(manifests),
        "human_samples": coarse["samples"],
        "external_replication_samples": replication["samples"],
        "ccle_lines": ccle["aligned_samples"],
        "ccle_targets": ccle["targets"],
        "release_decision": synthesis["decision"],
        "raw_source_files_distributed": len(tracked_raw_files),
    }
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
