"""Fetch and verify the pinned Human-GEM source files."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path
from urllib import request

VERSION = "v2.0.0"
BASE_URL = f"https://raw.githubusercontent.com/SysBioChalmers/Human-GEM/{VERSION}/model"
FILES = {
    "Human-GEM.yml": "01a8ea826bdaa36511b9eed28eedcb5e73690333860856c682387c6380dd3fc6",
    "genes.tsv": "2a6058a157b3b9f3c958ba753e71d7a702d2716742972741d2ef36535ae9dff6",
    "metabolites.tsv": "9d58bb3a9cd217c47dc0d54ea9517b0a23e0603bbb60d5428551e0449af1eab3",
}
TASK_BASE_URL = (
    f"https://raw.githubusercontent.com/SysBioChalmers/Human-GEM/{VERSION}"
    "/data/metabolicTasks"
)
TASK_FILES = {
    "metabolicTasks_Essential.txt": (
        "d0c58503fc44d1eb49ea1b89f4cc90930a06bcd2c3972083c0986cfb126a0f18"
    ),
}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--destination", type=Path, default=Path("data/raw/human-gem/v2.0.0"))
    args = parser.parse_args()
    destination = args.destination.resolve()
    destination.mkdir(parents=True, exist_ok=True)

    for filename, expected_checksum in FILES.items():
        target = destination / filename
        with request.urlopen(f"{BASE_URL}/{filename}") as response:
            payload = response.read()
        observed_checksum = hashlib.sha256(payload).hexdigest()
        if observed_checksum != expected_checksum:
            raise RuntimeError(
                f"Checksum mismatch for {filename}: expected {expected_checksum}, "
                f"observed {observed_checksum}"
            )
        target.write_bytes(payload)
        print(f"verified {filename} sha256={observed_checksum}")

    task_destination = destination / "metabolicTasks"
    task_destination.mkdir(parents=True, exist_ok=True)
    for filename, expected_checksum in TASK_FILES.items():
        target = task_destination / filename
        with request.urlopen(f"{TASK_BASE_URL}/{filename}") as response:
            payload = response.read()
        observed_checksum = hashlib.sha256(payload).hexdigest()
        if observed_checksum != expected_checksum:
            raise RuntimeError(
                f"Checksum mismatch for {filename}: expected {expected_checksum}, "
                f"observed {observed_checksum}"
            )
        target.write_bytes(payload)
        print(f"verified metabolicTasks/{filename} sha256={observed_checksum}")


if __name__ == "__main__":
    main()
