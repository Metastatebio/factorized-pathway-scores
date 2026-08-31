"""Fetch and verify the locked ST002081/AN003790 mwTab input."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path
from urllib import request

URL = "https://www.metabolomicsworkbench.org/rest/study/analysis_id/AN003790/mwtab/txt"
EXPECTED_SHA256 = "e5d68bea4f6adf113f9bafa0e8bd333b949f289f0dc645a18435d542c3ae9c9b"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--destination",
        type=Path,
        default=Path("data/raw/st002081/ST002081_AN003790.txt"),
    )
    args = parser.parse_args()
    query = request.Request(URL, headers={"User-Agent": "factorized-pathway-scores/1.0"})
    with request.urlopen(query, timeout=180) as response:
        # Workbench currently appends two newline bytes that were not present in the pinned copy.
        # Normalizing trailing line endings preserves the exact analysis input across retrievals.
        payload = response.read().rstrip(b"\r\n")
    observed = hashlib.sha256(payload).hexdigest()
    if observed != EXPECTED_SHA256:
        raise RuntimeError(
            f"ST002081 checksum mismatch: expected {EXPECTED_SHA256}, observed {observed}"
        )
    destination = args.destination.resolve()
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_bytes(payload)
    print(f"verified {destination} sha256={observed}")


if __name__ == "__main__":
    main()

