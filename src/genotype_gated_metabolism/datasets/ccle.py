"""Reproducible adapter for the public 2019 CCLE multi-omics resource."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
from dataclasses import dataclass
from pathlib import Path
from urllib import request

import numpy as np
import pandas as pd

from .base import CrossOmicsDataset

CCLE_RELEASE = "CCLE-2019"
CCLE_BASE_URL = "https://data.broadinstitute.org/ccle"
CCLE_FILES = {
    "CCLE_metabolomics_20190502.csv": (
        "7c1d24aa575f4c58a29019026b5df8e6d1142a56925aba32ff3f1d1d5a7fd0ac"
    ),
    "CCLE_RNAseq_genes_rpkm_20180929.gct.gz": (
        "954c94233cac97a695567497dcd103bd01d342b507385dc7949777267cf2440a"
    ),
    "Cell_lines_annotations_20181226.txt": (
        "77648d1cada2f325ba2c049a5fb5408434f73b963c8e1148fb82ba2653e469e1"
    ),
}


def canonicalize_growth_medium(value: object) -> str:
    """Map heterogeneous CCLE culture-medium strings to auditable broad classes."""
    if pd.isna(value) or not str(value).strip():
        return "unknown"
    compact = re.sub(r"[^A-Z0-9]+", "", str(value).upper())
    if "DMEMF12" in compact or "DMEMHAMSF12" in compact:
        return "DMEM_F12"
    if "RPMI" in compact:
        return "RPMI"
    if "IMDM" in compact:
        return "IMDM"
    if "DMEM" in compact:
        return "DMEM"
    if "EMEM" in compact:
        return "EMEM"
    if "MCCOY" in compact:
        return "McCoy"
    if "LEIBOVITZL15" in compact or compact.startswith("L15"):
        return "L15"
    if "HAMSF12" in compact or compact.startswith("F12"):
        return "F12"
    if "MEM" in compact:
        return "MEM"
    return "other"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fetch_ccle(destination: Path, force: bool = False) -> dict[str, object]:
    """Download the pinned CCLE files and verify their content hashes."""
    destination.mkdir(parents=True, exist_ok=True)
    observed: dict[str, str] = {}
    for filename, expected_checksum in CCLE_FILES.items():
        target = destination / filename
        if force or not target.exists():
            with request.urlopen(f"{CCLE_BASE_URL}/{filename}") as response:
                payload = response.read()
            target.write_bytes(payload)
        checksum = _sha256(target)
        if checksum != expected_checksum:
            raise RuntimeError(
                f"Checksum mismatch for {filename}: expected {expected_checksum}, "
                f"observed {checksum}"
            )
        observed[filename] = checksum

    manifest: dict[str, object] = {
        "dataset": "Cancer Cell Line Encyclopedia",
        "release": CCLE_RELEASE,
        "source": CCLE_BASE_URL,
        "files": observed,
    }
    (destination / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


def _load_expression(
    path: Path, genes: list[str], strict: bool = True
) -> tuple[pd.DataFrame, list[str]]:
    requested = set(genes)
    if not requested:
        raise ValueError("At least one expression gene must be requested.")

    with gzip.open(path, "rt") as handle:
        version = handle.readline().strip()
        dimensions = handle.readline().strip()
        header = handle.readline().rstrip("\n").split("\t")
        if version != "#1.2" or len(header) < 3:
            raise ValueError("Unexpected CCLE GCT expression format.")
        samples = header[2:]
        rows: dict[str, list[np.ndarray]] = {}
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) != len(header):
                continue
            symbol = fields[1]
            if symbol in requested:
                rows.setdefault(symbol, []).append(np.asarray(fields[2:], dtype=float))

    missing = sorted(requested - set(rows))
    if missing and strict:
        raise KeyError(f"Genes absent from CCLE expression data: {', '.join(missing)}")
    values = {
        gene: np.log2(1.0 + np.mean(np.vstack(gene_rows), axis=0))
        for gene, gene_rows in rows.items()
    }
    expression = pd.DataFrame(values, index=pd.Index(samples, name="sample_id"))
    expression.attrs["gct_dimensions"] = dimensions
    expression.attrs["transform"] = "log2(RPKM + 1)"
    return expression.sort_index(axis=1), missing


def load_ccle(raw_dir: Path, genes: list[str]) -> CrossOmicsDataset:
    """Load selected CCLE expression features with metabolomics and annotations."""
    required_paths = {filename: raw_dir / filename for filename in CCLE_FILES}
    missing_files = [name for name, path in required_paths.items() if not path.exists()]
    if missing_files:
        raise FileNotFoundError(
            "Missing CCLE files; run genotype-metabolism-fetch-ccle first: "
            + ", ".join(missing_files)
        )
    for filename, expected in CCLE_FILES.items():
        observed = _sha256(required_paths[filename])
        if observed != expected:
            raise RuntimeError(f"CCLE source checksum changed for {filename}.")

    metabolomics_source = pd.read_csv(required_paths["CCLE_metabolomics_20190502.csv"])
    depmap_ids = (
        metabolomics_source.loc[:, ["CCLE_ID", "DepMap_ID"]]
        .rename(columns={"CCLE_ID": "sample_id", "DepMap_ID": "metabolomics_depmap_id"})
        .set_index("sample_id")
    )
    metabolomics = metabolomics_source.rename(columns={"CCLE_ID": "sample_id"}).set_index(
        "sample_id"
    )
    metabolomics = metabolomics.drop(columns=["DepMap_ID"], errors="ignore")
    metabolomics = metabolomics.apply(pd.to_numeric, errors="coerce")

    annotations = pd.read_csv(required_paths["Cell_lines_annotations_20181226.txt"], sep="\t")
    annotations = annotations.rename(
        columns={"CCLE_ID": "sample_id", "depMapID": "depmap_id"}
    ).set_index("sample_id")
    annotations = annotations.join(depmap_ids, how="left")
    annotations["depmap_id"] = annotations["depmap_id"].fillna(
        annotations["metabolomics_depmap_id"]
    )
    annotations["lineage"] = annotations["Site_Primary"].fillna("unknown").astype(str)
    annotations["growth_medium_raw"] = annotations["Growth.Medium"]
    annotations["growth_medium"] = annotations["Growth.Medium"].map(canonicalize_growth_medium)
    annotations["growth_medium_rpmi"] = np.where(
        annotations["growth_medium"].eq("RPMI"),
        1.0,
        np.where(annotations["growth_medium"].eq("unknown"), np.nan, 0.0),
    )

    expression, missing_genes = _load_expression(
        required_paths["CCLE_RNAseq_genes_rpkm_20180929.gct.gz"], genes, strict=False
    )
    dataset = CrossOmicsDataset(
        sample_metadata=annotations,
        blocks={"metabolomics": metabolomics, "transcriptomics": expression},
        provenance={
            "dataset": "CCLE",
            "release": CCLE_RELEASE,
            "metabolomics_transform": "source-provided cleaned log10 abundance",
            "transcriptomics_transform": expression.attrs["transform"],
            "requested_expression_genes": len(genes),
            "missing_expression_genes": missing_genes,
            "model_system": "cancer cell lines",
        },
    )
    return dataset.aligned(["metabolomics", "transcriptomics"])


@dataclass(frozen=True)
class CCLEPaths:
    """Resolved source paths retained for future adapter extensions."""

    raw_dir: Path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--destination", type=Path, default=Path("data/raw/ccle/2019"))
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    manifest = fetch_ccle(args.destination.resolve(), force=args.force)
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
