"""Read public study metadata from the Metabolomics Workbench REST service."""

from __future__ import annotations

import hashlib
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from html.parser import HTMLParser
from pathlib import Path
from urllib import error, request

import pandas as pd

DEFAULT_BASE_URL = "https://www.metabolomicsworkbench.org/rest"
DATA_METADATA_FIELDS = {
    "study_id",
    "analysis_id",
    "analysis_summary",
    "metabolite_name",
    "metabolite_id",
    "refmet_name",
    "units",
}


class _TableFieldParser(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self._in_cell = False
        self._cell_parts: list[str] = []
        self._row: list[str] = []
        self.fields: dict[str, str] = {}

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        del attrs
        if tag.lower() in {"td", "th"}:
            self._in_cell = True
            self._cell_parts = []
        elif tag.lower() == "br" and self._in_cell:
            self._cell_parts.append(" ")

    def handle_data(self, data: str) -> None:
        if self._in_cell:
            self._cell_parts.append(data)

    def handle_endtag(self, tag: str) -> None:
        lowered = tag.lower()
        if lowered in {"td", "th"} and self._in_cell:
            value = " ".join("".join(self._cell_parts).split())
            self._row.append(value)
            self._in_cell = False
            self._cell_parts = []
        elif lowered == "tr":
            if len(self._row) >= 2 and self._row[0]:
                key = self._row[0].rstrip(":").strip()
                value = " ".join(item for item in self._row[1:] if item)
                if key in self.fields and value:
                    self.fields[key] = f"{self.fields[key]}\n{value}"
                elif value:
                    self.fields[key] = value
            self._row = []


def parse_study_page_fields(payload: str) -> dict[str, str]:
    """Extract labelled repository fields from a public study metadata page."""
    parser = _TableFieldParser()
    parser.feed(payload)
    normalized: dict[str, str] = {}
    for label, value in parser.fields.items():
        key = "".join(character.lower() if character.isalnum() else "_" for character in label)
        key = "_".join(part for part in key.split("_") if part)
        normalized[key] = value
    return normalized


def parse_record_blocks(payload: str) -> list[dict[str, str]]:
    """Parse Workbench's blank-line-delimited key/value response format."""
    records: list[dict[str, str]] = []
    current: dict[str, str] = {}
    for raw_line in [*payload.splitlines(), ""]:
        line = raw_line.rstrip("\r\n")
        if not line.strip():
            if current:
                records.append(current)
                current = {}
            continue
        if "\t" not in line:
            raise ValueError(f"Malformed Metabolomics Workbench response line: {line!r}")
        key, value = line.split("\t", 1)
        if key in current:
            raise ValueError(f"Duplicate field {key!r} in one Workbench record.")
        current[key] = value
    return records


def parse_measurement_blocks(payload: str) -> list[dict[str, str]]:
    """Parse processed data while excluding ambiguous duplicate sample identifiers."""
    records: list[dict[str, str]] = []
    for block in payload.split("\n\n"):
        rows: list[tuple[str, str]] = []
        for line in block.splitlines():
            if not line.strip():
                continue
            if "\t" not in line:
                raise ValueError(f"Malformed Metabolomics Workbench response line: {line!r}")
            rows.append(tuple(line.split("\t", 1)))
        if not rows:
            continue
        metadata: dict[str, str] = {}
        sample_rows: list[tuple[str, str]] = []
        for key, value in rows:
            if key in DATA_METADATA_FIELDS:
                if key in metadata:
                    raise ValueError(f"Duplicate data metadata field {key!r}.")
                metadata[key] = value
            else:
                sample_rows.append((key, value))
        sample_counts = pd.Series([key for key, _ in sample_rows]).value_counts()
        duplicate_samples = set(sample_counts.loc[sample_counts.gt(1)].index)
        records.extend(
            {
                **metadata,
                "local_sample_id": sample_id,
                "value": value,
            }
            for sample_id, value in sample_rows
            if sample_id not in duplicate_samples
        )
    return records


class MetabolomicsWorkbenchClient:
    """Small retrying client for public Workbench study metadata."""

    def __init__(
        self,
        base_url: str = DEFAULT_BASE_URL,
        *,
        site_url: str | None = None,
        cache_dir: Path | None = None,
        refresh: bool = False,
        timeout: float = 90.0,
        attempts: int = 3,
        user_agent: str = "genotype-gated-metabolism/0.1",
    ) -> None:
        if attempts < 1:
            raise ValueError("attempts must be positive.")
        self.base_url = base_url.rstrip("/")
        self.site_url = (site_url or self.base_url.removesuffix("/rest")).rstrip("/")
        self.cache_dir = cache_dir
        self.refresh = refresh
        self.timeout = timeout
        self.attempts = attempts
        self.user_agent = user_agent

    def _text_url(self, url: str) -> str:
        cache_path = None
        if self.cache_dir is not None:
            cache_name = hashlib.sha256(url.encode()).hexdigest() + ".txt"
            cache_path = self.cache_dir / cache_name
            if cache_path.exists() and not self.refresh:
                return cache_path.read_text()
        last_error: Exception | None = None
        for attempt in range(self.attempts):
            try:
                query = request.Request(url, headers={"User-Agent": self.user_agent})
                with request.urlopen(query, timeout=self.timeout) as response:
                    payload = response.read().decode("utf-8", errors="replace")
                if cache_path is not None:
                    cache_path.parent.mkdir(parents=True, exist_ok=True)
                    cache_path.write_text(payload)
                return payload
            except (TimeoutError, error.URLError, error.HTTPError) as exc:
                last_error = exc
                if attempt + 1 < self.attempts:
                    time.sleep(0.5 * (2**attempt))
        raise RuntimeError(f"Metabolomics Workbench request failed: {url}") from last_error

    def _records(self, path: str) -> list[dict[str, str]]:
        url = f"{self.base_url}/{path.lstrip('/')}"
        return parse_record_blocks(self._text_url(url))

    def available(self) -> pd.DataFrame:
        """Return public project, study, and analysis identifiers."""
        return pd.DataFrame.from_records(
            self._records("study/study_id/ST/available/json"),
            columns=["project_id", "study_id", "analysis_id"],
        )

    def summaries(self, prefix: str) -> pd.DataFrame:
        """Return study summaries matching an accession prefix."""
        if not prefix.startswith("ST"):
            raise ValueError("A Workbench study prefix must begin with ST.")
        return pd.DataFrame.from_records(
            self._records(f"study/study_id/{prefix}/summary/json")
        )

    def catalog(
        self,
        *,
        prefixes: tuple[str, ...] | None = None,
        maximum_workers: int = 4,
    ) -> pd.DataFrame:
        """Fetch one de-duplicated public study catalog."""
        available = self.available()
        if available.empty:
            raise ValueError("Metabolomics Workbench returned an empty public catalog.")
        available = available.dropna(subset=["study_id"]).copy()
        available["study_id"] = available["study_id"].astype(str)
        if prefixes is None:
            prefixes = tuple(sorted({study_id[:5] for study_id in available["study_id"]}))
        if not prefixes:
            raise ValueError("No Workbench summary prefixes were selected.")

        summaries: list[pd.DataFrame] = []
        with ThreadPoolExecutor(max_workers=maximum_workers) as executor:
            futures = {executor.submit(self.summaries, prefix): prefix for prefix in prefixes}
            for future in as_completed(futures):
                block = future.result()
                if not block.empty:
                    summaries.append(block)
        if not summaries:
            raise ValueError("Workbench summary endpoints returned no studies.")
        catalog = pd.concat(summaries, ignore_index=True)
        if "study_id" not in catalog:
            raise ValueError("Workbench summaries do not contain study_id.")
        catalog = catalog.drop_duplicates("study_id", keep="first")

        availability = (
            available.groupby("study_id", as_index=False)
            .agg(
                project_id=("project_id", "first"),
                analyses_available=("analysis_id", "nunique"),
            )
            .reset_index(drop=True)
        )
        catalog = catalog.merge(availability, on="study_id", how="inner", validate="one_to_one")
        return catalog.sort_values("study_id", kind="stable").reset_index(drop=True)

    def factors(self, study_id: str) -> pd.DataFrame:
        """Return sample-level factors for one public study."""
        return pd.DataFrame.from_records(
            self._records(f"study/study_id/{study_id}/factors/json")
        )

    def metabolites(self, study_id: str) -> pd.DataFrame:
        """Return named metabolites for one public study."""
        return pd.DataFrame.from_records(
            self._records(f"study/study_id/{study_id}/metabolites/json")
        )

    def measurements(self, study_id: str) -> pd.DataFrame:
        """Return repository-processed metabolite measurements in long form."""
        url = f"{self.base_url}/study/study_id/{study_id}/data/json"
        long = pd.DataFrame.from_records(parse_measurement_blocks(self._text_url(url)))
        if long.empty:
            return pd.DataFrame(
                columns=[
                    "study_id",
                    "analysis_id",
                    "metabolite_name",
                    "metabolite_id",
                    "refmet_name",
                    "units",
                    "local_sample_id",
                    "value",
                ]
            )
        long["value"] = pd.to_numeric(long["value"], errors="coerce")
        return long.dropna(subset=["value"]).reset_index(drop=True)

    def details(self, study_id: str) -> dict[str, str]:
        """Return labelled fields omitted from the compact REST summary."""
        url = (
            f"{self.site_url}/data/DRCCMetadata.php?Mode=Study&StudyID={study_id}"
            "&StudyType=MS&ResultType=1"
        )
        details = parse_study_page_fields(self._text_url(url))
        details["study_page_url"] = url
        return details
