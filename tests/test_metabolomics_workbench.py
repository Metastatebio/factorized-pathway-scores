from __future__ import annotations

import pandas as pd
import pytest

from genotype_gated_metabolism.datasets.metabolomics_workbench import (
    MetabolomicsWorkbenchClient,
    parse_measurement_blocks,
    parse_record_blocks,
    parse_study_page_fields,
)


def test_parse_record_blocks_reads_blank_line_delimited_records() -> None:
    payload = "study_id\tST1\ntitle\tFirst\n\nstudy_id\tST2\ntitle\tSecond: detail\n"

    records = parse_record_blocks(payload)

    assert records == [
        {"study_id": "ST1", "title": "First"},
        {"study_id": "ST2", "title": "Second: detail"},
    ]


def test_parse_record_blocks_rejects_malformed_lines() -> None:
    with pytest.raises(ValueError, match="Malformed"):
        parse_record_blocks("study_id ST1")


def test_parse_study_page_fields_extracts_trial_design() -> None:
    payload = """
    <table><tr><td><b>Study Type</b></td><td>Randomized Controlled Trial</td></tr>
    <tr><td><b>Project DOI:</b></td><td>doi: 10.1/example</td></tr></table>
    """

    observed = parse_study_page_fields(payload)

    assert observed["study_type"] == "Randomized Controlled Trial"
    assert observed["project_doi"] == "doi: 10.1/example"


class _CatalogClient(MetabolomicsWorkbenchClient):
    def _records(self, path: str) -> list[dict[str, str]]:
        if path.endswith("available/json"):
            return [
                {"project_id": "PR1", "study_id": "ST000001", "analysis_id": "AN1"},
                {"project_id": "PR1", "study_id": "ST000001", "analysis_id": "AN2"},
                {"project_id": "PR2", "study_id": "ST001001", "analysis_id": "AN3"},
            ]
        if "ST000/summary" in path:
            return [
                {
                    "study_id": "ST000001",
                    "study_title": "Trial",
                    "species": "Homo sapiens",
                }
            ]
        if "ST001/summary" in path:
            return [
                {
                    "study_id": "ST001001",
                    "study_title": "Other",
                    "species": "Mus musculus",
                }
            ]
        raise AssertionError(path)


def test_catalog_joins_analysis_availability() -> None:
    observed = _CatalogClient().catalog(maximum_workers=2)

    assert list(observed["study_id"]) == ["ST000001", "ST001001"]
    assert observed.loc[0, "analyses_available"] == 2
    assert observed.loc[1, "project_id"] == "PR2"
    assert isinstance(observed, pd.DataFrame)


def test_parse_measurements_excludes_ambiguous_duplicate_sample_identifiers() -> None:
    payload = (
        "study_id\tST1\nmetabolite_name\tLeu\nunits\tmM\n"
        "sample_a\t0.12\nsample_b\t0.2\nsample_b\t0.3\n\n"
    )

    observed = parse_measurement_blocks(payload)

    assert observed == [
        {
            "study_id": "ST1",
            "metabolite_name": "Leu",
            "units": "mM",
            "local_sample_id": "sample_a",
            "value": "0.12",
        }
    ]
