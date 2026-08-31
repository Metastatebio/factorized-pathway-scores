import pandas as pd

from genotype_gated_metabolism.metabolite_mapping import (
    expand_panel_by_exact_model_names,
    map_metabolite_panel,
    normalize_metabolite_name,
)


def test_map_metabolite_panel_handles_unique_ambiguous_and_conflicting_ids() -> None:
    annotations = pd.DataFrame(
        [
            {"metsNoComp": "M1", "metHMDBID": "HMDB1", "metChEBIID": "CHEBI:1"},
            {"metsNoComp": "M2", "metHMDBID": "HMDB2", "metChEBIID": "CHEBI:2"},
            {"metsNoComp": "M3", "metHMDBID": "HMDB2", "metChEBIID": "CHEBI:3"},
        ]
    )
    panel = pd.DataFrame(
        [
            {
                "assay_metabolite_id": "unique",
                "display_name": "one",
                "hmdb_id": "HMDB1",
                "chebi_id": "CHEBI:1",
            },
            {
                "assay_metabolite_id": "ambiguous",
                "display_name": "two",
                "hmdb_id": "HMDB2",
                "chebi_id": "",
            },
            {
                "assay_metabolite_id": "conflict",
                "display_name": "three",
                "hmdb_id": "HMDB1",
                "chebi_id": "CHEBI:2",
            },
            {
                "assay_metabolite_id": "missing",
                "display_name": "four",
                "hmdb_id": "HMDB9",
                "chebi_id": "",
            },
        ]
    )

    mapped = map_metabolite_panel(panel, annotations).set_index("assay_metabolite_id")

    assert mapped.loc["unique", "mapping_status"] == "mapped"
    assert mapped.loc["unique", "human_gem_base_id"] == "M1"
    assert mapped.loc["ambiguous", "mapping_status"] == "ambiguous"
    assert mapped.loc["conflict", "mapping_status"] == "conflicting_identifiers"
    assert mapped.loc["missing", "mapping_status"] == "unmapped"


def test_normalize_metabolite_name_preserves_biochemical_qualifiers() -> None:
    assert normalize_metabolite_name("α-ketoglutarate") == "alphaketoglutarate"
    assert normalize_metabolite_name("L-serine") != normalize_metabolite_name("D-serine")


def test_expand_panel_uses_only_unique_exact_names_and_curated_rows_win() -> None:
    curated = pd.DataFrame(
        [
            {
                "assay_metabolite_id": "citrate",
                "display_name": "curated citrate",
                "hmdb_id": "HMDBC",
                "chebi_id": "",
            }
        ]
    )
    annotations = pd.DataFrame(
        [
            {"metsNoComp": "M1", "metHMDBID": "HMDBC", "metChEBIID": ""},
            {"metsNoComp": "M2", "metHMDBID": "HMDBA", "metChEBIID": "CHEBI:2"},
            {"metsNoComp": "M3", "metHMDBID": "HMDBD", "metChEBIID": "CHEBI:3"},
            {"metsNoComp": "M4", "metHMDBID": "HMDBD2", "metChEBIID": "CHEBI:4"},
        ]
    )
    expanded = expand_panel_by_exact_model_names(
        curated,
        ["citrate", "alpha-ketoglutarate", "D-serine", "ambiguous"],
        {
            "M1": "citrate",
            "M2": "alpha ketoglutarate",
            "M3": "ambiguous",
            "M4": "ambiguous",
        },
        annotations,
    ).set_index("assay_metabolite_id")

    assert set(expanded.index) == {"citrate", "alpha-ketoglutarate"}
    assert expanded.loc["citrate", "display_name"] == "curated citrate"
    assert expanded.loc["citrate", "mapping_origin"] == "curated_stable_id"
    assert expanded.loc["alpha-ketoglutarate", "hmdb_id"] == "HMDBA"
