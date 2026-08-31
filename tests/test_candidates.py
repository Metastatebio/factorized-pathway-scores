import pandas as pd
import pytest

from genotype_gated_metabolism.candidates import generate_candidate_catalog
from genotype_gated_metabolism.model import MetabolicNetwork, Reaction, parse_gpr_to_dnf


def test_parse_gpr_to_dnf_preserves_and_or_logic() -> None:
    alternatives = parse_gpr_to_dnf("(G1 and G2) or GALT")
    assert set(alternatives) == {frozenset({"G1", "G2"}), frozenset({"GALT"})}


def test_generate_candidates_uses_reaction_path_and_enzyme_alternatives() -> None:
    network = MetabolicNetwork(
        source="fixture",
        version="1",
        reactions=(
            Reaction(
                "R1",
                "A to X",
                frozenset({"A", "X"}),
                parse_gpr_to_dnf("(G1 and G2) or GALT"),
                "fixture pathway",
            ),
            Reaction(
                "R2",
                "X to B",
                frozenset({"X", "B"}),
                parse_gpr_to_dnf("G3"),
                "fixture pathway",
            ),
        ),
        metabolite_names={"A": "A", "X": "X", "B": "B"},
    )
    metabolite_map = pd.DataFrame(
        [
            {"assay_metabolite_id": "assay_a", "human_gem_base_id": "A", "display_name": "A"},
            {"assay_metabolite_id": "assay_b", "human_gem_base_id": "B", "display_name": "B"},
        ]
    )
    frequencies = {"G1": 0.2, "G2": 0.3, "G3": 0.4, "GALT": 0.1}

    catalog = generate_candidate_catalog(
        network,
        metabolite_map,
        frequencies,
        maximum_reaction_distance=1,
        maximum_paths_per_pair=10,
        maximum_currency_metabolite_degree=10,
        maximum_genes_per_signature=3,
        maximum_candidates=100,
        total_cohort_size=10_000,
        discovery_fraction=0.7,
        target_effect=0.5,
        residual_sd=1.0,
        discovery_alpha=0.01,
        replication_alpha=0.05,
        minimum_discovery_carriers=100,
        minimum_replication_carriers=50,
    )

    assert set(catalog["genes"]) == {"G1;G2;G3", "G3;GALT"}
    three_gene = catalog.loc[catalog["genes"] == "G1;G2;G3"].iloc[0]
    assert three_gene["reaction_ids"] == "R1;R2"
    assert three_gene["intermediate_metabolites"] == "X"
    assert three_gene["review_priority"] == "3_one_hop_same_subsystem"
    assert three_gene["signature_prevalence"] == pytest.approx(0.024)
    assert three_gene["expected_discovery_carriers"] == pytest.approx(168)


def test_generate_candidates_filters_promiscuous_reactions_and_subsystems() -> None:
    network = MetabolicNetwork(
        source="fixture",
        version="1",
        reactions=(
            Reaction(
                "specific",
                "A to B",
                frozenset({"A", "B"}),
                parse_gpr_to_dnf("GOOD"),
                "specific metabolism",
            ),
            Reaction(
                "wide",
                "bulk release",
                frozenset({"A", "B", "C", "D", "E"}),
                parse_gpr_to_dnf("WIDE"),
                "specific metabolism",
            ),
            Reaction(
                "degradation",
                "protein release",
                frozenset({"A", "B"}),
                parse_gpr_to_dnf("DEGRADE"),
                "Protein degradation",
            ),
        ),
        metabolite_names={item: item for item in "ABCDE"},
    )
    metabolite_map = pd.DataFrame(
        [
            {"assay_metabolite_id": "a", "human_gem_base_id": "A", "display_name": "A"},
            {"assay_metabolite_id": "b", "human_gem_base_id": "B", "display_name": "B"},
        ]
    )
    catalog = generate_candidate_catalog(
        network,
        metabolite_map,
        {},
        maximum_reaction_distance=0,
        maximum_paths_per_pair=10,
        maximum_currency_metabolite_degree=10,
        maximum_genes_per_signature=3,
        maximum_candidates=100,
        total_cohort_size=100,
        discovery_fraction=0.7,
        target_effect=0.5,
        residual_sd=1.0,
        discovery_alpha=0.01,
        replication_alpha=0.05,
        minimum_discovery_carriers=10,
        minimum_replication_carriers=5,
        maximum_metabolites_per_reaction=4,
        excluded_subsystem_patterns=["protein degradation"],
    )

    assert catalog["genes"].tolist() == ["GOOD"]
