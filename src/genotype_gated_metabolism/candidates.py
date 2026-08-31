"""Generate pathway-constrained gene-signature/metabolite-pair hypotheses."""

from __future__ import annotations

from collections import defaultdict, deque
from collections.abc import Sequence
from itertools import combinations, pairwise, product
from typing import Any

import numpy as np
import pandas as pd

from .model import MetabolicNetwork, Reaction
from .power import interaction_power_normal, joint_signature_prevalence


def build_reaction_adjacency(
    reactions: tuple[Reaction, ...], maximum_metabolite_degree: int
) -> tuple[dict[str, set[str]], set[str], dict[str, set[str]]]:
    """Connect reactions sharing non-currency metabolites."""
    metabolite_reactions: dict[str, set[str]] = defaultdict(set)
    for reaction in reactions:
        for metabolite in reaction.metabolites:
            metabolite_reactions[metabolite].add(reaction.reaction_id)

    currency_metabolites = {
        metabolite
        for metabolite, reaction_ids in metabolite_reactions.items()
        if len(reaction_ids) > maximum_metabolite_degree
    }
    adjacency: dict[str, set[str]] = {reaction.reaction_id: set() for reaction in reactions}
    for metabolite, reaction_ids in metabolite_reactions.items():
        if metabolite in currency_metabolites:
            continue
        for left, right in combinations(sorted(reaction_ids), 2):
            adjacency[left].add(right)
            adjacency[right].add(left)
    return adjacency, currency_metabolites, metabolite_reactions


def bounded_reaction_paths(
    starts: set[str],
    goals: set[str],
    adjacency: dict[str, set[str]],
    maximum_reaction_distance: int,
    maximum_paths: int,
) -> tuple[list[tuple[str, ...]], bool]:
    """Return deterministic simple paths containing at most distance+1 reactions."""
    queue = deque((start,) for start in sorted(starts))
    paths: list[tuple[str, ...]] = []
    truncated = False
    while queue:
        path = queue.popleft()
        current = path[-1]
        if current in goals:
            paths.append(path)
            if len(paths) >= maximum_paths:
                truncated = bool(queue)
                break
            continue
        if len(path) - 1 >= maximum_reaction_distance:
            continue
        for neighbor in sorted(adjacency.get(current, set())):
            if neighbor not in path:
                queue.append((*path, neighbor))
    return paths, truncated


def _signature_metrics(
    genes: frozenset[str],
    gene_frequencies: dict[str, float],
    total_cohort_size: int,
    discovery_fraction: float,
    target_effect: float,
    residual_sd: float,
    discovery_alpha: float,
    replication_alpha: float,
    minimum_discovery_carriers: int,
    minimum_replication_carriers: int,
) -> dict[str, Any]:
    missing = sorted(genes - set(gene_frequencies))
    if missing:
        return {
            "signature_prevalence": np.nan,
            "expected_discovery_carriers": np.nan,
            "expected_replication_carriers": np.nan,
            "discovery_power": np.nan,
            "replication_power": np.nan,
            "eligible": False,
            "ineligibility_reason": f"missing frequencies: {';'.join(missing)}",
        }

    prevalence = joint_signature_prevalence([gene_frequencies[gene] for gene in genes])
    discovery_size = int(np.floor(total_cohort_size * discovery_fraction))
    replication_size = total_cohort_size - discovery_size
    discovery_carriers = discovery_size * prevalence
    replication_carriers = replication_size * prevalence
    discovery_power = interaction_power_normal(
        discovery_size, prevalence, target_effect, residual_sd, discovery_alpha
    )
    replication_power = interaction_power_normal(
        replication_size, prevalence, target_effect, residual_sd, replication_alpha
    )
    failures: list[str] = []
    if discovery_carriers < minimum_discovery_carriers:
        failures.append("discovery carrier gate")
    if replication_carriers < minimum_replication_carriers:
        failures.append("replication carrier gate")
    if discovery_power < 0.80:
        failures.append("discovery power gate")
    if replication_power < 0.80:
        failures.append("replication power gate")
    return {
        "signature_prevalence": prevalence,
        "expected_discovery_carriers": discovery_carriers,
        "expected_replication_carriers": replication_carriers,
        "discovery_power": discovery_power,
        "replication_power": replication_power,
        "eligible": not failures,
        "ineligibility_reason": "; ".join(failures),
    }


def generate_candidate_catalog(
    network: MetabolicNetwork,
    measured_metabolites: pd.DataFrame,
    gene_frequencies: dict[str, float],
    *,
    maximum_reaction_distance: int,
    maximum_paths_per_pair: int,
    maximum_currency_metabolite_degree: int,
    maximum_genes_per_signature: int,
    maximum_candidates: int,
    total_cohort_size: int,
    discovery_fraction: float,
    target_effect: float,
    residual_sd: float,
    discovery_alpha: float,
    replication_alpha: float,
    minimum_discovery_carriers: int,
    minimum_replication_carriers: int,
    maximum_metabolites_per_reaction: int | None = None,
    excluded_subsystem_patterns: Sequence[str] = (),
) -> pd.DataFrame:
    """Generate topology-supported hypotheses from explicitly mapped assay metabolites."""
    required_columns = {"assay_metabolite_id", "human_gem_base_id", "display_name"}
    missing_columns = required_columns - set(measured_metabolites.columns)
    if missing_columns:
        raise ValueError(f"Metabolite map is missing: {', '.join(sorted(missing_columns))}")
    if measured_metabolites["assay_metabolite_id"].duplicated().any():
        raise ValueError("assay_metabolite_id values must be unique.")

    excluded_patterns = tuple(pattern.lower() for pattern in excluded_subsystem_patterns)
    eligible_reactions = tuple(
        reaction
        for reaction in network.reactions
        if (
            maximum_metabolites_per_reaction is None
            or len(reaction.metabolites) <= maximum_metabolites_per_reaction
        )
        and not any(pattern in reaction.subsystem.lower() for pattern in excluded_patterns)
    )
    reaction_by_id = {reaction.reaction_id: reaction for reaction in eligible_reactions}
    adjacency, currency_metabolites, metabolite_reactions = build_reaction_adjacency(
        eligible_reactions, maximum_currency_metabolite_degree
    )
    mapped_rows = measured_metabolites.to_dict(orient="records")
    network_metabolites = set().union(*(reaction.metabolites for reaction in network.reactions))
    unknown = sorted(
        {
            str(row["human_gem_base_id"])
            for row in mapped_rows
            if str(row["human_gem_base_id"]) not in network_metabolites
        }
    )
    if unknown:
        raise ValueError(f"Mapped metabolites absent from network: {', '.join(unknown)}")

    records_by_key: dict[tuple[str, str, tuple[str, ...]], dict[str, Any]] = {}
    for left, right in combinations(mapped_rows, 2):
        left_base = str(left["human_gem_base_id"])
        right_base = str(right["human_gem_base_id"])
        if left_base == right_base:
            continue
        paths, truncated = bounded_reaction_paths(
            metabolite_reactions[left_base],
            metabolite_reactions[right_base],
            adjacency,
            maximum_reaction_distance,
            maximum_paths_per_pair,
        )
        for path in paths:
            reactions = [reaction_by_id[reaction_id] for reaction_id in path]
            for enzyme_choices in product(*(reaction.gene_sets for reaction in reactions)):
                genes = frozenset().union(*enzyme_choices)
                if not 1 <= len(genes) <= maximum_genes_per_signature:
                    continue
                sorted_genes = tuple(sorted(genes))
                key = (
                    str(left["assay_metabolite_id"]),
                    str(right["assay_metabolite_id"]),
                    sorted_genes,
                )
                existing = records_by_key.get(key)
                if existing and int(existing["reaction_distance"]) <= len(path) - 1:
                    continue

                intermediates = set()
                for first, second in pairwise(reactions):
                    intermediates.update(first.metabolites & second.metabolites)
                intermediates -= currency_metabolites
                subsystems = sorted(
                    {reaction.subsystem for reaction in reactions if reaction.subsystem}
                )
                contains_transport = any(
                    "transport" in reaction.subsystem.lower() for reaction in reactions
                )
                if len(path) == 1 and not contains_transport:
                    review_priority = "1_direct_reaction"
                elif len(path) == 1:
                    review_priority = "2_direct_transport"
                elif len(path) == 2 and len(subsystems) == 1:
                    review_priority = "3_one_hop_same_subsystem"
                elif len(path) == 2:
                    review_priority = "4_one_hop_cross_subsystem"
                else:
                    review_priority = "5_two_hop_exploratory"
                metrics = _signature_metrics(
                    genes,
                    gene_frequencies,
                    total_cohort_size,
                    discovery_fraction,
                    target_effect,
                    residual_sd,
                    discovery_alpha,
                    replication_alpha,
                    minimum_discovery_carriers,
                    minimum_replication_carriers,
                )
                records_by_key[key] = {
                    "hypothesis_id": "",
                    "metabolite_a": left["assay_metabolite_id"],
                    "metabolite_a_name": left["display_name"],
                    "metabolite_b": right["assay_metabolite_id"],
                    "metabolite_b_name": right["display_name"],
                    "direction_status": "unresolved",
                    "genes": ";".join(sorted_genes),
                    "gene_count": len(sorted_genes),
                    "reaction_ids": ";".join(path),
                    "reaction_distance": len(path) - 1,
                    "intermediate_metabolites": ";".join(sorted(intermediates)),
                    "subsystems": ";".join(subsystems),
                    "contains_transport_reaction": contains_transport,
                    "review_priority": review_priority,
                    "pathway_source": f"{network.source} v{network.version}",
                    "path_search_truncated": truncated,
                    **metrics,
                }
                if len(records_by_key) > maximum_candidates:
                    raise RuntimeError(
                        f"Candidate count exceeded configured maximum ({maximum_candidates})."
                    )

    records = sorted(
        records_by_key.values(),
        key=lambda row: (
            str(row["review_priority"]),
            str(row["metabolite_a"]),
            str(row["metabolite_b"]),
            int(row["gene_count"]),
            str(row["genes"]),
        ),
    )
    for index, record in enumerate(records, start=1):
        record["hypothesis_id"] = f"GGM-{index:07d}"
    return pd.DataFrame.from_records(records)
