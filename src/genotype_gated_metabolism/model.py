"""Load a pinned Human-GEM release into the candidate-generation data model."""

from __future__ import annotations

import ast
import re
from collections import defaultdict
from dataclasses import dataclass
from itertools import product
from pathlib import Path

import pandas as pd
import yaml


@dataclass(frozen=True)
class Reaction:
    reaction_id: str
    name: str
    metabolites: frozenset[str]
    gene_sets: tuple[frozenset[str], ...]
    subsystem: str
    stoichiometry: tuple[tuple[str, float], ...] = ()
    lower_bound: float = -1000.0
    upper_bound: float = 1000.0


@dataclass(frozen=True)
class MetabolicNetwork:
    source: str
    version: str
    reactions: tuple[Reaction, ...]
    metabolite_names: dict[str, str]


def parse_gpr_to_dnf(
    expression: str, maximum_alternatives: int = 512
) -> tuple[frozenset[str], ...]:
    """Parse a Boolean GPR expression into alternative required-gene sets."""
    expression = expression.strip()
    if not expression:
        return (frozenset(),)
    normalized = re.sub(r"\bAND\b", "and", expression, flags=re.IGNORECASE)
    normalized = re.sub(r"\bOR\b", "or", normalized, flags=re.IGNORECASE)
    tree = ast.parse(normalized, mode="eval")

    def expand(node: ast.AST) -> list[frozenset[str]]:
        if isinstance(node, ast.Name):
            return [frozenset({node.id})]
        if isinstance(node, ast.BoolOp) and isinstance(node.op, ast.Or):
            alternatives = [item for child in node.values for item in expand(child)]
        elif isinstance(node, ast.BoolOp) and isinstance(node.op, ast.And):
            child_alternatives = [expand(child) for child in node.values]
            alternatives = [frozenset().union(*items) for items in product(*child_alternatives)]
        else:
            raise TypeError(f"Unsupported GPR expression element: {ast.dump(node)}")
        if len(alternatives) > maximum_alternatives:
            raise ValueError(
                f"GPR expands to more than {maximum_alternatives} enzyme alternatives."
            )
        return alternatives

    alternatives = expand(tree.body)
    return tuple(sorted(set(alternatives), key=lambda item: (len(item), tuple(sorted(item)))))


def load_human_gem(
    model_path: Path,
    genes_path: Path,
    metabolites_path: Path,
    maximum_gpr_alternatives: int = 512,
) -> MetabolicNetwork:
    """Load Human-GEM YAML and annotation TSVs with compartment-collapsed metabolites."""
    for path in (model_path, genes_path, metabolites_path):
        if not path.exists():
            raise FileNotFoundError(f"Missing Human-GEM input: {path}")

    genes = pd.read_csv(genes_path, sep="\t", dtype=str).fillna("")
    gene_symbols = dict(zip(genes["genes"], genes["geneSymbols"], strict=True))

    metabolite_annotations = pd.read_csv(metabolites_path, sep="\t", dtype=str).fillna("")
    metabolite_to_base = dict(
        zip(
            metabolite_annotations["mets"],
            metabolite_annotations["metsNoComp"],
            strict=True,
        )
    )

    sections = dict(yaml.safe_load(model_path.read_bytes()))
    metadata = sections["metaData"]
    raw_metabolites = {dict(item)["id"]: dict(item) for item in sections["metabolites"]}
    metabolite_names: dict[str, str] = {}
    for metabolite_id, base_id in metabolite_to_base.items():
        raw = raw_metabolites.get(metabolite_id, {})
        name = str(raw.get("name", "")).strip()
        if name and base_id not in metabolite_names:
            metabolite_names[base_id] = name

    reactions: list[Reaction] = []
    for raw_reaction in sections["reactions"]:
        row = dict(raw_reaction)
        compartment_metabolites = dict(row["metabolites"])
        base_metabolites = frozenset(
            metabolite_to_base[metabolite_id]
            for metabolite_id in compartment_metabolites
            if metabolite_id in metabolite_to_base
        )
        gene_sets_ensembl = parse_gpr_to_dnf(
            str(row.get("gene_reaction_rule", "")), maximum_gpr_alternatives
        )
        gene_sets_symbols = tuple(
            frozenset((gene_symbols.get(gene_id) or gene_id) for gene_id in gene_set)
            for gene_set in gene_sets_ensembl
        )
        collapsed_stoichiometry: dict[str, float] = defaultdict(float)
        for metabolite_id, coefficient in compartment_metabolites.items():
            if metabolite_id in metabolite_to_base:
                collapsed_stoichiometry[metabolite_to_base[metabolite_id]] += float(coefficient)
        stoichiometry = tuple(
            sorted(
                (metabolite_id, coefficient)
                for metabolite_id, coefficient in collapsed_stoichiometry.items()
                if not abs(coefficient) < 1e-12
            )
        )
        reactions.append(
            Reaction(
                reaction_id=str(row["id"]),
                name=str(row.get("name", "")),
                metabolites=base_metabolites,
                gene_sets=gene_sets_symbols,
                subsystem=str(row.get("subsystem", "")),
                stoichiometry=stoichiometry,
                lower_bound=float(row.get("lower_bound", -1000.0)),
                upper_bound=float(row.get("upper_bound", 1000.0)),
            )
        )

    return MetabolicNetwork(
        source="Human-GEM",
        version=str(metadata["version"]),
        reactions=tuple(reactions),
        metabolite_names=metabolite_names,
    )
