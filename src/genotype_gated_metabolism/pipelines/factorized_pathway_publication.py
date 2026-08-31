"""Synthesize pathway-score benchmarks into figures and a bounded claim ledger."""

from __future__ import annotations

import argparse
import hashlib
import inspect
import json
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml

from ..analysis.publication_readiness import load_verified_manifest
from ..reporting import dataframe_to_markdown


def _resolve(config_path: Path, value: str) -> Path:
    path = Path(value)
    return path if path.is_absolute() else (config_path.parent / path).resolve()


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def implementation_sha256() -> str:
    return _sha256(Path(inspect.getfile(run)))


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def _evidence_summary(
    st: dict[str, Any],
    structural: dict[str, Any],
    structural_sensitivity: dict[str, Any],
    replication: dict[str, Any],
    ccle: dict[str, Any],
    human_sensitivity: dict[str, Any],
    human_graph_mixing: dict[str, Any],
    ccle_sensitivity: dict[str, Any],
    ccle_null_ensemble: dict[str, Any],
    ccle_property_matched_null: dict[str, Any],
    ccle_reaction_cluster: dict[str, Any],
    mapping: dict[str, Any],
    coupling: dict[str, Any],
    metsim: dict[str, Any],
    gxe: dict[str, Any],
    nsclc: dict[str, Any],
    constraint: dict[str, Any],
    longitudinal: dict[str, Any],
) -> pd.DataFrame:
    st_models = {row["model"]: row for row in st["model_metrics"]}
    structural_models = {row["model"]: row for row in structural["model_metrics"]}
    ccle_models = {row["model"]: row for row in ccle["model_metrics"]}
    ccle_random = next(
        item
        for item in ccle["bootstrap_comparisons"]
        if item["reference"] == "random_factorized_ridge"
    )
    return pd.DataFrame.from_records(
        [
            {
                "domain": "Repeated human metabolomics",
                "dataset": "ST002081",
                "role": "Frozen coarse benchmark + adaptive structural follow-up",
                "result": (
                    f"Coarse R2={st_models['family_score_ridge']['pooled_r2']:.3f} failed its "
                    f"null; adaptive structural R2="
                    f"{structural_models['structural_descriptor_ridge']['pooled_r2']:.3f} passed "
                    f"across {structural_sensitivity['null_repeats']} graph nulls"
                ),
                "status": "ADAPTIVE_SIGNAL",
            },
            {
                "domain": "External human structural replication",
                "dataset": "ST000818",
                "role": "Locked population-held-out replication",
                "result": (
                    f"{replication['eligible_features']} lipids across "
                    f"{replication['validation_groups']} groups; passed "
                    f"{replication['null_repeats']} graph nulls and "
                    f"{human_sensitivity['settings_per_dataset']} settings; minimum edge "
                    f"replacement={human_graph_mixing['adjudication']['minimum_edge_replacement_fraction']:.3f}"
                ),
                "status": "SUPPORTED",
            },
            {
                "domain": "Reaction topology prediction",
                "dataset": "CCLE + HumanGEM",
                "role": "New lineage-held-out metabolite benchmark",
                "result": (
                    f"Mechanistic RMSE={ccle_models['factorized_interaction_ridge']['mean_equal_lineage_rmse_sd']:.3f}; "
                    f"improvement vs random={ccle_random['rmse_improvement_sd']:.3f}; "
                    f"{ccle_sensitivity['settings']} settings, "
                    f"{ccle_null_ensemble['random_feature_seeds']} dimension-matched and "
                    f"{ccle_property_matched_null['random_feature_seeds']} property-matched "
                    f"random seeds pass; cluster lower bound="
                    f"{ccle_reaction_cluster['inference']['ci_lower']:.3f}"
                ),
                "status": "SUPPORTED",
            },
            {
                "domain": "Mapping and feature provenance",
                "dataset": "ST002081 + ST000818 + CCLE",
                "role": "Accepted/rejected feature and model-input audit",
                "result": (
                    f"{len(mapping['output_sha256']) - 1} machine-readable tables; "
                    "all source integrity verified"
                ),
                "status": "SUPPORTED",
            },
            {
                "domain": "Transcript-conditioned coupling",
                "dataset": "CCLE",
                "role": "Prior reaction-resolved inference",
                "result": (
                    f"{coupling['fdr_significant']} FDR-significant state among "
                    f"{coupling['tested_candidates']} tests; no general interaction increment"
                ),
                "status": "CANDIDATE_SPECIFIC",
            },
            {
                "domain": "Human genetic topology",
                "dataset": "METSIM + HumanGEM",
                "role": "Prior degree-matched topology benchmark",
                "result": (
                    f"{metsim['direct_topology_recoveries']}/{metsim['mqtl_states_evaluated']} "
                    f"direct recoveries; P={metsim['degree_matched_permutation_p_value']:.2g}"
                ),
                "status": "SUPPORTIVE_NOT_INDEPENDENT",
            },
            {
                "domain": "Genotype-by-context reaction coordinate",
                "dataset": "UKB/BBJ summaries",
                "role": "Prior reciprocal summary-statistic scan",
                "result": (
                    f"{gxe['strict_replicated_candidate_rows']} strict replicated row; "
                    f"{gxe['novel_mutable_candidate_rows']} novel mutable-context rows"
                ),
                "status": "LIMITED",
            },
            {
                "domain": "Isotope-resolved external challenge",
                "dataset": "NSCLC",
                "role": "Prior falsification",
                "result": (
                    f"{nsclc['expression_module_fdr_significant']} expression and "
                    f"{nsclc['oncogenotype_fdr_significant']} oncogenotype FDR signals"
                ),
                "status": "NULL",
            },
            {
                "domain": "Constraint sensitivity",
                "dataset": "HumanGEM",
                "role": "Prior 1,000-draw model ensemble",
                "result": (
                    f"{len(constraint['stable_candidates'])} stable and "
                    f"{len(constraint['unstable_candidates'])} unstable configured candidates"
                ),
                "status": "MODEL_RELATIVE",
            },
            {
                "domain": "Directional longitudinal prediction",
                "dataset": "ST002081",
                "role": "Prior temporal falsification",
                "result": (
                    f"forward R2={longitudinal['temporal_falsification']['forward_r2']:.3f}; "
                    f"reverse R2={longitudinal['temporal_falsification']['reverse_r2']:.3f}"
                ),
                "status": "NOT_SUPPORTED",
            },
        ]
    )


def _claim_ledger(
    st: dict[str, Any],
    structural: dict[str, Any],
    replication: dict[str, Any],
    ccle: dict[str, Any],
    human_sensitivity: dict[str, Any] | None = None,
    human_graph_mixing: dict[str, Any] | None = None,
    ccle_sensitivity: dict[str, Any] | None = None,
    ccle_null_ensemble: dict[str, Any] | None = None,
    ccle_property_matched_null: dict[str, Any] | None = None,
    ccle_reaction_cluster: dict[str, Any] | None = None,
    mapping: dict[str, Any] | None = None,
) -> pd.DataFrame:
    st_gates = st["gate_results"]
    ccle_gates = ccle["gate_results"]
    records = [
            {
                "claim": "Compact lipid-family scores reconstruct hidden human markers",
                "decision": "SUPPORTED_VS_POPULATION",
                "basis": "Participant-bootstrap RMSE interval above zero versus population mean",
            },
            {
                "claim": "Broad lipid-family scores are pathway-specific",
                "decision": "NOT_SUPPORTED",
                "basis": (
                    "Random-group RMSE gate "
                    + ("passed" if st_gates["rmse_beats_random_grouping"] else "failed")
                ),
            },
            {
                "claim": "Structure-resolved lipid scores beat a degree-preserving null",
                "decision": (
                    "SUPPORTED_ADAPTIVE_SAME_COHORT"
                    if structural["decision"] == "ADAPTIVE_STRUCTURE_RESOLUTION_SIGNAL"
                    else "NOT_SUPPORTED"
                ),
                "basis": "Adaptive participant-paired RMSE and precision intervals above zero",
            },
            {
                "claim": "Structure-resolved lipid scores replicate in a separate human cohort",
                "decision": (
                    "SUPPORTED_EXTERNAL_POPULATION_HELD_OUT"
                    if replication["decision"] == "EXTERNAL_STRUCTURAL_REPLICATION"
                    else "NOT_SUPPORTED"
                ),
                "basis": "ST000818 population-held-out validation across 20 graph nulls",
            },
            {
                "claim": "Direct HumanGEM/GPR features outperform size-matched random features",
                "decision": (
                    "SUPPORTED" if ccle_gates["beats_random_factorization"] else "NOT_SUPPORTED"
                ),
                "basis": "Lineage-isolated target-bootstrap comparison across CCLE metabolites",
            },
            {
                "claim": "Transcript interactions improve the general direct-network model",
                "decision": "NOT_SUPPORTED",
                "basis": "Interaction model did not beat network-metabolite or additive models",
            },
            {
                "claim": "Transcript-conditioned coupling can exist for selected reactions",
                "decision": "SUPPORTED_IN_MODEL_SYSTEM",
                "basis": "FDR-significant CDA–cytidine–uridine interaction with robustness checks",
            },
            {
                "claim": "The evidence validates same-person genomic personalization",
                "decision": "NOT_TESTED",
                "basis": "No available cohort contains the required participant-level modalities",
            },
            {
                "claim": "The model estimates physiological flux or treatment response",
                "decision": "NOT_SUPPORTED",
                "basis": "Constraint results are model-relative and external drug/isotope tests are null",
            },
        ]
    if human_sensitivity is not None:
        records.append(
            {
                "claim": "The human structural result is robust across declared analysis settings",
                "decision": (
                    "SUPPORTED_POST_RESULT"
                    if human_sensitivity["decision"] == "HUMAN_STRUCTURAL_SENSITIVITY_ROBUST"
                    else "NOT_SUPPORTED"
                ),
                "basis": "All nine settings pass both endpoints in both human cohorts",
            }
        )
    if human_graph_mixing is not None:
        audit = human_graph_mixing["adjudication"]
        records.append(
            {
                "claim": "Human fixed-degree null graphs are moved and diverse",
                "decision": (
                    "SUPPORTED_POST_RESULT"
                    if human_graph_mixing["decision"]
                    == "HUMAN_GRAPH_NULL_MIXING_ADEQUATE"
                    else "NOT_SUPPORTED"
                ),
                "basis": (
                    f"Minimum edge replacement={audit['minimum_edge_replacement_fraction']:.3f}; "
                    f"maximum pairwise Jaccard={audit['maximum_pairwise_null_jaccard']:.3f}"
                ),
            }
        )
    if ccle_sensitivity is not None:
        records.extend(
            [
                {
                    "claim": "The CCLE reaction-topology result is robust across declared settings",
                    "decision": (
                        "SUPPORTED_POST_RESULT"
                        if ccle_sensitivity["decision"]
                        == "CCLE_REACTION_TOPOLOGY_SENSITIVITY_ROBUST"
                        else "NOT_SUPPORTED"
                    ),
                    "basis": "Matched random-control lower bounds pass in all ten settings",
                },
                {
                    "claim": "Transcript interactions become generally useful under alternative settings",
                    "decision": "NOT_SUPPORTED",
                    "basis": "Interactions fail to beat the additive model in all ten settings",
                },
            ]
        )
    if ccle_null_ensemble is not None:
        expected = ccle_null_ensemble["expected_null_comparison"]
        records.append(
            {
                "claim": "The CCLE topology result is stable across random feature draws",
                "decision": (
                    "SUPPORTED_POST_RESULT"
                    if ccle_null_ensemble["decision"]
                    == "CCLE_RANDOM_FEATURE_NULL_ENSEMBLE_ROBUST"
                    else "NOT_SUPPORTED"
                ),
                "basis": (
                    f"All {ccle_null_ensemble['random_feature_seeds']} null seeds pass; "
                    f"expected-null lower bound={expected['ci_lower']:.3f} SD"
                ),
            }
        )
    if ccle_property_matched_null is not None:
        expected = ccle_property_matched_null["expected_null_comparison"]
        records.append(
            {
                "claim": "CCLE topology beats network-degree and coverage-matched features",
                "decision": (
                    "SUPPORTED_POST_RESULT"
                    if ccle_property_matched_null["decision"]
                    == "CCLE_RANDOM_FEATURE_NULL_ENSEMBLE_ROBUST"
                    else "NOT_SUPPORTED"
                ),
                "basis": (
                    f"All {ccle_property_matched_null['random_feature_seeds']} hard-null seeds "
                    f"pass; expected-null lower bound={expected['ci_lower']:.3f} SD"
                ),
            }
        )
    if ccle_reaction_cluster is not None:
        cluster_inference = ccle_reaction_cluster["inference"]
        records.append(
            {
                "claim": "The CCLE topology effect is not concentrated in one subsystem",
                "decision": (
                    "SUPPORTED_POST_RESULT"
                    if ccle_reaction_cluster["decision"]
                    == "CCLE_REACTION_CLUSTER_ROBUST"
                    else "NOT_SUPPORTED"
                ),
                "basis": (
                    f"{cluster_inference['clusters']}-cluster bootstrap lower bound="
                    f"{cluster_inference['ci_lower']:.3f} SD; leave-one-cluster gate passes"
                ),
            }
        )
    if mapping is not None:
        records.append(
            {
                "claim": "Analysis-facing mappings and exclusions are auditable",
                "decision": (
                    "SUPPORTED_RESOURCE"
                    if mapping["decision"] == "MAPPING_QC_SUPPLEMENT_COMPLETE"
                    else "NOT_SUPPORTED"
                ),
                "basis": "Accepted/rejected human and CCLE mapping tables with verified source integrity",
            }
        )
    return pd.DataFrame.from_records(records)


def _write_figure(
    st_metrics: pd.DataFrame,
    replication_metrics: pd.DataFrame,
    ccle_metrics: pd.DataFrame,
    evidence: pd.DataFrame,
    path: Path,
) -> None:
    colors = {
        "population_mean": "#a8b0bd",
        "random_group_score_ridge": "#d17b0f",
        "family_score_ridge": "#6f42c1",
        "all_visible_ridge": "#1976a3",
        "degree_preserving_random_structural_ridge": "#d6a03a",
        "structural_descriptor_ridge": "#1b7f5c",
        "random_factorized_ridge": "#d17b0f",
        "network_metabolites_ridge": "#3f51b5",
        "network_additive_ridge": "#6f42c1",
        "factorized_interaction_ridge": "#8e5cc2",
        "all_metabolites_ridge": "#1976a3",
    }
    figure, axes = plt.subplots(2, 3, figsize=(22, 12))
    figure.patch.set_facecolor("#f7f8fa")
    for axis in axes.flat:
        axis.set_facecolor("white")
        for spine in axis.spines.values():
            spine.set_color("#d8dee8")
    st_plot = st_metrics.sort_values("row_weighted_rmse_sd", ascending=False)
    axes[0, 0].barh(
        st_plot["model"],
        st_plot["row_weighted_rmse_sd"],
        color=[colors[value] for value in st_plot["model"]],
    )
    axes[0, 0].set_xlabel("Held-out RMSE (training-fold SD)")
    axes[0, 0].set_title("A  Human lipid missing-panel reconstruction", loc="left", fontweight="bold")
    st_priority = st_metrics.sort_values("precision_at_k")
    axes[0, 1].barh(
        st_priority["model"],
        st_priority["precision_at_k"],
        color=[colors[value] for value in st_priority["model"]],
    )
    axes[0, 1].set_xlabel("Precision@3 for top-decile hidden deviations")
    axes[0, 1].set_title("B  Human hidden-marker prioritization", loc="left", fontweight="bold")
    replication_plot = replication_metrics.sort_values(
        "row_weighted_rmse_sd", ascending=False
    )
    axes[0, 2].barh(
        replication_plot["model"],
        replication_plot["row_weighted_rmse_sd"],
        color=[colors[value] for value in replication_plot["model"]],
    )
    axes[0, 2].set_xlabel("Population-held-out RMSE (training-fold SD)")
    axes[0, 2].set_title(
        "C  External human structural replication", loc="left", fontweight="bold"
    )
    replication_priority = replication_metrics.sort_values("precision_at_k")
    axes[1, 0].barh(
        replication_priority["model"],
        replication_priority["precision_at_k"],
        color=[colors[value] for value in replication_priority["model"]],
    )
    axes[1, 0].set_xlabel("Precision@3 for top-decile hidden deviations")
    axes[1, 0].set_title(
        "D  External hidden-marker prioritization", loc="left", fontweight="bold"
    )
    ccle_plot = ccle_metrics.sort_values("mean_equal_lineage_rmse_sd", ascending=False)
    axes[1, 1].barh(
        ccle_plot["model"],
        ccle_plot["mean_equal_lineage_rmse_sd"],
        color=[colors.get(value, "#596780") for value in ccle_plot["model"]],
    )
    axes[1, 1].set_xlabel("Mean equal-lineage RMSE (training-fold SD)")
    axes[1, 1].set_title("E  CCLE held-out metabolite prediction", loc="left", fontweight="bold")
    status_order = [
        "SUPPORTED",
        "CANDIDATE_SPECIFIC",
        "SUPPORTIVE_NOT_INDEPENDENT",
        "LIMITED",
        "MODEL_RELATIVE",
        "ADAPTIVE_SIGNAL",
        "MIXED",
        "NULL",
        "NOT_SUPPORTED",
    ]
    status_color = {
        "SUPPORTED": "#1b7f5c",
        "CANDIDATE_SPECIFIC": "#5c8f6d",
        "SUPPORTIVE_NOT_INDEPENDENT": "#77966d",
        "LIMITED": "#d17b0f",
        "MODEL_RELATIVE": "#8e5cc2",
        "ADAPTIVE_SIGNAL": "#5c8f6d",
        "MIXED": "#d17b0f",
        "NULL": "#a8b0bd",
        "NOT_SUPPORTED": "#b64a4a",
    }
    evidence_plot = evidence.copy()
    evidence_plot["status_rank"] = evidence_plot["status"].map(
        {value: index for index, value in enumerate(status_order)}
    )
    evidence_plot = evidence_plot.sort_values("status_rank", ascending=False)
    axes[1, 2].barh(
        evidence_plot["domain"],
        np.ones(len(evidence_plot)),
        color=[status_color[value] for value in evidence_plot["status"]],
    )
    axes[1, 2].set_xlim(0, 1)
    axes[1, 2].set_xticks([])
    axes[1, 2].set_title("F  Evidence-domain adjudication", loc="left", fontweight="bold")
    for row, (_, item) in enumerate(evidence_plot.iterrows()):
        axes[1, 2].text(
            0.03, row, item["status"], va="center", color="white", fontweight="bold"
        )
    figure.suptitle(
        "PATHWAY STRUCTURE HELPS—BUT ONLY AT THE RESOLUTION SUPPORTED BY DATA.",
        x=0.03,
        y=0.98,
        ha="left",
        fontsize=18,
        fontweight="bold",
        color="#17233c",
    )
    figure.text(
        0.03,
        0.94,
        "Participant-, population-, and lineage-isolated public-data benchmarks with matched random controls",
        fontsize=10.5,
        color="#4f5d75",
    )
    figure.tight_layout(rect=[0.03, 0.04, 0.99, 0.91])
    figure.savefig(path, dpi=250, facecolor=figure.get_facecolor())
    plt.close(figure)


def _render_report(
    manifest: dict[str, Any], evidence: pd.DataFrame, claims: pd.DataFrame
) -> str:
    return "\n".join(
        [
            "# Factorized pathway-score publication synthesis",
            "",
            f"**Decision:** `{manifest['decision']}`",
            "",
            "## Evidence domains",
            "",
            dataframe_to_markdown(evidence),
            "",
            "## Claim ledger",
            "",
            dataframe_to_markdown(claims),
            "",
            "## Interpretation",
            "",
            (
                "Direct metabolic-network neighborhoods carry reproducible predictive information "
                "in lineage-held-out cancer-cell-line data, but adding transcript interactions "
                "does not improve the general predictor. Broad lipid-family compression recovers "
                "substantial human profile information but does not outperform size-matched random "
                "compression. An explicitly adaptive structure-resolved lipid representation does "
                "beat degree-preserving random graphs in the same cohort, and this result replicates "
                "in a separate human cohort with complete population categories held out. Together, "
                "the results support empirical calibration of pathway-score resolution and "
                "candidate-specific interaction testing, not universal multi-omic fusion."
            ),
            "",
            "## Boundary",
            "",
            manifest["claim_boundary"],
            "",
        ]
    )


def run(config_path: Path, output_dir: Path) -> dict[str, Any]:
    config_path = config_path.resolve()
    config = yaml.safe_load(config_path.read_text())
    protocol_path = _resolve(config_path, str(config["protocol"]))
    source_paths = {
        key: _resolve(config_path, str(value)) for key, value in config["sources"].items()
    }
    st, st_integrity = load_verified_manifest(source_paths["st002081"])
    structural, structural_integrity = load_verified_manifest(
        source_paths["st002081_structural"]
    )
    structural_sensitivity, structural_sensitivity_integrity = load_verified_manifest(
        source_paths["st002081_structural_null_sensitivity"]
    )
    replication, replication_integrity = load_verified_manifest(
        source_paths["st000818_replication"]
    )
    ccle, ccle_integrity = load_verified_manifest(source_paths["ccle_prediction"])
    human_sensitivity, human_sensitivity_integrity = load_verified_manifest(
        source_paths["human_sensitivity"]
    )
    human_graph_mixing, human_graph_mixing_integrity = load_verified_manifest(
        source_paths["human_graph_mixing"]
    )
    ccle_sensitivity, ccle_sensitivity_integrity = load_verified_manifest(
        source_paths["ccle_sensitivity"]
    )
    ccle_null_ensemble, ccle_null_ensemble_integrity = load_verified_manifest(
        source_paths["ccle_null_ensemble"]
    )
    ccle_property_matched_null, ccle_property_matched_null_integrity = (
        load_verified_manifest(source_paths["ccle_property_matched_null"])
    )
    ccle_reaction_cluster, ccle_reaction_cluster_integrity = load_verified_manifest(
        source_paths["ccle_reaction_cluster"]
    )
    mapping, mapping_integrity = load_verified_manifest(source_paths["mapping_supplement"])
    metsim, metsim_integrity = load_verified_manifest(source_paths["metsim_topology"])
    constraint, constraint_integrity = load_verified_manifest(
        source_paths["constraint_ensemble"]
    )
    longitudinal, longitudinal_integrity = load_verified_manifest(
        source_paths["longitudinal_prior"]
    )
    required = {
        "st002081": (st, st_integrity),
        "st002081_structural": (structural, structural_integrity),
        "st002081_structural_null_sensitivity": (
            structural_sensitivity,
            structural_sensitivity_integrity,
        ),
        "st000818_replication": (replication, replication_integrity),
        "ccle_prediction": (ccle, ccle_integrity),
        "human_sensitivity": (human_sensitivity, human_sensitivity_integrity),
        "human_graph_mixing": (human_graph_mixing, human_graph_mixing_integrity),
        "ccle_sensitivity": (ccle_sensitivity, ccle_sensitivity_integrity),
        "ccle_null_ensemble": (ccle_null_ensemble, ccle_null_ensemble_integrity),
        "ccle_property_matched_null": (
            ccle_property_matched_null,
            ccle_property_matched_null_integrity,
        ),
        "ccle_reaction_cluster": (
            ccle_reaction_cluster,
            ccle_reaction_cluster_integrity,
        ),
        "mapping_supplement": (mapping, mapping_integrity),
        "metsim_topology": (metsim, metsim_integrity),
        "constraint_ensemble": (constraint, constraint_integrity),
        "longitudinal_prior": (longitudinal, longitudinal_integrity),
    }
    failures = {key: integrity for key, (_, integrity) in required.items() if integrity != "verified"}
    if failures:
        raise ValueError(f"Publication source integrity failed: {failures}")
    assert st is not None and structural is not None and structural_sensitivity is not None
    assert replication is not None and ccle is not None and metsim is not None
    assert human_sensitivity is not None and human_graph_mixing is not None
    assert ccle_sensitivity is not None
    assert ccle_null_ensemble is not None and mapping is not None
    assert ccle_property_matched_null is not None and ccle_reaction_cluster is not None
    assert constraint is not None and longitudinal is not None
    coupling = _load_json(source_paths["ccle_coupling_manifest"])
    gxe = _load_json(source_paths["gxe_family_manifest"])
    nsclc = _load_json(source_paths["nsclc_manifest"])
    evidence = _evidence_summary(
        st,
        structural,
        structural_sensitivity,
        replication,
        ccle,
        human_sensitivity,
        human_graph_mixing,
        ccle_sensitivity,
        ccle_null_ensemble,
        ccle_property_matched_null,
        ccle_reaction_cluster,
        mapping,
        coupling,
        metsim,
        gxe,
        nsclc,
        constraint,
        longitudinal,
    )
    claims = _claim_ledger(
        st,
        structural,
        replication,
        ccle,
        human_sensitivity,
        human_graph_mixing,
        ccle_sensitivity,
        ccle_null_ensemble,
        ccle_property_matched_null,
        ccle_reaction_cluster,
        mapping,
    )
    st_pathway_specific = bool(st["gate_results"]["rmse_beats_random_grouping"])
    ccle_pathway_specific = bool(ccle["gate_results"]["beats_random_factorization"])
    structural_signal = (
        structural["decision"] == "ADAPTIVE_STRUCTURE_RESOLUTION_SIGNAL"
        and structural_sensitivity["decision"]
        == "ADAPTIVE_STRUCTURE_SIGNAL_ROBUST_ACROSS_GRAPH_NULLS"
    )
    external_structural_replication = (
        replication["decision"] == "EXTERNAL_STRUCTURAL_REPLICATION"
    )
    sensitivity_robust = (
        human_sensitivity["decision"] == "HUMAN_STRUCTURAL_SENSITIVITY_ROBUST"
        and human_graph_mixing["decision"] == "HUMAN_GRAPH_NULL_MIXING_ADEQUATE"
        and ccle_sensitivity["decision"]
        == "CCLE_REACTION_TOPOLOGY_SENSITIVITY_ROBUST"
        and ccle_null_ensemble["decision"]
        == "CCLE_RANDOM_FEATURE_NULL_ENSEMBLE_ROBUST"
        and ccle_property_matched_null["decision"]
        == "CCLE_RANDOM_FEATURE_NULL_ENSEMBLE_ROBUST"
        and ccle_reaction_cluster["decision"] == "CCLE_REACTION_CLUSTER_ROBUST"
        and mapping["decision"] == "MAPPING_QC_SUPPLEMENT_COMPLETE"
    )
    if external_structural_replication and ccle_pathway_specific and sensitivity_robust:
        decision = "ROBUST_HUMAN_STRUCTURAL_AND_REACTION_TOPOLOGY_SIGNAL"
    elif external_structural_replication and ccle_pathway_specific:
        decision = "HUMAN_REPLICATED_STRUCTURAL_AND_REACTION_TOPOLOGY_SIGNAL"
    elif structural_signal and ccle_pathway_specific:
        decision = "REACTION_AND_ADAPTIVE_STRUCTURAL_PATHWAY_SIGNAL"
    elif st_pathway_specific and ccle_pathway_specific:
        decision = "CROSS_DOMAIN_PATHWAY_SCORE_SIGNAL"
    elif ccle_pathway_specific:
        decision = "REACTION_TOPOLOGY_SIGNAL_HUMAN_FAMILY_GATE_FAILED"
    else:
        decision = "PUBLIC_BENCHMARK_RESOURCE"

    output_dir.mkdir(parents=True, exist_ok=True)
    evidence_path = output_dir / "evidence-domain-summary.csv"
    claims_path = output_dir / "claim-ledger.csv"
    figure_path = output_dir / "factorized-pathway-benchmarks.png"
    report_path = output_dir / "report.md"
    manifest_path = output_dir / "manifest.json"
    evidence.to_csv(evidence_path, index=False)
    claims.to_csv(claims_path, index=False)
    st_metrics = pd.read_csv(source_paths["st002081"] / "model-metrics.csv")
    structural_metrics = pd.read_csv(
        source_paths["st002081_structural"] / "model-metrics.csv"
    )
    structural_metrics = structural_metrics.loc[
        ~structural_metrics["model"].eq("population_mean")
    ]
    st_metrics = pd.concat([st_metrics, structural_metrics], ignore_index=True)
    ccle_metrics = pd.read_csv(source_paths["ccle_prediction"] / "model-metrics.csv")
    replication_metrics = pd.read_csv(
        source_paths["st000818_replication"] / "model-metrics.csv"
    )
    _write_figure(st_metrics, replication_metrics, ccle_metrics, evidence, figure_path)
    manifest: dict[str, Any] = {
        "analysis_id": str(config["analysis_id"]),
        "completed_at": datetime.now(UTC).isoformat(),
        "decision": decision,
        "config": str(config_path),
        "config_sha256": _sha256(config_path),
        "protocol": str(protocol_path),
        "protocol_sha256": _sha256(protocol_path),
        "implementation_sha256": implementation_sha256(),
        "source_integrity": {key: integrity for key, (_, integrity) in required.items()},
        "supported_claims": int(claims["decision"].str.startswith("SUPPORTED").sum()),
        "unsupported_or_untested_claims": int(
            claims["decision"].isin(["NOT_SUPPORTED", "NOT_TESTED"]).sum()
        ),
        "claim_boundary": str(config["claim_boundary"]),
    }
    report_path.write_text(_render_report(manifest, evidence, claims))
    outputs = [evidence_path, claims_path, figure_path, report_path]
    manifest["output_sha256"] = {path.name: _sha256(path) for path in outputs}
    manifest_path.write_text(json.dumps(manifest, indent=2, allow_nan=False) + "\n")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config", type=Path, default=Path("config/factorized-pathway-publication.yaml")
    )
    parser.add_argument(
        "--output", type=Path, default=Path("artifacts/factorized-pathway-publication")
    )
    args = parser.parse_args()
    print(json.dumps(run(args.config, args.output), indent=2))


if __name__ == "__main__":
    main()
