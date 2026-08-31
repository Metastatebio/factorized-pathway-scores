"""Run a pathway-constrained public CCLE coupling proof-of-concept."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from sklearn.model_selection import train_test_split
from statsmodels.stats.multitest import multipletests

from ..analysis.coupling import (
    fit_coupling,
    permutation_p_value,
    refit_without_most_influential,
)
from ..analysis.prediction import compare_interaction_prediction
from ..analysis.robustness import (
    leave_one_group_out_summary,
    winsorize_within_groups,
)
from ..candidates import generate_candidate_catalog
from ..datasets.ccle import fetch_ccle, load_ccle
from ..features.signatures import (
    build_expression_signature,
    parse_gene_signature,
    zscore_within_groups,
)
from ..metabolite_mapping import (
    expand_panel_by_exact_model_names,
    map_metabolite_panel,
    read_human_gem_metabolites,
)
from ..model import load_human_gem


def _resolve(project_root: Path, configured_path: str) -> Path:
    path = Path(configured_path)
    return path if path.is_absolute() else project_root / path


def _file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _candidate_catalog(
    config: dict[str, object], project_root: Path
) -> tuple[pd.DataFrame, pd.DataFrame]:
    network_config = config["network"]
    dataset_config = config["dataset"]
    candidate_config = config["candidate_generation"]
    network = load_human_gem(
        _resolve(project_root, network_config["model"]),
        _resolve(project_root, network_config["genes"]),
        _resolve(project_root, network_config["metabolites"]),
        int(network_config["maximum_gpr_alternatives"]),
    )
    panel = pd.read_csv(_resolve(project_root, dataset_config["assay_panel"]), dtype=str)
    annotations = read_human_gem_metabolites(_resolve(project_root, network_config["metabolites"]))
    if bool(dataset_config.get("expand_exact_name_matches", False)):
        raw_metabolomics = _resolve(project_root, dataset_config["raw_dir"]) / (
            "CCLE_metabolomics_20190502.csv"
        )
        assay_columns = pd.read_csv(raw_metabolomics, nrows=0).columns
        panel = expand_panel_by_exact_model_names(
            panel,
            [column for column in assay_columns if column not in {"CCLE_ID", "DepMap_ID"}],
            network.metabolite_names,
            annotations,
        )
    mapping = map_metabolite_panel(panel, annotations)
    mapped = mapping.loc[mapping["mapping_status"].eq("mapped")].copy()
    if mapped.empty:
        raise ValueError("No CCLE assay metabolites mapped uniquely to Human-GEM.")

    catalog = generate_candidate_catalog(
        network,
        mapped,
        gene_frequencies={},
        maximum_reaction_distance=int(network_config["maximum_reaction_distance"]),
        maximum_paths_per_pair=int(network_config["maximum_paths_per_metabolite_pair"]),
        maximum_currency_metabolite_degree=int(
            network_config["maximum_currency_metabolite_degree"]
        ),
        maximum_genes_per_signature=int(candidate_config["maximum_genes_per_signature"]),
        maximum_candidates=int(network_config["maximum_candidates"]),
        total_cohort_size=928,
        discovery_fraction=0.70,
        target_effect=0.30,
        residual_sd=1.0,
        discovery_alpha=0.0001,
        replication_alpha=0.05,
        minimum_discovery_carriers=20,
        minimum_replication_carriers=10,
        maximum_metabolites_per_reaction=(
            int(network_config["maximum_metabolites_per_reaction"])
            if network_config.get("maximum_metabolites_per_reaction") is not None
            else None
        ),
        excluded_subsystem_patterns=network_config.get("excluded_subsystem_patterns", []),
    )
    included = set(candidate_config["included_review_priorities"])
    catalog = catalog.loc[catalog["review_priority"].isin(included)].copy()
    return mapping, catalog


def _analysis_frame(
    dataset,
    candidate: pd.Series,
    aggregation: str,
    group_column: str,
    minimum_group_size: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    genes = parse_gene_signature(str(candidate["genes"]))
    signature = build_expression_signature(
        dataset.blocks["transcriptomics"], genes, aggregation=aggregation
    )
    raw_frame = pd.DataFrame(
        {
            "outcome": dataset.blocks["metabolomics"][str(candidate["metabolite_b"])],
            "exposure": dataset.blocks["metabolomics"][str(candidate["metabolite_a"])],
            "modifier": signature,
            "group": dataset.sample_metadata[group_column].astype(str),
        }
    ).dropna()
    usable_groups = raw_frame["group"].value_counts()
    usable_groups = usable_groups[usable_groups >= minimum_group_size].index
    raw_frame = raw_frame.loc[raw_frame["group"].isin(usable_groups)].copy()
    frame = _standardize_frame(raw_frame)
    return frame, raw_frame


def _standardize_frame(frame: pd.DataFrame) -> pd.DataFrame:
    standardized = frame.copy()
    for column in ("outcome", "exposure", "modifier"):
        standardized[column] = zscore_within_groups(standardized[column], standardized["group"])
    return standardized.dropna()


def _stability_split_results(
    frame: pd.DataFrame,
    holdout_fraction: float,
    random_seed: int,
    control_group_slopes: bool,
) -> dict[str, object]:
    counts = frame["group"].value_counts()
    stratify = frame["group"] if counts.min() >= 4 else None
    train_index, holdout_index = train_test_split(
        np.arange(len(frame)),
        test_size=holdout_fraction,
        random_state=random_seed,
        stratify=stratify,
    )
    train = fit_coupling(
        frame.iloc[train_index],
        "outcome",
        "exposure",
        "modifier",
        group="group",
        control_group_slopes=control_group_slopes,
    )
    holdout = fit_coupling(
        frame.iloc[holdout_index],
        "outcome",
        "exposure",
        "modifier",
        group="group",
        control_group_slopes=control_group_slopes,
    )
    return {
        "stability_train_n": train.n,
        "stability_train_estimate": train.estimate,
        "stability_train_p_value": train.p_value,
        "stability_holdout_n": holdout.n,
        "stability_holdout_estimate": holdout.estimate,
        "stability_holdout_p_value": holdout.p_value,
        "holdout_same_direction": bool(train.estimate * holdout.estimate > 0),
    }


def _repeated_stability_results(
    frame: pd.DataFrame,
    holdout_fraction: float,
    random_seed: int,
    repetitions: int,
    reference_estimate: float,
    control_group_slopes: bool,
) -> dict[str, float | int]:
    holdout_estimates: list[float] = []
    holdout_p_values: list[float] = []
    for repetition in range(repetitions):
        result = _stability_split_results(
            frame,
            holdout_fraction,
            random_seed + repetition,
            control_group_slopes,
        )
        holdout_estimates.append(float(result["stability_holdout_estimate"]))
        holdout_p_values.append(float(result["stability_holdout_p_value"]))
    estimates = np.asarray(holdout_estimates)
    p_values = np.asarray(holdout_p_values)
    return {
        "stability_repetitions": repetitions,
        "stability_holdout_sign_agreement": float(np.mean(estimates * reference_estimate > 0)),
        "stability_holdout_p_lt_0_05_rate": float(np.mean(p_values < 0.05)),
        "stability_holdout_median_estimate": float(np.median(estimates)),
    }


def run_pipeline(config_path: Path, output_dir: Path, fetch: bool = True) -> dict[str, object]:
    config_path = config_path.resolve()
    project_root = config_path.parent.parent
    config = yaml.safe_load(config_path.read_text())
    dataset_config = config["dataset"]
    analysis_config = config["analysis"]
    random_seed = int(config["study"]["random_seed"])
    raw_dir = _resolve(project_root, dataset_config["raw_dir"])
    if fetch:
        fetch_ccle(raw_dir)

    mapping, catalog = _candidate_catalog(config, project_root)
    if catalog.empty:
        raise ValueError("No direct-reaction candidates were generated for the CCLE panel.")
    genes = sorted(
        {gene for signature in catalog["genes"] for gene in parse_gene_signature(str(signature))}
    )
    dataset = load_ccle(raw_dir, genes)

    records: list[dict[str, object]] = []
    frames: dict[str, pd.DataFrame] = {}
    raw_frames: dict[str, pd.DataFrame] = {}
    for _, candidate in catalog.iterrows():
        hypothesis_id = str(candidate["hypothesis_id"])
        try:
            frame, raw_frame = _analysis_frame(
                dataset,
                candidate,
                aggregation=str(analysis_config["expression_signature_aggregation"]),
                group_column=str(dataset_config["group_column"]),
                minimum_group_size=int(dataset_config["minimum_group_size"]),
            )
            result = fit_coupling(
                frame,
                outcome="outcome",
                exposure="exposure",
                modifier="modifier",
                group="group",
                covariance_type=str(analysis_config["covariance_type"]),
                control_group_slopes=bool(analysis_config["control_group_slopes"]),
            )
        except (ValueError, KeyError) as error:
            records.append(
                {
                    **candidate.to_dict(),
                    "analysis_status": "not_tested",
                    "analysis_error": str(error),
                }
            )
            continue
        records.append(
            {
                **candidate.to_dict(),
                **result.to_dict(),
                "analysis_status": "tested",
                "analysis_error": "",
            }
        )
        frames[hypothesis_id] = frame
        raw_frames[hypothesis_id] = raw_frame

    results = pd.DataFrame.from_records(records)
    tested_mask = results["analysis_status"].eq("tested")
    results["fdr_q_value"] = np.nan
    if tested_mask.any():
        results.loc[tested_mask, "fdr_q_value"] = multipletests(
            results.loc[tested_mask, "p_value"], method="fdr_bh"
        )[1]
    results["permutation_p_value"] = np.nan
    results["baseline_r2"] = np.nan
    results["interaction_r2"] = np.nan
    results["delta_r2"] = np.nan
    results["baseline_rmse"] = np.nan
    results["interaction_rmse"] = np.nan
    results["delta_rmse"] = np.nan

    ranked = results.loc[tested_mask].sort_values("p_value")
    follow_up_ranked = ranked.loc[
        ranked["fdr_q_value"] <= float(analysis_config["follow_up_q_value"])
    ]
    if follow_up_ranked.empty:
        follow_up_ranked = ranked
    permutation_rows = follow_up_ranked.head(int(analysis_config["permutation_top_k"])).index
    predictive_rows = follow_up_ranked.head(int(analysis_config["predictive_top_k"])).index
    for rank, index in enumerate(permutation_rows):
        hypothesis_id = str(results.at[index, "hypothesis_id"])
        frame = frames[hypothesis_id]
        results.at[index, "permutation_p_value"] = permutation_p_value(
            frame,
            outcome="outcome",
            exposure="exposure",
            modifier="modifier",
            observed_estimate=float(results.at[index, "estimate"]),
            group="group",
            permutations=int(analysis_config["permutations"]),
            random_seed=random_seed + rank,
            control_group_slopes=bool(analysis_config["control_group_slopes"]),
        )
    for rank, index in enumerate(predictive_rows):
        hypothesis_id = str(results.at[index, "hypothesis_id"])
        frame = raw_frames[hypothesis_id]
        prediction = compare_interaction_prediction(
            frame,
            outcome="outcome",
            exposure="exposure",
            modifier="modifier",
            group="group",
            folds=int(analysis_config["cross_validation_folds"]),
            random_seed=random_seed + rank,
            model_name=str(analysis_config["prediction_model"]),
            model_parameters=dict(analysis_config["prediction_model_parameters"]),
        )
        for key, value in prediction.to_dict().items():
            if key not in {"n", "folds"}:
                results.at[index, key] = value

    numeric_robustness_columns = [
        "stability_train_n",
        "stability_train_estimate",
        "stability_train_p_value",
        "stability_holdout_n",
        "stability_holdout_estimate",
        "stability_holdout_p_value",
        "stability_repetitions",
        "stability_holdout_sign_agreement",
        "stability_holdout_p_lt_0_05_rate",
        "stability_holdout_median_estimate",
        "winsorized_estimate",
        "winsorized_p_value",
        "influence_removed_n",
        "maximum_cooks_distance",
        "influence_trimmed_estimate",
        "influence_trimmed_p_value",
        "leave_one_group_out_runs",
        "leave_one_group_out_sign_agreement",
        "leave_one_group_out_min_estimate",
        "leave_one_group_out_max_estimate",
    ]
    for column in numeric_robustness_columns:
        results[column] = np.nan
    results["holdout_same_direction"] = pd.Series(pd.NA, index=results.index, dtype="boolean")
    follow_up = results.loc[
        tested_mask & (results["fdr_q_value"] <= float(analysis_config["false_discovery_rate"]))
    ].index
    if follow_up.empty:
        follow_up = ranked.head(min(5, len(ranked))).index
    for rank, index in enumerate(follow_up):
        hypothesis_id = str(results.at[index, "hypothesis_id"])
        raw_frame = raw_frames[hypothesis_id]
        frame = frames[hypothesis_id]
        split = _stability_split_results(
            raw_frame,
            holdout_fraction=float(analysis_config["holdout_fraction"]),
            random_seed=random_seed + rank,
            control_group_slopes=bool(analysis_config["control_group_slopes"]),
        )
        for key, value in split.items():
            results.at[index, key] = value

        winsorized = winsorize_within_groups(
            raw_frame,
            ["outcome", "exposure", "modifier"],
            "group",
            tail_fraction=float(analysis_config["winsorization_tail_fraction"]),
        )
        winsorized_result = fit_coupling(
            _standardize_frame(winsorized),
            "outcome",
            "exposure",
            "modifier",
            group="group",
            covariance_type=str(analysis_config["covariance_type"]),
            control_group_slopes=bool(analysis_config["control_group_slopes"]),
        )
        results.at[index, "winsorized_estimate"] = winsorized_result.estimate
        results.at[index, "winsorized_p_value"] = winsorized_result.p_value

        influence = refit_without_most_influential(
            frame,
            "outcome",
            "exposure",
            "modifier",
            group="group",
            exclusion_fraction=float(analysis_config["influence_exclusion_fraction"]),
            covariance_type=str(analysis_config["covariance_type"]),
            control_group_slopes=bool(analysis_config["control_group_slopes"]),
        )
        results.at[index, "influence_removed_n"] = influence.removed_n
        results.at[index, "maximum_cooks_distance"] = influence.maximum_cooks_distance
        results.at[index, "influence_trimmed_estimate"] = influence.estimate
        results.at[index, "influence_trimmed_p_value"] = influence.p_value

        leave_one_out = leave_one_group_out_summary(
            frame,
            "outcome",
            "exposure",
            "modifier",
            "group",
            reference_estimate=float(results.at[index, "estimate"]),
            covariance_type=str(analysis_config["covariance_type"]),
            control_group_slopes=bool(analysis_config["control_group_slopes"]),
        )
        for key, value in leave_one_out.items():
            results.at[index, key] = value

        repeated = _repeated_stability_results(
            raw_frame,
            holdout_fraction=float(analysis_config["holdout_fraction"]),
            random_seed=random_seed + rank * 100,
            repetitions=int(analysis_config["stability_repetitions"]),
            reference_estimate=float(results.at[index, "estimate"]),
            control_group_slopes=bool(analysis_config["control_group_slopes"]),
        )
        for key, value in repeated.items():
            results.at[index, key] = value

    results["proof_of_concept_pass"] = (
        (results["fdr_q_value"] <= float(analysis_config["false_discovery_rate"]))
        & (results["permutation_p_value"] <= 0.05)
        & (results["delta_r2"] > 0)
        & (results["winsorized_p_value"] <= 0.05)
        & (results["winsorized_estimate"] * results["estimate"] > 0)
        & (results["influence_trimmed_p_value"] <= 0.05)
        & (results["influence_trimmed_estimate"] * results["estimate"] > 0)
        & (
            results["leave_one_group_out_sign_agreement"]
            >= float(analysis_config["minimum_leave_one_group_out_sign_agreement"])
        )
        & (
            results["stability_holdout_sign_agreement"]
            >= float(analysis_config["minimum_holdout_sign_agreement"])
        )
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    mapping_path = output_dir / "ccle_metabolite_mapping.csv"
    catalog_path = output_dir / "ccle_candidate_hypotheses.csv"
    result_path = output_dir / "ccle_coupling_results.csv"
    mapping.to_csv(mapping_path, index=False)
    catalog.to_csv(catalog_path, index=False)
    results.sort_values(["analysis_status", "p_value"], na_position="last").to_csv(
        result_path, index=False
    )
    manifest: dict[str, object] = {
        "study": config["study"]["name"],
        "config": str(config_path),
        "config_sha256": _file_sha256(config_path),
        "dataset_provenance": dict(dataset.provenance),
        "aligned_samples": len(dataset.sample_ids),
        "mapped_metabolites": int(mapping["mapping_status"].eq("mapped").sum()),
        "direct_candidates": len(catalog),
        "tested_candidates": int(tested_mask.sum()),
        "fdr_significant": int(
            (results["fdr_q_value"] <= float(analysis_config["false_discovery_rate"])).sum()
        ),
        "proof_of_concept_passes": int(results["proof_of_concept_pass"].sum()),
        "model_system_boundary": (
            "CCLE cancer-cell-line transcript state is not human germline genotype evidence."
        ),
        "stability_boundary": (
            "Holdout analyses are post-selection sensitivity checks, not independent "
            "replication cohorts."
        ),
        "outputs": {
            "mapping": str(mapping_path.resolve()),
            "candidates": str(catalog_path.resolve()),
            "results": str(result_path.resolve()),
        },
    }
    manifest_path = output_dir / "ccle_proof_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=Path("config/ccle-proof.yaml"))
    parser.add_argument("--output-dir", type=Path, default=Path("artifacts/ccle-proof"))
    parser.add_argument("--skip-fetch", action="store_true")
    args = parser.parse_args()
    manifest = run_pipeline(args.config, args.output_dir.resolve(), fetch=not args.skip_fetch)
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
