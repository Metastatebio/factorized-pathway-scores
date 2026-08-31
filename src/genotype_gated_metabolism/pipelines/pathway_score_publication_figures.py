"""Render journal-oriented figures for the pathway-score manuscript."""

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
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

COLORS = {
    "population": "#9AA5B5",
    "random": "#D6810B",
    "coarse": "#7450C8",
    "structure": "#16856A",
    "network": "#3659B8",
    "additive": "#6546B3",
    "interaction": "#9A62C7",
    "all": "#187FA6",
    "negative": "#B74D4D",
    "ink": "#16233A",
    "grid": "#DCE2EA",
    "paper": "#FFFFFF",
}


def _resolve(config_path: Path, value: str) -> Path:
    path = Path(value)
    return path if path.is_absolute() else (config_path.parent / path).resolve()


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _semantic_source_sha256(path: Path) -> str:
    """Hash stable manifest content without run-local timestamps and absolute paths."""
    if path.name != "manifest.json":
        return _sha256(path)
    payload = json.loads(path.read_text())
    for key in ("completed_at", "config", "protocol"):
        payload.pop(key, None)
    normalized = json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)
    return hashlib.sha256(normalized.encode()).hexdigest()


def _style() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.titlesize": 10,
            "axes.labelsize": 9,
            "axes.edgecolor": COLORS["grid"],
            "axes.linewidth": 0.8,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "figure.facecolor": COLORS["paper"],
            "axes.facecolor": COLORS["paper"],
            "savefig.facecolor": COLORS["paper"],
            "pdf.fonttype": 42,
        }
    )


def _panel(axis: plt.Axes, label: str, title: str) -> None:
    axis.set_title(f"{label}  {title}", loc="left", color=COLORS["ink"], fontweight="bold")
    axis.spines[["top", "right"]].set_visible(False)
    axis.grid(axis="x", color=COLORS["grid"], linewidth=0.6, alpha=0.8)
    axis.set_axisbelow(True)


def _save(figure: plt.Figure, base: Path, *, dpi: int) -> list[Path]:
    png = base.with_suffix(".png")
    pdf = base.with_suffix(".pdf")
    figure.savefig(png, dpi=dpi, bbox_inches="tight")
    figure.savefig(
        pdf,
        bbox_inches="tight",
        metadata={
            "Creator": "factorized-pathway-scores 1.0.0",
            "CreationDate": None,
            "ModDate": None,
        },
    )
    plt.close(figure)
    return [png, pdf]


def _study_design(path: Path, *, dpi: int) -> list[Path]:
    figure, axis = plt.subplots(figsize=(13, 6.2))
    axis.set_xlim(0, 1)
    axis.set_ylim(0, 1)
    axis.axis("off")
    columns = [0.015, 0.265, 0.515, 0.765]
    widths = [0.205, 0.205, 0.205, 0.22]
    headers = ["Public data", "Biological holdout", "Matched falsification", "Qualified claim"]
    for x, width, header in zip(columns, widths, headers, strict=True):
        axis.text(
            x + width / 2,
            0.945,
            header,
            ha="center",
            va="center",
            fontsize=10,
            color=COLORS["ink"],
            fontweight="bold",
        )
    lanes = [
        (
            0.72,
            "Frozen human test",
            "ST002081\n1,539 samples · 112 people",
            "Complete participant\nholdout",
            "Lipid families vs\nsize-matched random groups",
            "Predictive compression\n≠ family specificity",
            COLORS["coarse"],
        ),
        (
            0.47,
            "Locked human replication",
            "ST000818\n450 people · 15 populations",
            "Complete population-\ncategory holdout",
            "Chemical descriptors vs\nfixed-degree graph nulls",
            "Replicated structure-aware\nassay representation",
            COLORS["structure"],
        ),
        (
            0.22,
            "Reaction-topology test",
            "CCLE + HumanGEM\n913 lines · 60 targets",
            "Complete cancer-lineage\nholdout",
            "Direct neighborhoods vs\nproperty-matched features",
            "Compact topology signal;\nno general interaction gain",
            COLORS["network"],
        ),
    ]
    for y, lane, source, holdout, null, claim, color in lanes:
        axis.text(0.005, y + 0.105, lane, ha="left", va="bottom", color=color, fontweight="bold")
        for index, text in enumerate([source, holdout, null, claim]):
            box = FancyBboxPatch(
                (columns[index], y - 0.055),
                widths[index],
                0.14,
                boxstyle="round,pad=0.012,rounding_size=0.012",
                linewidth=1.2,
                edgecolor=color,
                facecolor="#F8FAFC",
            )
            axis.add_patch(box)
            axis.text(
                columns[index] + widths[index] / 2,
                y + 0.015,
                text,
                ha="center",
                va="center",
                color=COLORS["ink"],
                linespacing=1.25,
            )
            if index < 3:
                axis.add_patch(
                    FancyArrowPatch(
                        (columns[index] + widths[index] + 0.008, y + 0.015),
                        (columns[index + 1] - 0.008, y + 0.015),
                        arrowstyle="-|>",
                        mutation_scale=12,
                        linewidth=1.0,
                        color="#6D7888",
                    )
                )
    footer = FancyBboxPatch(
        (0.16, 0.015),
        0.68,
        0.08,
        boxstyle="round,pad=0.012,rounding_size=0.012",
        linewidth=0,
        facecolor="#EAF1F7",
    )
    axis.add_patch(footer)
    axis.text(
        0.5,
        0.055,
        "Specificity requires out-of-group prediction and superiority to an equally difficult null.",
        ha="center",
        va="center",
        color=COLORS["ink"],
        fontweight="bold",
    )
    return _save(figure, path / "figure-1-study-design", dpi=dpi)


def _human_figure(sources: dict[str, Path], path: Path, *, dpi: int) -> list[Path]:
    coarse = pd.read_csv(sources["st002081_metrics"])
    structural = pd.read_csv(sources["st002081_structural_metrics"])
    discovery_nulls = pd.read_csv(sources["st002081_nulls"])
    replication_nulls = pd.read_csv(sources["st000818_nulls"])
    sensitivity = pd.read_csv(sources["human_sensitivity"])
    mixing = json.loads(sources["human_graph_mixing"].read_text())["adjudication"]

    figure, axes = plt.subplots(2, 2, figsize=(12.2, 8.4), constrained_layout=True)
    selected = pd.concat(
        [
            coarse.loc[
                coarse["model"].isin(
                    ["population_mean", "family_score_ridge", "random_group_score_ridge", "all_visible_ridge"]
                )
            ],
            structural.loc[
                structural["model"].isin(
                    ["degree_preserving_random_structural_ridge", "structural_descriptor_ridge"]
                )
            ],
        ],
        ignore_index=True,
    )
    labels = {
        "population_mean": "Population mean",
        "family_score_ridge": "Lipid families",
        "random_group_score_ridge": "Size-matched random groups",
        "degree_preserving_random_structural_ridge": "Fixed-degree graph null",
        "structural_descriptor_ridge": "Chemical descriptors",
        "all_visible_ridge": "All visible metabolites",
    }
    colors = {
        "population_mean": COLORS["population"],
        "family_score_ridge": COLORS["coarse"],
        "random_group_score_ridge": COLORS["random"],
        "degree_preserving_random_structural_ridge": "#D9A536",
        "structural_descriptor_ridge": COLORS["structure"],
        "all_visible_ridge": COLORS["all"],
    }
    selected = selected.sort_values("row_weighted_rmse_sd", ascending=False)
    axes[0, 0].barh(
        [labels[value] for value in selected["model"]],
        selected["row_weighted_rmse_sd"],
        color=[colors[value] for value in selected["model"]],
    )
    axes[0, 0].set_xlabel("Participant-held-out RMSE (training-fold SD)")
    _panel(axes[0, 0], "A", "Prediction is not the same as specificity")

    for axis, data, panel_label, title, color in [
        (axes[0, 1], discovery_nulls, "B", "Discovery fixed-degree null ensemble", COLORS["structure"]),
        (axes[1, 0], replication_nulls, "C", "Locked population-held-out replication", COLORS["network"]),
    ]:
        ordered = data.sort_values("null_seed").reset_index(drop=True)
        y = np.arange(1, len(ordered) + 1)
        x = ordered["rmse_improvement_sd"].to_numpy()
        lower = x - ordered["rmse_ci_lower"].to_numpy()
        upper = ordered["rmse_ci_upper"].to_numpy() - x
        axis.errorbar(
            x,
            y,
            xerr=np.vstack([lower, upper]),
            fmt="o",
            ms=3.8,
            color=color,
            ecolor="#9BA7B6",
            elinewidth=0.8,
            capsize=1.8,
        )
        axis.axvline(0, color=COLORS["negative"], linewidth=1, linestyle="--")
        axis.set_ylabel("Graph-null realization")
        axis.set_xlabel("RMSE improvement over fixed-degree null (SD)")
        axis.set_yticks([1, 5, 10, 15, 20])
        _panel(axis, panel_label, title)

    settings = list(dict.fromkeys(sensitivity["setting"].astype(str)))
    positions = np.arange(len(settings))
    for offset, (dataset, color, marker) in enumerate(
        [("ST002081", COLORS["structure"], "o"), ("ST000818", COLORS["network"], "s")]
    ):
        block = sensitivity.loc[sensitivity["dataset"].eq(dataset)].set_index("setting").loc[settings]
        x = positions + (-0.12 if offset == 0 else 0.12)
        axes[1, 1].errorbar(
            x,
            block["rmse_improvement_sd"],
            yerr=np.vstack(
                [
                    block["rmse_improvement_sd"] - block["rmse_ci_lower"],
                    block["rmse_ci_upper"] - block["rmse_improvement_sd"],
                ]
            ),
            fmt=marker,
            ms=4,
            color=color,
            ecolor=color,
            alpha=0.85,
            capsize=2,
            label=dataset,
        )
    axes[1, 1].axhline(0, color=COLORS["negative"], linewidth=1, linestyle="--")
    axes[1, 1].set_xticks(positions, [value.replace("_", "\n") for value in settings], rotation=35, ha="right")
    axes[1, 1].set_ylabel("RMSE improvement over fixed-degree null (SD)")
    axes[1, 1].legend(frameon=False, loc="upper right")
    axes[1, 1].text(
        0.01,
        0.98,
        f"Null mixing: ≥{mixing['minimum_edge_replacement_fraction']:.1%} edges replaced\nmax pairwise Jaccard={mixing['maximum_pairwise_null_jaccard']:.3f}",
        transform=axes[1, 1].transAxes,
        ha="left",
        va="top",
        fontsize=8,
        color=COLORS["ink"],
    )
    _panel(axes[1, 1], "D", "Human sensitivity and null quality")
    return _save(figure, path / "figure-2-human-validation", dpi=dpi)


def _ccle_figure(sources: dict[str, Path], path: Path, *, dpi: int) -> list[Path]:
    metrics = pd.read_csv(sources["ccle_metrics"])
    sensitivity = pd.read_csv(sources["ccle_sensitivity"])
    dimension = pd.read_csv(sources["ccle_dimension_nulls"])
    matched = pd.read_csv(sources["ccle_property_nulls"])
    targets = pd.read_csv(sources["ccle_target_clusters"])
    mapping = pd.read_csv(sources["mapping_summary"])

    figure, axes = plt.subplots(2, 3, figsize=(15, 8.6), constrained_layout=True)
    model_order = [
        "population_mean",
        "random_factorized_ridge",
        "network_metabolites_ridge",
        "factorized_interaction_ridge",
        "network_additive_ridge",
        "all_metabolites_ridge",
    ]
    model_labels = {
        "population_mean": "Population mean",
        "random_factorized_ridge": "Random features",
        "network_metabolites_ridge": "Direct metabolites",
        "factorized_interaction_ridge": "Metabolites + GPR + interactions",
        "network_additive_ridge": "Metabolites + additive GPR",
        "all_metabolites_ridge": "All mapped metabolites",
    }
    model_colors = {
        "population_mean": COLORS["population"],
        "random_factorized_ridge": COLORS["random"],
        "network_metabolites_ridge": COLORS["network"],
        "factorized_interaction_ridge": COLORS["interaction"],
        "network_additive_ridge": COLORS["additive"],
        "all_metabolites_ridge": COLORS["all"],
    }
    block = metrics.set_index("model").loc[model_order].iloc[::-1]
    axes[0, 0].barh(
        [model_labels[value] for value in block.index],
        block["mean_equal_lineage_rmse_sd"],
        color=[model_colors[value] for value in block.index],
    )
    axes[0, 0].set_xlabel("Mean lineage-held-out RMSE (SD)")
    _panel(axes[0, 0], "A", "Reaction-neighborhood model ablation")

    rng = np.random.default_rng(42)
    for position, (data, label, color) in enumerate(
        [
            (dimension, "Dimension\nmatched", COLORS["random"]),
            (matched, "Degree + coverage\nmatched", COLORS["network"]),
        ]
    ):
        values = data["improvement_sd"].to_numpy()
        axes[0, 1].scatter(
            np.full(len(values), position) + rng.normal(0, 0.035, len(values)),
            values,
            s=22,
            color=color,
            alpha=0.8,
        )
        axes[0, 1].plot([position - 0.18, position + 0.18], [values.mean()] * 2, color=COLORS["ink"], linewidth=2)
    axes[0, 1].axhline(0, color=COLORS["negative"], linewidth=1, linestyle="--")
    axes[0, 1].set_xticks([0, 1], ["Dimension\nmatched", "Degree + coverage\nmatched"])
    axes[0, 1].set_ylabel("Random-minus-topology RMSE (SD)")
    _panel(axes[0, 1], "B", "All 40 random-feature realizations favor topology")

    subsystem = (
        targets.groupby("subsystem", as_index=False)
        .agg(targets=("target", "size"), mean_effect_sd=("effect_sd", "mean"))
        .sort_values("mean_effect_sd")
    )
    axes[0, 2].barh(
        subsystem["subsystem"],
        subsystem["mean_effect_sd"],
        color=np.where(subsystem["mean_effect_sd"].ge(0), COLORS["structure"], COLORS["negative"]),
    )
    axes[0, 2].axvline(0, color=COLORS["ink"], linewidth=0.8)
    axes[0, 2].set_xlabel("Mean property-matched effect (SD)")
    axes[0, 2].tick_params(axis="y", labelsize=6.5)
    _panel(axes[0, 2], "C", "Effects span 20 dominant HumanGEM subsystems")

    ordered = sensitivity.sort_values("rmse_improvement_vs_random_sd")
    y = np.arange(len(ordered))
    point = ordered["rmse_improvement_vs_random_sd"].to_numpy()
    axes[1, 0].errorbar(
        point,
        y,
        xerr=np.vstack([point - ordered["random_ci_lower"], ordered["random_ci_upper"] - point]),
        fmt="o",
        color=COLORS["network"],
        ecolor="#8E9BAA",
        capsize=2,
    )
    axes[1, 0].axvline(0, color=COLORS["negative"], linewidth=1, linestyle="--")
    axes[1, 0].set_yticks(y, [value.replace("_", " ") for value in ordered["setting"]])
    axes[1, 0].set_xlabel("Topology improvement over random (SD)")
    _panel(axes[1, 0], "D", "Topology passes all ten analysis settings")

    ordered = sensitivity.sort_values("rmse_improvement_vs_additive_sd")
    y = np.arange(len(ordered))
    point = ordered["rmse_improvement_vs_additive_sd"].to_numpy()
    axes[1, 1].errorbar(
        point,
        y,
        xerr=np.vstack([point - ordered["additive_ci_lower"], ordered["additive_ci_upper"] - point]),
        fmt="o",
        color=COLORS["interaction"],
        ecolor="#8E9BAA",
        capsize=2,
    )
    axes[1, 1].axvline(0, color=COLORS["ink"], linewidth=1)
    axes[1, 1].set_yticks(y, [value.replace("_", " ") for value in ordered["setting"]])
    axes[1, 1].set_xlabel("Interaction improvement over additive GPR (SD)")
    _panel(axes[1, 1], "E", "Interactions lose to additive GPR in every setting")

    mapping = mapping.copy()
    mapping["acceptance"] = mapping["accepted"] / mapping["total"]
    y = np.arange(len(mapping))
    axes[1, 2].barh(y, mapping["accepted"], color=COLORS["structure"], label="Accepted")
    axes[1, 2].barh(y, mapping["rejected"], left=mapping["accepted"], color="#D9DEE6", label="Rejected")
    axes[1, 2].set_yticks(y, mapping["domain"])
    axes[1, 2].set_xlabel("Features or signatures")
    axes[1, 2].legend(frameon=False, fontsize=8)
    _panel(axes[1, 2], "F", "Mapping denominators bound generalization")
    return _save(figure, path / "figure-3-ccle-validation", dpi=dpi)


def run(config_path: Path, output_dir: Path) -> dict[str, Any]:
    config_path = config_path.resolve()
    config = yaml.safe_load(config_path.read_text())
    sources = {
        key: _resolve(config_path, str(value))
        for key, value in config["sources"].items()
    }
    for key, source in sources.items():
        if not source.exists():
            raise FileNotFoundError(f"Figure source {key} does not exist: {source}")
    _style()
    output_dir.mkdir(parents=True, exist_ok=True)
    dpi = int(config["render"]["dpi"])
    outputs = []
    outputs.extend(_study_design(output_dir, dpi=dpi))
    outputs.extend(_human_figure(sources, output_dir, dpi=dpi))
    outputs.extend(_ccle_figure(sources, output_dir, dpi=dpi))
    project_root = config_path.parent.parent
    source_index = pd.DataFrame(
        [
            {
                "source": key,
                "path": path.relative_to(project_root).as_posix(),
                "sha256": _semantic_source_sha256(path),
            }
            for key, path in sorted(sources.items())
        ]
    )
    source_index_path = output_dir / "figure-source-index.csv"
    source_index.to_csv(source_index_path, index=False)
    outputs.append(source_index_path)
    manifest: dict[str, Any] = {
        "analysis_id": str(config["analysis_id"]),
        "completed_at": datetime.now(UTC).isoformat(),
        "config": str(config_path),
        "config_sha256": _sha256(config_path),
        "implementation_sha256": _sha256(Path(inspect.getfile(run))),
        "figures": 3,
        "source_files": len(sources),
        "output_sha256": {value.name: _sha256(value) for value in outputs},
    }
    manifest_path = output_dir / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, allow_nan=False) + "\n")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config", type=Path, default=Path("config/pathway-score-publication-figures.yaml")
    )
    parser.add_argument(
        "--output", type=Path, default=Path("artifacts/pathway-score-publication-figures")
    )
    args = parser.parse_args()
    print(json.dumps(run(args.config, args.output), indent=2))


if __name__ == "__main__":
    main()
