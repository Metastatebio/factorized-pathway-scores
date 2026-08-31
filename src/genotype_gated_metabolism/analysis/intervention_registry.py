"""Conservative metadata audit for public human intervention metabolomics studies."""

from __future__ import annotations

import re
from collections.abc import Iterable, Mapping

import pandas as pd

from ..metabolite_mapping import normalize_metabolite_name

SYSTEMIC_SAMPLE_TERMS = ("blood", "plasma", "serum", "urine")
TREATMENT_FACTOR_TERMS = (
    "arm",
    "diet",
    "dose",
    "group",
    "intervention",
    "status",
    "supplement",
    "treatment",
)
TIME_FACTOR_TERMS = ("day", "hour", "period", "time", "visit", "week")


def parse_factor_string(value: object) -> dict[str, str]:
    """Parse `Factor:Value | Factor:Value` metadata without losing inner colons."""
    if pd.isna(value):
        return {}
    factors: dict[str, str] = {}
    for item in str(value).split(" | "):
        name, separator, level = item.partition(":")
        if not separator or not name.strip():
            continue
        factors[name.strip()] = level.strip()
    return factors


def _factor_levels(factors: pd.DataFrame) -> dict[str, set[str]]:
    levels: dict[str, set[str]] = {}
    if "factors" not in factors:
        return levels
    for encoded in factors["factors"]:
        for name, value in parse_factor_string(encoded).items():
            levels.setdefault(name, set()).add(value)
    return levels


def _matching_factor_names(levels: Mapping[str, set[str]], terms: Iterable[str]) -> list[str]:
    lowered_terms = tuple(term.lower() for term in terms)
    return [
        name
        for name in levels
        if any(term in name.lower() for term in lowered_terms)
    ]


def title_is_intervention_candidate(title: object, keywords: Iterable[str]) -> bool:
    """Return whether a title merits detailed intervention metadata retrieval."""
    normalized = str(title).lower()
    return any(keyword.lower() in normalized for keyword in keywords)


def classify_intervention(title: object) -> str:
    """Assign one broad intervention class from the study title."""
    text = str(title).lower()
    if re.search(r"\b(?:fasting|fasted)\b|time[ -]restricted", text):
        return "fasting_or_timing"
    classes = (
        ("weight_loss", ("weight loss", "calorie restriction", "caloric restriction")),
        ("exercise", ("exercise", "training", "physical activity")),
        ("diet_composition", ("diet", "feeding", "food", "protein", "carbohydrate", "fat")),
        ("supplement", ("supplement", "vitamin", "amino acid", "probiotic", "prebiotic")),
        ("pharmacologic", ("drug", "aspirin", "treatment", "therapy", "placebo")),
        ("surgery", ("surgery", "bypass", "gastrectomy")),
    )
    for label, terms in classes:
        if any(term in text for term in terms):
            return label
    return "other_intervention"


def normalize_trial_key(
    project_doi: object,
    project_id: object,
    study_id: object,
) -> str:
    """Return a conservative project-level key for independent-trial counting."""
    doi = str(project_doi).strip().lower()
    doi = re.sub(r"^(?:doi:\s*|https?://(?:dx\.)?doi\.org/)", "", doi).strip()
    if doi:
        return f"doi:{doi}"
    project = str(project_id).strip().lower()
    if project:
        return f"project:{project}"
    return f"study:{str(study_id).strip().lower()}"


def _joined_unique(values: Iterable[object]) -> str:
    observed = {
        item.strip()
        for value in values
        for item in str(value).split(";")
        if item.strip()
    }
    return ";".join(sorted(observed))


def aggregate_trial_registry(
    studies: pd.DataFrame,
    *,
    core_markers: set[str],
) -> pd.DataFrame:
    """Collapse assay substudies into conservative independent project/trial units."""
    if studies.empty:
        return pd.DataFrame()
    required = {"study_id", "project_doi", "project_id", "eligibility_status"}
    missing = sorted(required - set(studies.columns))
    if missing:
        raise ValueError(f"Trial aggregation is missing columns: {', '.join(missing)}")
    frame = studies.copy()
    frame["trial_key"] = frame.apply(
        lambda row: normalize_trial_key(
            row["project_doi"], row["project_id"], row["study_id"]
        ),
        axis=1,
    )
    records: list[dict[str, object]] = []
    eligibility_rank = {
        "priority_protocol_review": 0,
        "protocol_review_required": 1,
        "insufficient_metadata": 2,
    }
    for trial_key, group in frame.groupby("trial_key", sort=True):
        ranked = group.assign(
            _eligibility_rank=group["eligibility_status"].map(eligibility_rank).fillna(3)
        ).sort_values(
            ["_eligibility_rank", "ppm1k_core_marker_count", "sample_count", "study_id"],
            ascending=[True, False, False, True],
            kind="stable",
        )
        representative = ranked.iloc[0]
        metadata_ready = bool(group["metadata_ready"].astype(bool).any())
        randomization_evidence = bool(
            group["randomization_metadata_evidence"].astype(bool).any()
        )
        if not metadata_ready:
            eligibility = "insufficient_metadata"
        elif randomization_evidence:
            eligibility = "priority_protocol_review"
        else:
            eligibility = "protocol_review_required"
        markers = {
            marker
            for encoded in group["ppm1k_markers"]
            for marker in str(encoded).split(";")
            if marker
        }
        project_dois = [str(value).strip() for value in group["project_doi"] if str(value).strip()]
        project_ids = [str(value).strip() for value in group["project_id"] if str(value).strip()]
        records.append(
            {
                "trial_key": trial_key,
                "project_doi": min(project_dois) if project_dois else "",
                "project_id": min(project_ids) if project_ids else "",
                "representative_study_id": representative["study_id"],
                "study_ids": _joined_unique(group["study_id"]),
                "study_count": len(group),
                "multiple_analytical_substudies": len(group) > 1,
                "study_titles": " | ".join(sorted(set(group["study_title"].astype(str)))),
                "intervention_classes": _joined_unique(group["intervention_class"]),
                "sample_count_max": int(group["sample_count"].max()),
                "total_subjects_max": int(group["total_subjects"].max()),
                "named_metabolites_max": int(group["named_metabolites"].max()),
                "systemic_sample": bool(group["systemic_sample"].astype(bool).any()),
                "treatment_levels": _joined_unique(group["treatment_levels"]),
                "time_levels": _joined_unique(group["time_levels"]),
                "control_detected": bool(group["control_detected"].astype(bool).any()),
                "randomization_metadata_evidence": randomization_evidence,
                "metadata_ready": metadata_ready,
                "eligibility_status": eligibility,
                "raw_data_available": _joined_unique(group["raw_data_available"]),
                "ppm1k_marker_count": len(markers),
                "ppm1k_core_marker_count": len(markers & core_markers),
                "ppm1k_markers": ";".join(sorted(markers)),
                "study_page_url": representative["study_page_url"],
                "study_url": representative["study_url"],
                "hard_limit": (
                    "A shared repository DOI/project is counted once. This does not prove that "
                    "all assay substudies share participants or one randomized contrast."
                ),
            }
        )
    return pd.DataFrame.from_records(records).sort_values(
        "trial_key", kind="stable"
    ).reset_index(drop=True)


def summarize_design(
    study: Mapping[str, object],
    factors: pd.DataFrame,
    metabolites: pd.DataFrame,
    *,
    minimum_samples: int,
    minimum_named_metabolites: int,
    randomization_keywords: Iterable[str],
    protocol_text: object = "",
) -> dict[str, object]:
    """Summarize what study metadata establish and what still needs review."""
    title = str(study.get("study_title", ""))
    levels = _factor_levels(factors)
    treatment_names = _matching_factor_names(levels, TREATMENT_FACTOR_TERMS)
    time_names = _matching_factor_names(levels, TIME_FACTOR_TERMS)
    treatment_levels = sorted(
        set().union(*(levels[name] for name in treatment_names)) if treatment_names else set()
    )
    time_levels = sorted(
        set().union(*(levels[name] for name in time_names)) if time_names else set()
    )
    sample_sources = sorted(
        {
            str(value).strip()
            for value in factors.get("sample_source", pd.Series(dtype=object)).dropna()
            if str(value).strip()
        }
    )
    systemic = any(
        term in source.lower() for source in sample_sources for term in SYSTEMIC_SAMPLE_TERMS
    )
    raw_names = metabolites.get(
        "metabolite_name", pd.Series("", index=metabolites.index, dtype=object)
    ).fillna("").astype(str)
    refmet_names = metabolites.get(
        "refmet_name", pd.Series("", index=metabolites.index, dtype=object)
    ).fillna("").astype(str)
    preferred_names = refmet_names.mask(refmet_names.eq(""), raw_names)
    named_metabolites = int(preferred_names.loc[preferred_names.ne("")].nunique())
    raw_sample_count = pd.to_numeric(study.get("number_of_samples", 0), errors="coerce")
    sample_count = int(raw_sample_count) if pd.notna(raw_sample_count) else 0
    title_randomization_evidence = any(
        keyword.lower() in title.lower() for keyword in randomization_keywords
    )
    repository_randomization_evidence = any(
        keyword.lower() in str(protocol_text).lower() for keyword in randomization_keywords
    )
    randomization_evidence = bool(
        title_randomization_evidence or repository_randomization_evidence
    )
    control_detected = any(
        any(term in level.lower() for term in ("control", "placebo", "usual", "sham"))
        for level in treatment_levels
    )
    metadata_ready = bool(
        systemic
        and sample_count >= minimum_samples
        and named_metabolites >= minimum_named_metabolites
        and len(treatment_levels) >= 2
    )
    if not metadata_ready:
        eligibility = "insufficient_metadata"
    elif randomization_evidence:
        eligibility = "priority_protocol_review"
    else:
        eligibility = "protocol_review_required"
    return {
        "study_id": study.get("study_id", ""),
        "study_title": title,
        "intervention_class": classify_intervention(title),
        "sample_count": sample_count,
        "sample_sources": ";".join(sample_sources),
        "systemic_sample": systemic,
        "factor_names": ";".join(sorted(levels)),
        "treatment_factor_names": ";".join(treatment_names),
        "treatment_levels": ";".join(treatment_levels),
        "treatment_level_count": len(treatment_levels),
        "time_factor_names": ";".join(time_names),
        "time_levels": ";".join(time_levels),
        "time_level_count": len(time_levels),
        "control_detected": control_detected,
        "title_randomization_evidence": title_randomization_evidence,
        "repository_randomization_evidence": repository_randomization_evidence,
        "randomization_metadata_evidence": randomization_evidence,
        "named_metabolites": named_metabolites,
        "metadata_ready": metadata_ready,
        "eligibility_status": eligibility,
        "hard_limit": (
            "Repository metadata do not establish randomization, participant linkage, "
            "analysis-population integrity, or endpoint usability; protocol review is required."
        ),
    }


def marker_matches(
    study_id: str,
    metabolites: pd.DataFrame,
    marker_aliases: Mapping[str, Iterable[str]],
) -> pd.DataFrame:
    """Return exact normalized-name matches to a frozen marker alias dictionary."""
    alias_index: dict[str, str] = {}
    for marker, aliases in marker_aliases.items():
        for alias in {marker, *aliases}:
            normalized = normalize_metabolite_name(alias)
            if normalized and normalized in alias_index and alias_index[normalized] != marker:
                raise ValueError(f"Marker alias {alias!r} maps to multiple frozen markers.")
            alias_index[normalized] = marker
    records: list[dict[str, object]] = []
    for row in metabolites.to_dict(orient="records"):
        source_names = [row.get("refmet_name", ""), row.get("metabolite_name", "")]
        matched = {
            alias_index[normalized]
            for value in source_names
            if (normalized := normalize_metabolite_name(value)) in alias_index
        }
        for marker in sorted(matched):
            records.append(
                {
                    "study_id": study_id,
                    "marker": marker,
                    "analysis_id": row.get("analysis_id", ""),
                    "source_metabolite_name": row.get("metabolite_name", ""),
                    "refmet_name": row.get("refmet_name", ""),
                }
            )
    return pd.DataFrame.from_records(
        records,
        columns=[
            "study_id",
            "marker",
            "analysis_id",
            "source_metabolite_name",
            "refmet_name",
        ],
    ).drop_duplicates(["study_id", "marker", "analysis_id"])
