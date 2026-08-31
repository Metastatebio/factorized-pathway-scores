"""Leakage-controlled predictive comparison for interaction features."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any

import numpy as np
import pandas as pd
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.model_selection import KFold, StratifiedKFold

from .estimators import build_regressor


@dataclass(frozen=True)
class PredictionComparison:
    n: int
    folds: int
    baseline_r2: float
    interaction_r2: float
    delta_r2: float
    baseline_rmse: float
    interaction_rmse: float
    delta_rmse: float

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


def _prediction_design(
    frame: pd.DataFrame, exposure: str, modifier: str, group: str | None
) -> tuple[pd.DataFrame, pd.DataFrame]:
    baseline = frame[[exposure, modifier]].copy()
    interaction = baseline.copy()
    interaction[f"{exposure}:x:{modifier}"] = frame[exposure] * frame[modifier]
    if group is not None:
        dummies = pd.get_dummies(frame[group].astype(str), prefix=group, drop_first=False)
        dummies = dummies.astype(float)
        baseline = pd.concat([baseline, dummies], axis=1)
        interaction = pd.concat([interaction, dummies], axis=1)
        for column in dummies:
            baseline[f"{exposure}:x:{column}"] = frame[exposure] * dummies[column]
            baseline[f"{modifier}:x:{column}"] = frame[modifier] * dummies[column]
            interaction[f"{exposure}:x:{column}"] = baseline[f"{exposure}:x:{column}"]
            interaction[f"{modifier}:x:{column}"] = baseline[f"{modifier}:x:{column}"]
    return baseline.astype(float), interaction.astype(float)


def compare_interaction_prediction(
    data: pd.DataFrame,
    outcome: str,
    exposure: str,
    modifier: str,
    group: str | None = None,
    folds: int = 5,
    random_seed: int = 20260829,
    model_name: str = "ridge",
    model_parameters: dict[str, Any] | None = None,
) -> PredictionComparison:
    """Compare out-of-fold prediction with and without the A:G feature."""
    required = [outcome, exposure, modifier, *([group] if group else [])]
    frame = data.loc[:, required].replace([np.inf, -np.inf], np.nan).dropna().copy()
    if len(frame) < max(40, folds * 5):
        raise ValueError("Insufficient complete samples for predictive evaluation.")
    baseline, interaction = _prediction_design(frame, exposure, modifier, group)
    target = frame[outcome].to_numpy(dtype=float)

    if group is not None and frame[group].value_counts().min() >= folds:
        splitter = StratifiedKFold(n_splits=folds, shuffle=True, random_state=random_seed)
        splits = splitter.split(baseline, frame[group].astype(str))
    else:
        splitter = KFold(n_splits=folds, shuffle=True, random_state=random_seed)
        splits = splitter.split(baseline)

    baseline_prediction = np.full(len(frame), np.nan)
    interaction_prediction = np.full(len(frame), np.nan)
    for fold, (train, test) in enumerate(splits):
        baseline_model = build_regressor(
            model_name, model_parameters, random_seed=random_seed + fold
        )
        interaction_model = build_regressor(
            model_name, model_parameters, random_seed=random_seed + fold
        )
        baseline_model.fit(baseline.iloc[train], target[train])
        interaction_model.fit(interaction.iloc[train], target[train])
        baseline_prediction[test] = baseline_model.predict(baseline.iloc[test])
        interaction_prediction[test] = interaction_model.predict(interaction.iloc[test])

    baseline_r2 = float(r2_score(target, baseline_prediction))
    interaction_r2 = float(r2_score(target, interaction_prediction))
    baseline_rmse = float(np.sqrt(mean_squared_error(target, baseline_prediction)))
    interaction_rmse = float(np.sqrt(mean_squared_error(target, interaction_prediction)))
    return PredictionComparison(
        n=len(frame),
        folds=folds,
        baseline_r2=baseline_r2,
        interaction_r2=interaction_r2,
        delta_r2=interaction_r2 - baseline_r2,
        baseline_rmse=baseline_rmse,
        interaction_rmse=interaction_rmse,
        delta_rmse=baseline_rmse - interaction_rmse,
    )
