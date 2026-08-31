"""Configuration-driven scikit-learn regressors for cross-omic prediction."""

from __future__ import annotations

from collections.abc import Mapping

from sklearn.base import RegressorMixin
from sklearn.ensemble import HistGradientBoostingRegressor, RandomForestRegressor
from sklearn.linear_model import ElasticNet, Ridge
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler


def build_regressor(
    name: str,
    parameters: Mapping[str, object] | None = None,
    random_seed: int = 20260829,
) -> RegressorMixin:
    """Build a supported regressor from a stable name and parameter mapping."""
    normalized = name.strip().lower().replace("-", "_")
    options = dict(parameters or {})
    if normalized == "ridge":
        return make_pipeline(StandardScaler(), Ridge(**options))
    if normalized == "elastic_net":
        options.setdefault("random_state", random_seed)
        options.setdefault("max_iter", 10000)
        return make_pipeline(StandardScaler(), ElasticNet(**options))
    if normalized == "random_forest":
        options.setdefault("random_state", random_seed)
        options.setdefault("n_jobs", -1)
        return RandomForestRegressor(**options)
    if normalized == "hist_gradient_boosting":
        options.setdefault("random_state", random_seed)
        return HistGradientBoostingRegressor(**options)
    supported = "ridge, elastic_net, random_forest, hist_gradient_boosting"
    raise ValueError(f"Unknown regressor '{name}'. Supported regressors: {supported}.")
