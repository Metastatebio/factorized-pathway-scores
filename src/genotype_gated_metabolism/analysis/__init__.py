"""Statistical and predictive cross-omics analyses."""

from .coupling import CouplingResult, fit_coupling, permutation_p_value
from .prediction import PredictionComparison, compare_interaction_prediction

__all__ = [
    "CouplingResult",
    "PredictionComparison",
    "compare_interaction_prediction",
    "fit_coupling",
    "permutation_p_value",
]
