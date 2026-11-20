"""BrainTrace numerical helpers used by the Ensembl 97 reproduction route."""
from __future__ import annotations
from typing import Any, Sequence
import numpy as np
import pandas as pd

NUMERICAL_SIGNAL_TOLERANCE = 1e-12


def _finite_standard_deviation(values: np.ndarray) -> float | None:
    if not values.size or not np.isfinite(values).all():
        return None
    scale = float(np.max(np.abs(values)))
    result = 0.0 if scale == 0 else float(np.std(values / scale) * scale)
    return result if np.isfinite(result) else None


def panel_signal_evidence(
    genes: Sequence[str],
    observed: pd.Series,
    scoring_values: np.ndarray,
) -> dict[str, Any]:
    """Describe actual panel observations separately from zero-filled scores.

    This only rejects undefined or numerically constant Pearson evidence. It is
    not a calibrated coverage, detection or out-of-distribution threshold.
    """
    panel = pd.Index(genes).astype(str)
    present = panel.isin(observed.index)
    values = observed.reindex(panel[present]).to_numpy(dtype=float)
    vector = np.asarray(scoring_values, dtype=float)
    finite = np.isfinite(values)
    diagnostics: dict[str, Any] = {
        "n_panel_genes": int(len(panel)),
        "n_observed_genes": int(present.sum()),
        "n_finite_observed_genes": int(finite.sum()),
        "n_nonzero_observed_genes": int(np.count_nonzero(values[finite])),
        "observed_standard_deviation": _finite_standard_deviation(values),
        "scoring_standard_deviation": _finite_standard_deviation(vector),
    }
    reason = None
    if not present.any():
        reason = "no_observed_panel_genes"
    elif not finite.all():
        reason = "nonfinite_observed_panel_values"
    elif len(values) < 2:
        reason = "fewer_than_two_observed_panel_genes"
    elif float(np.ptp(values)) <= NUMERICAL_SIGNAL_TOLERANCE:
        reason = "constant_observed_panel_values"
    elif not np.isfinite(vector).all():
        reason = "nonfinite_scoring_values"
    elif vector.size < 2 or float(np.ptp(vector)) <= NUMERICAL_SIGNAL_TOLERANCE:
        reason = "constant_scoring_vector"
    diagnostics["evaluable"] = reason is None
    diagnostics["reason"] = reason
    return diagnostics


def ranking_is_undifferentiated(scores: np.ndarray) -> bool:
    values = np.asarray(scores, dtype=float)
    return (
        values.size < 2
        or not np.isfinite(values).all()
        or float(np.ptp(values)) <= NUMERICAL_SIGNAL_TOLERANCE
    )
