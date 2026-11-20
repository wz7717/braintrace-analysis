"""BrainTrace numerical helpers used by the Ensembl 97 reproduction route."""
from __future__ import annotations
import numpy as np

def build_region_reference(
    values: np.ndarray,
    regions: list[str],
    training: dict[str, np.ndarray],
) -> np.ndarray:
    return np.column_stack(
        [values[:, training[region]].mean(axis=1, dtype=np.float64) for region in regions]
    ).astype(np.float32)


def zscore(scores: np.ndarray) -> np.ndarray:
    scores = np.asarray(scores, dtype=np.float64)
    std = float(scores.std())
    if std <= 1e-12:
        return np.zeros_like(scores)
    return (scores - float(scores.mean())) / std
