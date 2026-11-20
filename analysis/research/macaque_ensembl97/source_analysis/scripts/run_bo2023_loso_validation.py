"""BrainTrace numerical helpers used by the Ensembl 97 reproduction route."""
from __future__ import annotations
import numpy as np

def correlation_scores(reference: np.ndarray, sample: np.ndarray, rows: np.ndarray | None = None) -> np.ndarray:
    from core.models import trace_corr
    ref = reference if rows is None else reference[rows, :]
    vec = sample if rows is None else sample[rows]
    return trace_corr(ref, vec)
