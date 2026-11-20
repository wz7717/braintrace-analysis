"""BrainTrace numerical helpers used by the Ensembl 97 reproduction route."""
from __future__ import annotations
import numpy as np

def trace_corr(A: np.ndarray, b: np.ndarray) -> np.ndarray:
    A = np.asarray(A, dtype=float)
    b = np.asarray(b, dtype=float)
    if A.ndim != 2 or b.ndim != 1 or A.shape[0] != b.size:
        raise ValueError("Pearson reference and query dimensions do not match")
    if not np.isfinite(A).all() or not np.isfinite(b).all():
        raise ValueError("Pearson reference and query must contain only finite values")
    b0 = b - b.mean()
    A0 = A - A.mean(axis=0, keepdims=True)
    num = (A0 * b0[:, None]).sum(axis=0)
    den = np.sqrt((A0 ** 2).sum(axis=0) * (b0 ** 2).sum() + 1e-12)
    corr = num / den
    if not np.isfinite(corr).all():
        raise ValueError("Pearson calculation produced nonfinite scores")
    return corr
