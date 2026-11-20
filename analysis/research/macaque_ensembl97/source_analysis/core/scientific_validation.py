"""BrainTrace numerical helpers used by the Ensembl 97 reproduction route."""
from __future__ import annotations
import numpy as np
import pandas as pd
from core.signal_evidence import panel_signal_evidence, ranking_is_undifferentiated

def signal_reason(genes: list[str], observed: pd.Series, vector: np.ndarray, scores: np.ndarray) -> str | None:
    signal = panel_signal_evidence(genes, observed, vector)
    if not signal['evaluable']:
        return str(signal['reason'])
    if ranking_is_undifferentiated(scores):
        return 'scores_do_not_distinguish_candidates'
    return None
