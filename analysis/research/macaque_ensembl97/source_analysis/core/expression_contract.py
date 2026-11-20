"""BrainTrace numerical helpers used by the Ensembl 97 reproduction route."""
from __future__ import annotations
from typing import Any
import pandas as pd

def validate_sample_id(value: Any) -> str:
    """Do not silently trim IDs: the same value must identify every stored row."""
    if value is None or pd.isna(value):
        raise ValueError("sample_id is required.")
    value = str(value)
    if not value or value != value.strip():
        raise ValueError("sample_id must be non-empty and have no leading/trailing whitespace.")
    return value
