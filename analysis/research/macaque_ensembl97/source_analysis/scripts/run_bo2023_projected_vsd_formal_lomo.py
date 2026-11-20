"""BrainTrace numerical helpers used by the Ensembl 97 reproduction route."""
from __future__ import annotations
from pathlib import Path
import pandas as pd
from core.bo2023_metadata import assert_unique_region_network_mapping, normalize_bo2023_network_labels, normalize_bo2023_region_labels
from core.expression_contract import validate_sample_id

def read_metadata(
    path: Path,
    sheet: str,
    region_col: str,
    network_col: str,
    monkey_col: str,
    label_map: dict[str, str] | None = None,
) -> pd.DataFrame:
    info = pd.read_excel(path, sheet_name=sheet, usecols=["No.", region_col, network_col, monkey_col], dtype=str)
    required = ["No.", region_col, network_col, monkey_col]
    for column in required:
        if info[column].isna().any() or info[column].astype(str).str.strip().eq("").any():
            raise ValueError(f"Metadata {column} contains missing or empty values; no samples were removed.")
    info["No."] = info["No."].map(validate_sample_id)
    for _, rows in info.groupby("No.", sort=False):
        if rows[required].astype(str).apply(lambda col: col.str.strip()).nunique().gt(1).any():
            raise ValueError("Conflicting metadata for duplicate sample_id; cannot define donor folds.")
    info["sample_id"] = info["No."].astype(str).str.strip()
    info["region_id"] = info[region_col].astype(str).str.strip()
    if label_map:
        info["region_id"] = info["region_id"].map(lambda value: label_map.get(value, value))
    info = normalize_bo2023_region_labels(info, columns=["region_id"])
    info["endpoint_label"] = info[network_col].astype(str).str.strip()
    info["monkey_id"] = info[monkey_col].astype(str).str.strip()
    info = normalize_bo2023_network_labels(info, network_col="endpoint_label")
    assert_unique_region_network_mapping(info, network_col="endpoint_label")
    info = info.drop_duplicates("sample_id")
    return info.copy()
