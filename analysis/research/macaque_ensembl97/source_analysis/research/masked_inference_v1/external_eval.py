"""BrainTrace numerical helpers used by the Ensembl 97 reproduction route."""
from __future__ import annotations
from .common import require_hash
from pathlib import Path
import numpy as np

def _stream_observed_counts(coverage, spec: dict, mapping: list[dict]) -> tuple[list[str], np.ndarray, np.ndarray, list[str]]:
    """All original counts checked; return only retained lnCPM in RAM."""
    if spec["cohort_id"] not in {"AHBA_counts:H0351_2001", "AHBA_counts:H0351_2002", "GSE189919_counts"}:
        raise ValueError("Numeric parsing/inference is restricted to AHBA and GSE counts")
    if spec["unit"] != "raw_counts" or spec.get("complete_original_library_certified") is not True:
        raise ValueError("A certified complete supplied raw-count library is required")
    path = Path(spec["path"]); require_hash(path, spec["sha256"])
    axis, assays = coverage.read_axis(path, delimiter=spec["delimiter"], zip_member=spec.get("zip_member"), sample_ids=spec.get("sample_ids"))
    if (len(axis), len(assays)) != (spec["expected_genes"], spec["expected_profiles"]):
        raise ValueError("Count-library axis differs from frozen specification")
    if spec.get("aggregation"):
        if spec["aggregation"] != "AHBA_donor_tissue_technical_counts_sum":
            raise ValueError("Only fixed AHBA technical-count summation is allowed")
        profiles, groups = coverage.ahba_groups(spec["metadata_rows"], assays)
    else:
        profiles, groups = list(assays), None
    row_to_feature = {record["platform_row"]: i for i, record in enumerate(mapping) if not record["reason"]}
    if len(row_to_feature) != sum(not record["reason"] for record in mapping):
        raise ValueError("A real target row cannot serve more than one source feature")
    panel_counts = np.full((len(mapping), len(profiles)), np.nan, dtype=np.float64)
    totals = np.zeros(len(profiles), dtype=np.float64)
    with coverage.table_rows(path, spec["delimiter"], spec.get("zip_member")) as rows:
        if spec.get("zip_member") is None and next(rows)[1:] != assays:
            raise ValueError("Sample axis changed during streaming")
        seen = 0
        for seen, row in enumerate(rows, 1):
            if seen > len(axis) or len(row) != len(assays) + 1 or row[0] != axis[seen - 1]:
                raise ValueError("Count-library row order/width changed")
            try:
                values = np.asarray([float(v) for v in row[1:]], dtype=np.float64)
            except (ValueError, OverflowError) as error:
                raise ValueError("Unparseable original count-library value") from error
            if not np.isfinite(values).all() or np.any(values < 0):
                raise ValueError("Nonfinite/negative original count-library value")
            collapsed = values if groups is None else np.asarray([values[indices].sum(dtype=np.float64) for indices in groups], dtype=np.float64)
            totals += collapsed
            if not np.isfinite(collapsed).all() or not np.isfinite(totals).all():
                raise ValueError("Count/library aggregation overflow")
            if row[0] in row_to_feature:
                panel_counts[row_to_feature[row[0]]] = collapsed
        if seen != len(axis):
            raise ValueError("Count-library row axis shortened")
    available = np.asarray([not record["reason"] for record in mapping], dtype=bool)
    if not np.isfinite(panel_counts[available]).all():
        raise ValueError("Mapped real count rows absent during streaming")
    x = np.full((len(profiles), len(mapping)), np.nan, dtype=np.float32)
    library_reasons = []
    for j, total in enumerate(totals):
        if total <= 0:
            library_reasons.append("zero_full_library")
            continue
        library_reasons.append("")
        with np.errstate(over="ignore", invalid="ignore"):
            x[j, available] = np.log1p(panel_counts[available, j] / total * 1e6).astype(np.float32)
    require_hash(path, spec["sha256"])
    return profiles, x, totals, library_reasons
