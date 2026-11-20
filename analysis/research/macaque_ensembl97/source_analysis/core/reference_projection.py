"""BrainTrace numerical helpers used by the Ensembl 97 reproduction route."""
from __future__ import annotations
import re
from pathlib import Path
import pandas as pd

def read_gene_map(path: Path) -> pd.DataFrame:
    mapping = pd.read_csv(path, usecols=["Gene.stable.ID", "Gene.name"], dtype=str)
    mapping = mapping.rename(columns={"Gene.stable.ID": "gene_id", "Gene.name": "gene_symbol"})
    if mapping["gene_id"].isna().any() or mapping["gene_id"].str.strip().eq("").any():
        raise ValueError("Gene map contains missing gene IDs")
    mapping["gene_id"] = mapping["gene_id"].str.strip()
    mapping["gene_symbol"] = mapping["gene_symbol"].map(clean_excel_date_gene_symbol)
    mapping.loc[mapping["gene_symbol"].isin(["", "nan", "None"]), "gene_symbol"] = pd.NA
    mapping = mapping.drop_duplicates(["gene_id", "gene_symbol"])
    if mapping.groupby("gene_id")["gene_symbol"].nunique(dropna=False).gt(1).any():
        raise ValueError("Gene map assigns conflicting symbols to a gene ID")
    return mapping


EXCEL_DATE_SYMBOL_RE = re.compile(r"^\d{4}-(\d{2})-(\d{2})(?: 00:00:00)?$")


def clean_excel_date_gene_symbol(value: object) -> str:
    text = "" if pd.isna(value) else str(value).strip()
    match = EXCEL_DATE_SYMBOL_RE.match(text)
    if not match:
        return text
    month = int(match.group(1))
    day = int(match.group(2))
    prefix = {3: "MARCH", 9: "SEPT", 12: "DEC"}.get(month)
    if not prefix:
        return text
    return f"{prefix}{day}"
