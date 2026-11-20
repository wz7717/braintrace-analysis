"""BrainTrace numerical helpers used by the Ensembl 97 reproduction route."""
from __future__ import annotations
from __future__ import annotations
import csv
import gzip
import hashlib
import io
import re
import zipfile
from collections import Counter, defaultdict
from contextlib import contextmanager
from pathlib import Path
from typing import Iterable, Sequence

HUMAN_ID = re.compile(r"ENSG[0-9]+(?:\.[0-9]+)?\Z")


MACAQUE_ID = re.compile(r"ENSMFAG[0-9]+\Z")


ORTHO_COLUMNS = ("Gene stable ID", "Gene name", "Human gene stable ID", "Human gene name", "Human homology type", "Human orthology confidence [0 low, 1 high]")


ANNOTATION_COLUMNS = ("well_id", "main_structure", "sub_structure", "ontology_structure_id", "ontology_structure_acronym", "hemisphere", "brain")


class SchemaError(ValueError):
    """A schema/unit/identity violation; no implicit repair is allowed."""


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def require_hash(path: Path, expected: str) -> None:
    if not re.fullmatch(r"[0-9a-f]{64}", expected) or sha256(path) != expected:
        raise SchemaError(f"SHA256 mismatch: {path}")


def canonical_human_id(token: str) -> str | None:
    """Only a terminal dot + numeric ENSG version is removable; no case folding."""
    return token.split(".", 1)[0] if HUMAN_ID.fullmatch(token) else None


class Resolver:
    def __init__(self, author_label_ids: dict[str, set[str]], orthology_rows: Iterable[dict]):
        self.author = author_label_ids
        self.edges = defaultdict(set)
        self.reverse_edges = defaultdict(set)
        self.relations = defaultdict(set)
        self.names = defaultdict(set)
        for row in orthology_rows:
            sid, _, hid, name, relation, confidence = (str(row.get(k, "")) for k in ORTHO_COLUMNS)
            if hid and canonical_human_id(hid) != hid:
                raise SchemaError("Fixed orthology human IDs must be unversioned ENSG IDs")
            if name and hid:
                self.names[name].add(hid)
            if confidence == "1" and "ortholog" in relation and hid:
                self.edges[sid].add(hid)
                self.reverse_edges[hid].add(sid)
                self.relations[(sid, hid)].add(relation)

    def source(self, label: str) -> dict:
        result = {"source_feature": label, "source_author_id": "", "human_stable_id": "", "platform_row": "", "reason": ""}
        ids = self.author.get(label, set())
        if len(ids) != 1:
            result["reason"] = "source_label_missing_author_id" if not ids else "source_label_ambiguous_author_ids"
            return result
        sid = next(iter(ids)); result["source_author_id"] = sid
        if not MACAQUE_ID.fullmatch(sid):
            result["reason"] = "invalid_author_stable_id"; return result
        targets = self.edges.get(sid, set())
        if len(targets) != 1:
            result["reason"] = "no_high_confidence_ortholog" if not targets else "source_orthology_multitarget"
            return result
        hid = next(iter(targets)); result["human_stable_id"] = hid
        if self.relations[(sid, hid)] != {"ortholog_one2one"}:
            result["reason"] = "not_strict_one2one"; return result
        if self.reverse_edges[hid] != {sid}:
            result["reason"] = "human_orthology_multisource"
        return result

    def platform(self, panel: Sequence[str], row_ids: Sequence[str]) -> list[dict]:
        if len(set(row_ids)) != len(row_ids) or any(not x or x != x.strip() for x in row_ids):
            raise SchemaError("Blank/whitespace/duplicate platform gene rows; no aggregation")
        canonical_rows = defaultdict(list)
        symbols = defaultdict(list)
        for token in row_ids:
            stable = canonical_human_id(token)
            if stable:
                canonical_rows[stable].append(token)
            else:
                symbols[token].append(token)
        if any(len(v) > 1 for v in canonical_rows.values()):
            raise SchemaError("Version stripping creates a human stable-ID row collision")
        rows_for_human = defaultdict(set)
        for hid, tokens in canonical_rows.items():
            rows_for_human[hid].update(tokens)
        for name, tokens in symbols.items():
            humans = self.names.get(name, set())
            if len(humans) == 1:
                rows_for_human[next(iter(humans))].update(tokens)
        result = []
        for label in panel:
            item = self.source(label)
            if not item["reason"]:
                candidates = rows_for_human.get(item["human_stable_id"], set())
                if len(candidates) != 1:
                    item["reason"] = "platform_human_id_absent_or_ambiguous_symbol" if not candidates else "platform_multiple_rows_for_human_id"
                else:
                    item["platform_row"] = next(iter(candidates))
            result.append(item)
        reused = Counter(item["platform_row"] for item in result if not item["reason"])
        for item in result:
            if not item["reason"] and reused[item["platform_row"]] != 1:
                item["reason"] = "platform_row_reused_for_source_features"
                item["platform_row"] = ""
        return result


@contextmanager
def table_rows(path: Path, delimiter: str = "\t", zip_member: str | None = None):
    if zip_member:
        with zipfile.ZipFile(path) as archive:
            members = [x for x in archive.namelist() if x == zip_member or x.endswith("/" + zip_member)]
            if len(members) != 1:
                raise SchemaError(f"ZIP member not unique: {zip_member}")
            with archive.open(members[0]) as stream:
                with io.TextIOWrapper(stream, encoding="utf-8-sig", newline="") as text:
                    yield csv.reader(text, delimiter=delimiter)
    else:
        opener = gzip.open if path.suffix == ".gz" else open
        with opener(path, "rt", encoding="utf-8-sig", newline="") as text:
            yield csv.reader(text, delimiter=delimiter)


def read_axis(path: Path, *, delimiter: str, zip_member: str | None = None, sample_ids: Sequence[str] | None = None) -> tuple[list[str], list[str]]:
    with table_rows(path, delimiter, zip_member) as rows:
        samples = list(sample_ids) if sample_ids is not None else next(rows)[1:]
        if not samples or any(not s or s != s.strip() for s in samples) or len(set(samples)) != len(samples):
            raise SchemaError("Blank/duplicate/whitespace sample identifiers")
        genes = []
        for row in rows:
            if len(row) != len(samples) + 1:
                raise SchemaError("Ragged count table")
            genes.append(row[0])
    return genes, samples


def ahba_groups(metadata_rows: Sequence[dict], sample_ids: Sequence[str]) -> tuple[list[str], list[list[int]]]:
    metadata = {r["ahba_sample_uid"]: r for r in metadata_rows}
    if len(metadata) != len(metadata_rows) or set(metadata) != set(sample_ids):
        raise SchemaError("AHBA metadata/count assay identities differ")
    grouped = defaultdict(list)
    for index, sample in enumerate(sample_ids):
        row = metadata[sample]
        if not row.get("ahba_donor") or not row.get("sample_name"):
            raise SchemaError("AHBA donor/tissue missing")
        grouped[(row["ahba_donor"], row["sample_name"])].append(index)
    names, positions = [], []
    for (donor, tissue), indices in sorted(grouped.items()):
        group = [metadata[sample_ids[i]] for i in indices]
        if any(len({r.get(c, "") for r in group}) != 1 for c in ANNOTATION_COLUMNS):
            raise SchemaError("Conflicting anatomy within AHBA technical group")
        flags = [r.get("replicate_sample", "") for r in group]
        if flags.count("No") != 1 or any(f not in {"No", "Yes"} for f in flags):
            raise SchemaError("AHBA technical group requires exactly one primary")
        names.append(f"{donor}|{tissue}"); positions.append(indices)
    return names, positions


def prepare_ahba_specs(spec: dict) -> list[dict]:
    """Read only bound metadata/header; count values are read by audit_table."""
    metadata_path = Path(spec["metadata"]["path"])
    require_hash(metadata_path, spec["metadata"]["sha256"])
    with metadata_path.open(encoding="utf-8-sig", newline="") as h:
        metadata = list(csv.DictReader(h))
    result = []
    for archive in spec["archives"]:
        path = Path(archive["path"]); require_hash(path, archive["sha256"])
        with table_rows(path, ",", "SampleAnnot.csv") as rows:
            header = next(rows)
            annotations = [dict(zip(header, row)) for row in rows]
        donor = archive["donor"]
        sample_ids = [f"{donor}|{r['RNAseq_sample_name']}" for r in annotations]
        donor_meta = [r for r in metadata if r["ahba_donor"] == donor]
        by_uid = {r["ahba_sample_uid"]: r for r in donor_meta}
        if set(by_uid) != set(sample_ids):
            raise SchemaError("AHBA source SampleAnnot identities mismatch bound metadata")
        for uid, annotation in zip(sample_ids, annotations):
            if any(by_uid[uid].get(k) != annotation.get(k) for k in ["sample_name", "replicate_sample", *ANNOTATION_COLUMNS]):
                raise SchemaError("AHBA ZIP versus canonical annotation mismatch")
        result.append({**archive, "cohort_id": f"AHBA_counts:{donor}", "unit": "raw_counts", "complete_original_library_certified": True, "delimiter": ",", "zip_member": "RNAseqCounts.csv", "sample_ids": sample_ids, "metadata_rows": donor_meta, "aggregation": "AHBA_donor_tissue_technical_counts_sum", "expected_genes": 22318, "expected_profiles": 121})
    if len(metadata) != 242 or len(result) != 2:
        raise SchemaError("AHBA fixed cohort requires exactly two donors / 242 assays")
    require_hash(metadata_path, spec["metadata"]["sha256"])
    return result
