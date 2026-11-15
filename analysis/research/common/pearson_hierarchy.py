"""Frozen donor-equal Pearson Network-to-Exact arithmetic.

The caller supplies fold training rows and truth-free query identities. This
module opens no files and uses no held-out labels to select or score references.
"""

from collections import defaultdict
from dataclasses import dataclass
import hashlib
import math

import numpy as np


class StructuralAbort(ValueError):
    """A fold cannot satisfy the frozen reference or hierarchy contract."""


class QueryFailure(ValueError):
    """One executed query has no valid Pearson score vector."""


def _canonical(values):
    return sorted(values, key=lambda value: value.encode("utf-8"))


def _centered(values):
    if not all(math.isfinite(value) for value in values):
        raise QueryFailure("QUERY_NONFINITE_VALUES")
    try:
        mean = math.fsum(values) / len(values)
        centered = tuple(value - mean for value in values)
        norm = math.fsum(value * value for value in centered)
    except (OverflowError, ValueError) as error:
        raise QueryFailure("QUERY_UNDEFINED_REDUCTION") from error
    if not math.isfinite(norm) or norm <= 0:
        raise QueryFailure("QUERY_ZERO_OR_INVALID_NORM")
    return centered, norm


def rank_scores(values, labels, means):
    """Return raw Pearson scores in descending order with exact UTF-8 ties."""
    query, query_norm = _centered(values)
    scores = []
    for label, reference in zip(labels, means, strict=True):
        centered, reference_norm = _centered(reference)
        try:
            denominator_squared = query_norm * reference_norm
            numerator = math.fsum(a * b for a, b in zip(query, centered, strict=True))
            if not math.isfinite(denominator_squared) or denominator_squared <= 0:
                raise QueryFailure("QUERY_UNDEFINED_SCORE")
            score = numerator / math.sqrt(denominator_squared)
        except (OverflowError, ValueError) as error:
            raise QueryFailure("QUERY_UNDEFINED_SCORE") from error
        if not math.isfinite(score):
            raise QueryFailure("QUERY_UNDEFINED_SCORE")
        scores.append({"label": label, "score": score})
    return sorted(scores, key=lambda item: (-item["score"], item["label"].encode("utf-8")))


def _level_means(matrix, train_rows, field):
    cells = defaultdict(lambda: defaultdict(list))
    for column, row in enumerate(train_rows):
        cells[row[field]][row["donor"]].append(column)
    labels = _canonical(cells)
    means = []
    support = []
    for label in labels:
        donors = _canonical(cells[label])
        profiles = []
        for donor in donors:
            columns = cells[label][donor]
            profiles.append([
                math.fsum(float(values[column]) for column in columns) / len(columns)
                for values in matrix
            ])
        means.append([
            math.fsum(profile[gene] for profile in profiles) / len(profiles)
            for gene in range(matrix.shape[0])
        ])
        support.append({
            "label": label,
            "train_specimens": sum(len(cells[label][donor]) for donor in donors),
            "train_donors": len(donors),
            "contributing_donors": donors,
        })
    result = np.asarray(means, dtype=np.float64)
    for label, values in zip(labels, result, strict=True):
        try:
            _centered(values)
        except QueryFailure as error:
            raise StructuralAbort(f"INVALID_TRAINING_REFERENCE:{field}:{label}:{error}") from error
    return labels, result, support


@dataclass
class Reference:
    gene_indices: list
    gene_ids: list
    network_labels: list
    fine_labels: list
    network_means: np.ndarray
    fine_means: np.ndarray
    network_support: list
    fine_support: list

    def metadata(self):
        """Compact reference identities, without exporting expression matrices."""
        result = {"G": self.gene_ids, "gene_indices": self.gene_indices}
        for level in ("network", "fine"):
            means = getattr(self, level + "_means")
            result[level] = {
                "labels": getattr(self, level + "_labels"),
                "support": getattr(self, level + "_support"),
                "centroid_shape": list(means.shape),
                "centroid_sha256_le_float64_C": hashlib.sha256(
                    np.asarray(means, dtype="<f8").tobytes(order="C")
                ).hexdigest(),
            }
        return result


def build_reference(x, genes, train_rows, parents):
    """Build each level from this fold's training specimens and donors only."""
    if len(genes) != len(set(genes)) or x.ndim != 2 or x.shape[0] != len(genes):
        raise StructuralAbort("GENE_AXIS_MISMATCH")
    rows = sorted(train_rows, key=lambda row: row["sample_id"].encode("utf-8"))
    if not rows or len({row["sample_id"] for row in rows}) != len(rows):
        raise StructuralAbort("INVALID_TRAIN_ROSTER")
    if set(row["truth"] for row in rows) != set(parents.values()):
        raise StructuralAbort("MISSING_NETWORK_TRAINING_SUPPORT")
    if any(parents.get(row["region"]) != row["truth"] for row in rows):
        raise StructuralAbort("FINE_PARENT_MISMATCH")
    matrix = np.asarray(x[:, [row["column"] for row in rows]], dtype=np.float64)
    eligible = np.isfinite(matrix).all(axis=1) & (matrix.max(axis=1) > matrix.min(axis=1))
    indices = sorted(np.flatnonzero(eligible).tolist(), key=lambda index: genes[index].encode("utf-8"))
    if len(indices) < 2:
        raise StructuralAbort("INSUFFICIENT_G")
    selected = matrix[indices]
    try:
        networks, network_means, network_support = _level_means(selected, rows, "truth")
        fine, fine_means, fine_support = _level_means(selected, rows, "region")
    except (OverflowError, ValueError) as error:
        raise StructuralAbort(f"INVALID_TRAINING_REDUCTION:{error}") from error
    return Reference(indices, [genes[index] for index in indices], networks, fine,
                     network_means, fine_means, network_support, fine_support)


def score_query(reference, full_gene_values, identity, parents):
    """Route and rank one query using its own Network Top3, without truth."""
    result = {key: identity[key] for key in ("sample_id", "column", "donor", "fold")}
    result.update({
        "Network_ranked_label_scores": None, "Network_Top3": None,
        "C": None, "K": None, "Exact_ranked_label_scores": None,
        "Network_status": "QUERY_FAILURE", "Exact_status": "UPSTREAM_FAILURE",
        "failure_reason": None,
    })
    values = np.asarray(full_gene_values, dtype=np.float64)[reference.gene_indices]
    try:
        network_scores = rank_scores(values, reference.network_labels, reference.network_means)
    except QueryFailure as error:
        result["failure_reason"] = str(error)
        return result
    result["Network_ranked_label_scores"] = network_scores
    result["Network_Top3"] = [entry["label"] for entry in network_scores[:3]]
    result["Network_status"] = "OK"
    pool_indices = [index for index, fine in enumerate(reference.fine_labels)
                    if parents[fine] in result["Network_Top3"]]
    pool = [reference.fine_labels[index] for index in pool_indices]
    if len(result["Network_Top3"]) != 3 or len(pool) < 3:
        raise StructuralAbort("INVALID_OWN_TOP3_FINE_POOL")
    result["C"], result["K"] = pool, len(pool)
    try:
        result["Exact_ranked_label_scores"] = rank_scores(
            values, pool, reference.fine_means[pool_indices]
        )
        result["Exact_status"] = "OK"
    except QueryFailure as error:
        result["Exact_status"] = "QUERY_FAILURE"
        result["failure_reason"] = str(error)
    return result
