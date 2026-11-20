"""BrainTrace numerical helpers used by the Ensembl 97 reproduction route."""
from __future__ import annotations
from __future__ import annotations
from collections import Counter
from typing import Any
import numpy as np
import pandas as pd
from scripts.run_bo2023_loso_validation import correlation_scores

class UnionFind:
    def __init__(self, items: list[str]) -> None:
        self.parent = {item: item for item in items}
        self.members = {item: {item} for item in items}

    def find(self, item: str) -> str:
        root = item
        while self.parent[root] != root:
            root = self.parent[root]
        while self.parent[item] != item:
            parent = self.parent[item]
            self.parent[item] = root
            item = parent
        return root

    def union(self, left: str, right: str, max_size: int) -> bool:
        a = self.find(left)
        b = self.find(right)
        if a == b:
            return True
        if len(self.members[a]) + len(self.members[b]) > max_size:
            return False
        if len(self.members[a]) < len(self.members[b]):
            a, b = b, a
        self.parent[b] = a
        self.members[a].update(self.members[b])
        del self.members[b]
        return True


def centroid_reference(
    values: np.ndarray,
    candidates: list[str],
    training: dict[str, np.ndarray],
    rows: np.ndarray,
) -> np.ndarray:
    return np.column_stack(
        [values[rows[:, None], training[region]].mean(axis=1, dtype=np.float64) for region in candidates]
    )


def build_resolution_groups(
    values: np.ndarray,
    candidates: list[str],
    training: dict[str, np.ndarray],
    region_network: dict[str, str | None],
    rows: np.ndarray,
    min_resolution_samples: int,
    min_merge_samples: int,
    min_pair_errors: int,
    min_confusion_rate: float,
    similarity_threshold: float,
    merge_similarity_threshold: float,
    max_group_size: int,
) -> tuple[dict[str, dict[str, Any]], pd.DataFrame]:
    annotations: dict[str, dict[str, Any]] = {}
    audit_rows: list[dict[str, Any]] = []
    reference = centroid_reference(values, candidates, training, rows)
    positions = {region: j for j, region in enumerate(candidates)}

    for region in candidates:
        annotations[region] = {
            "region_id": region,
            "network_id": region_network.get(region),
            "resolution_group": region,
            "group_members": [region],
            "resolution_tier": "high_resolution",
            "resolution_reasons": [],
            "training_samples": int(len(training[region])),
            "nearest_centroid_corr": None,
        }
        if len(training[region]) < min_resolution_samples:
            annotations[region]["resolution_reasons"].append(f"training_n<{min_resolution_samples}")
        if region_network.get(region) is None:
            annotations[region]["resolution_reasons"].append("cross_network_region")

    network_values = sorted({network for network in region_network.values() if network is not None})
    for network in network_values:
        network_regions = sorted([region for region in candidates if region_network.get(region) == network])
        if len(network_regions) < 2:
            continue
        network_positions = [positions[region] for region in network_regions]
        network_ref = reference[:, network_positions]
        for j, region in enumerate(network_regions):
            competitors = [i for i in range(len(network_regions)) if i != j]
            nearest = float(np.max(correlation_scores(network_ref[:, competitors], network_ref[:, j])))
            annotations[region]["nearest_centroid_corr"] = nearest
            if nearest >= similarity_threshold:
                annotations[region]["resolution_reasons"].append(
                    f"nearest_centroid_corr>={similarity_threshold:.2f}"
                )

        confusions: Counter[tuple[str, str]] = Counter()
        for truth in network_regions:
            truth_indices = training[truth]
            if len(truth_indices) <= 1:
                continue
            truth_pos = network_regions.index(truth)
            for sample_idx in truth_indices:
                local_ref = network_ref.copy()
                local_ref[:, truth_pos] = values[rows[:, None], truth_indices[truth_indices != sample_idx]].mean(
                    axis=1, dtype=np.float64
                )
                scores = correlation_scores(local_ref, values[rows, sample_idx])
                predicted = network_regions[int(np.argmax(scores))]
                if predicted != truth:
                    confusions[(truth, predicted)] += 1

        pair_candidates: list[dict[str, Any]] = []
        for left_i, left in enumerate(network_regions):
            for right in network_regions[left_i + 1 :]:
                errors = int(confusions[(left, right)] + confusions[(right, left)])
                pair_n = int(len(training[left]) + len(training[right]))
                rate = float(errors / pair_n) if pair_n else 0.0
                corr = float(correlation_scores(
                    reference[:, [positions[left], positions[right]]], reference[:, positions[left]]
                )[1])
                eligible = (
                    len(training[left]) >= min_merge_samples
                    and len(training[right]) >= min_merge_samples
                    and errors >= min_pair_errors
                    and rate >= min_confusion_rate
                    and corr >= merge_similarity_threshold
                )
                pair_candidates.append(
                    {
                        "network_id": network,
                        "left_region": left,
                        "right_region": right,
                        "pair_errors": errors,
                        "confusion_rate": rate,
                        "centroid_correlation": corr,
                        "eligible_for_merge": bool(eligible),
                    }
                )
        uf = UnionFind(network_regions)
        for pair in sorted(
            pair_candidates,
            key=lambda x: (-int(x["eligible_for_merge"]), -x["pair_errors"], -x["centroid_correlation"]),
        ):
            merged = False
            if pair["eligible_for_merge"]:
                merged = uf.union(str(pair["left_region"]), str(pair["right_region"]), max_group_size)
            pair["merged"] = bool(merged and pair["eligible_for_merge"])
            audit_rows.append(pair)

        groups: dict[str, list[str]] = {}
        for region in network_regions:
            groups.setdefault(uf.find(region), []).append(region)
        for members in groups.values():
            members = sorted(members)
            if len(members) <= 1:
                continue
            group_id = f"{network}::{' + '.join(members)}"
            for region in members:
                annotations[region]["resolution_group"] = group_id
                annotations[region]["group_members"] = members
                annotations[region]["resolution_reasons"].append("fold_local_merged_confusion_group")

    for annotation in annotations.values():
        reasons = list(dict.fromkeys(annotation["resolution_reasons"]))
        annotation["resolution_reasons"] = reasons
        annotation["resolution_tier"] = "low_resolution" if reasons else "high_resolution"
    return annotations, pd.DataFrame(audit_rows)
