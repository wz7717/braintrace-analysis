"""Compare the two saved-fold lower-ranking recipes at the same Network gate."""

from collections import Counter
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from time import perf_counter

import numpy as np

from common.pearson_hierarchy import QueryFailure, _centered, _level_means

SOURCE_ARM = "floor_clip_hybrid"
ENDPOINTS = ("Exact_Top1", "Exact_Top3", "Group_Top1", "Group_Top3")
METHODS = ("BrainTrace", "Pearson")


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8"))


def write_csv(path, rows):
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0])
        writer.writeheader()
        writer.writerows(rows)


def get_source_arm(path):
    with path.open(encoding="utf-8", newline="") as stream:
        return {
            row["sample_id"]: row
            for row in csv.DictReader(stream)
            if row["arm"] == SOURCE_ARM
        }


def build_group_map(rule, candidates, gate, train_index_set, source_regions, source_networks):
    """Recover the exact shared Group label partition from the original fold."""
    assert rule["beam"] == sorted(gate)
    assert rule["candidates"] == candidates
    groups = {}
    for region in candidates:
        ann = rule["annotations"][region]
        source_indices = rule["training_sample_indices"][region]
        assert set(source_indices) <= train_index_set
        assert all(source_regions[i] == region and source_networks[i] == ann["network_id"]
                   for i in source_indices)
        groups[ann["resolution_group"]] = sorted(ann["group_members"])
    assert sorted(region for members in groups.values() for region in members) == candidates
    assert all(region in groups[rule["annotations"][region]["resolution_group"]]
               for region in candidates)
    return groups


def original_hits(row, truth, group_map):
    exact = json.loads(row["Exact_ranking_json"])
    group = json.loads(row["Group_ranking_json"])
    groups = json.loads(row["Group_members_json"])
    assert all(sorted(members) == group_map[label]
               for label, members in zip(group, groups, strict=True))
    hits = [
        int(truth in exact[:1]),
        int(truth in exact[:3]),
        int(any(truth in members for members in groups[:1])),
        int(any(truth in members for members in groups[:3])),
    ]
    assert hits == [int(row[endpoint.replace("_Top", "_hit")]) for endpoint in ENDPOINTS]
    return dict(zip(ENDPOINTS, hits, strict=True)), exact, group


def score_queries(source_npz, folds_json, predictions_csv, reference_paths, rules_index):
    # np.load opens only the requested names. The VSD member remains unopened.
    with np.load(source_npz, allow_pickle=False) as source:
        x = source["x"]
        genes = source["genes"]
        samples = source["samples"]
        donors = source["donors"]
        regions = source["regions"]
        networks = source["labels"]

    original = get_source_arm(predictions_csv)
    rules = {
        (entry["fold"], tuple(sorted(entry["beam"]))): (entry["rule_id"], rules_index.parent / entry["path"])
        for entry in read_json(rules_index)
    }
    predictions = []
    fold_reference = []
    for fold in read_json(folds_json):
        number = fold["fold"]
        donor = fold["heldout_donor"]
        train = sorted(fold["train_indices"], key=lambda i: samples[i].encode("utf-8"))
        train_index_set = set(train)
        query = fold["test_indices"]
        assert donor not in donors[train] and set(donors[query]) == {donor}
        assert samples[query].tolist() == fold["test_sample_ids"]
        original_reference = read_json(reference_paths[str(number)])
        gene_indices = original_reference["gene_indices"]
        assert genes[gene_indices].tolist() == original_reference["G"]

        rows = [{"region": str(regions[i]), "donor": str(donors[i])} for i in train]
        labels, means, support = _level_means(x[np.ix_(gene_indices, train)], rows, "region")
        fine_hash = hashlib.sha256(np.asarray(means, dtype="<f8").tobytes(order="C")).hexdigest()
        assert fine_hash == original_reference["fine"]["centroid_sha256_le_float64_C"]
        centered_references = {label: _centered(mean) for label, mean in zip(labels, means, strict=True)}
        fold_reference.append({
            "fold": number, "heldout_donor": donor, "train_count": len(train),
            "query_count": len(query), "eligible_genes": len(gene_indices),
            "fine_reference_sha256_le_float64_C": fine_hash, "fine_support": support,
        })
        group_cache = {}
        for i in query:
            sample = str(samples[i])
            truth = str(regions[i])
            previous = original[sample]
            assert (int(previous["fold"]), previous["donor"], previous["truth_region"],
                    previous["truth_network"]) == (number, donor, truth, networks[i])
            gate = json.loads(previous["beam_json"])
            candidates = json.loads(previous["candidate_regions_json"])
            rule_key = tuple(sorted(gate))
            rule_id = None
            if candidates and rule_key not in group_cache:
                rule_id, rule_path = rules[(number, rule_key)]
                groups = build_group_map(
                    read_json(rule_path), candidates, gate, train_index_set, regions, networks
                )
                group_cache[rule_key] = rule_id, groups
            if candidates:
                rule_id, groups = group_cache[rule_key]
            else:
                groups = {}
            assert len(groups) == int(previous["n_candidate_groups"])
            assert set(candidates) <= centered_references.keys()
            original_outcomes, original_exact, original_groups = original_hits(previous, truth, groups)

            # Cached reference center/norm and the prior math.fsum scorer give the
            # same Pearson arithmetic as the 11 November own-gate baseline.
            exact_scores, group_scores = {}, {}
            failure = previous["Network_failure"] if not candidates else ""
            try:
                if candidates:
                    values, query_norm = _centered(x[gene_indices, i].astype(np.float64))
                    for label in candidates:
                        ref, ref_norm = centered_references[label]
                        dot = math.fsum(a * b for a, b in zip(values, ref, strict=True))
                        exact_scores[label] = dot / math.sqrt(query_norm * ref_norm)
                    group_scores = {
                        label: max(exact_scores[region] for region in members)
                        for label, members in groups.items()
                    }
            except QueryFailure as error:
                failure = str(error)

            pearson_exact = sorted(exact_scores, key=lambda label: (-exact_scores[label], label.encode("utf-8")))
            pearson_groups = sorted(group_scores, key=lambda label: (-group_scores[label], label.encode("utf-8")))
            pearson_outcomes = dict(zip(ENDPOINTS, [
                int(truth in pearson_exact[:1]),
                int(truth in pearson_exact[:3]),
                int(any(truth in groups[label] for label in pearson_groups[:1])),
                int(any(truth in groups[label] for label in pearson_groups[:3])),
            ], strict=True))

            in_pool = truth in candidates
            supported = truth in centered_references
            assert supported == (previous["fine_supported"] == "True")
            assert in_pool or not any(original_outcomes.values()) and not any(pearson_outcomes.values())
            status = ("TRAIN_UNSUPPORTED" if not supported else "IN_POOL" if in_pool else
                      "OUTSIDE_POOL" if candidates else "UPSTREAM_FAILURE")
            predictions.append({
                "source_col0": i, "sample_id": sample, "fold": number, "donor": donor,
                "truth_network": str(networks[i]), "truth_region": truth,
                "pool_status": status, "network_top3": gate,
                "candidate_regions": candidates,
                "group_rule": rule_id,
                "groups": groups, "truth_group": next((label for label, members in groups.items()
                                                      if truth in members), None),
                "BrainTrace_Exact_rank": original_exact, "BrainTrace_Group_rank": original_groups,
                "BrainTrace_Exact_failure": previous["Exact_failure"],
                "BrainTrace_Group_failure": previous["Group_failure"],
                "Pearson_failure": failure,
                "Pearson_Exact_rank_scores": [(label, exact_scores[label]) for label in pearson_exact],
                "Pearson_Group_rank_scores": [(label, group_scores[label]) for label in pearson_groups],
                "BrainTrace": original_outcomes, "Pearson": pearson_outcomes,
            })
        print(f"Fold {number} ({donor}), {len(query)} matched queries", flush=True)
    predictions.sort(key=lambda row: row["source_col0"])
    assert len(predictions) == len(original) == len(samples)
    return predictions, fold_reference


def summarize(predictions):
    summary = []
    per_donor = []
    for subset, records in (
        ("all_819", predictions),
        ("shared_in_pool", [row for row in predictions if row["pool_status"] == "IN_POOL"]),
    ):
        donors = sorted({row["donor"] for row in records})
        for endpoint in ENDPOINTS:
            cells = []
            for donor in donors:
                queries = [row for row in records if row["donor"] == donor]
                cell = {
                    "subset": subset, "endpoint": endpoint, "donor": donor, "n": len(queries),
                    **{method + "_hits": sum(row[method][endpoint] for row in queries) for method in METHODS},
                    "BrainTrace_only": sum(row["BrainTrace"][endpoint] > row["Pearson"][endpoint] for row in queries),
                    "Pearson_only": sum(row["Pearson"][endpoint] > row["BrainTrace"][endpoint] for row in queries),
                }
                for method in METHODS:
                    cell[method + "_proportion"] = cell[method + "_hits"] / len(queries)
                cells.append(cell)
            per_donor.extend(cells)
            aggregate = {"subset": subset, "endpoint": endpoint, "n": len(records), "donors": len(donors)}
            for method in METHODS:
                aggregate[method + "_hits"] = sum(cell[method + "_hits"] for cell in cells)
                aggregate[method + "_donor_equal"] = math.fsum(
                    cell[method + "_proportion"] for cell in cells
                ) / len(donors)
            aggregate["BrainTrace_minus_Pearson_pp"] = (
                aggregate["BrainTrace_donor_equal"] - aggregate["Pearson_donor_equal"]
            ) * 100
            aggregate["BrainTrace_only"] = sum(cell["BrainTrace_only"] for cell in cells)
            aggregate["Pearson_only"] = sum(cell["Pearson_only"] for cell in cells)
            summary.append(aggregate)
    return summary, per_donor


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    cfg = read_json(args.config)
    directory = args.config.resolve().parent
    references = {str(fold): directory / file for fold, file in cfg["frozen_references"].items()}
    output = args.out.resolve()
    output.mkdir(parents=True, exist_ok=True)
    start = perf_counter()
    predictions, fold_reference = score_queries(
        directory / cfg["source_npz"], directory / cfg["folds_json"],
        directory / cfg["predictions_csv"], references, directory / cfg["group_rules_index"],
    )
    summary, per_donor = summarize(predictions)
    with (output / "MATCHED_PREDICTIONS.jsonl").open("w", encoding="utf-8", newline="\n") as stream:
        for row in predictions:
            stream.write(json.dumps(row, ensure_ascii=False, allow_nan=False) + "\n")
    write_csv(output / "MATCHED_ENDPOINTS.csv", summary)
    write_csv(output / "BY_DONOR.csv", per_donor)
    run = {
        "status": "COMPLETE",
        "source_arm": SOURCE_ARM,
        "queries": len(predictions),
        "donors": len(fold_reference),
        "actual_fold_beams": len({row["group_rule"] for row in predictions if row["group_rule"]}),
        "pool_status": Counter(row["pool_status"] for row in predictions),
        "failures": {field: Counter(row[field] for row in predictions if row[field])
                     for field in ("BrainTrace_Exact_failure", "BrainTrace_Group_failure", "Pearson_failure")},
        "pearson_exact_tie_queries": sum(
            len({score for _, score in row["Pearson_Exact_rank_scores"]}) < len(row["Pearson_Exact_rank_scores"])
            for row in predictions
        ),
        "pearson_group_tie_queries": sum(
            len({score for _, score in row["Pearson_Group_rank_scores"]}) < len(row["Pearson_Group_rank_scores"])
            for row in predictions
        ),
        "fold_reference": fold_reference,
        "elapsed_seconds": perf_counter() - start,
    }
    (output / "RUN_SUMMARY.json").write_text(json.dumps(run, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({key: value for key, value in run.items() if key != "fold_reference"}), flush=True)


if __name__ == "__main__":
    main()
