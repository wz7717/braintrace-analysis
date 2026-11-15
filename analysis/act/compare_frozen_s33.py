"""Compare a four-column ACT replay with the original frozen S33 results."""

import argparse
import csv
import hashlib
import json
from collections import defaultdict
from fractions import Fraction
from pathlib import Path


def read_csv(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def body_sha(receipt):
    return next((attempt["sha256"] for attempt in reversed(receipt["attempts"]) if "sha256" in attempt), None)


def donor_equal(predictions, group_field):
    by_donor = defaultdict(list)
    for row in predictions:
        by_donor[row[group_field]].append(int(row["top1"]))
    donor_rates = [Fraction(sum(indicators), len(indicators)) for indicators in by_donor.values()]
    return sum(donor_rates) / len(donor_rates)


def compare_scores(original, replay, columns):
    equal = 0
    differences = []
    for field in columns:
        left, right = original[field], replay[field]
        if left and right:
            difference = abs(float(left) - float(right))
            equal += difference == 0
            differences.append(difference)
        elif left == right:
            equal += 1
    return equal, max(differences, default=0.0)


def verify(config_path, output, original):
    config = read_json(config_path)
    old_protocol_path = original / "FROZEN_PROTOCOL.json"
    old_protocol = read_json(old_protocol_path)
    old_roster_path = original / old_protocol["frozen_artifacts"]["roster"]["file"]
    minimal_path = config_path.parent / config["roster"]["file"]
    background_path = config_path.parent / config["background"]["file"]
    old_roster = read_csv(old_roster_path)
    minimal = read_csv(minimal_path)
    original_rows = {row["GEO_GSM"]: row for row in read_csv(original / "PREDICTIONS.csv")}
    old_sources = {row["GEO_GSM"]: row for row in read_json(original / "TPM_SOURCE_RECEIPTS.json")}
    new_rows = read_csv(output / "PREDICTIONS.csv")
    sources = read_json(output / "SOURCE_RECEIPTS.json")
    old_summary = read_json(original / "RESULTS.json")
    new_summary = read_json(output / "RESULTS.json")
    classes = [row["id"] for row in old_protocol["class_order"]]
    logit_fields = [f"logit_{label}" for label in classes]
    score_fields = [f"score_{label}" for label in classes]
    status_fields = ("status", "failure_reason")
    prediction_fields = ("truth_rank", "predicted_class_id", "predicted_name", "top1")
    truth_fields = ("region", "truth_class_id", "truth_name")

    source_to_group = {}
    group_to_source = {}
    checks = []
    for i, (frozen, short, new, source) in enumerate(
        zip(old_roster, minimal, new_rows, sources["tpm"], strict=True), start=1
    ):
        key = f"S{i:02d}"
        old = original_rows[frozen["GEO_GSM"]]
        old_source = old_sources[frozen["GEO_GSM"]]
        old_group = frozen["allen_donor_id"]
        group = short["donor_group"]
        one_to_one = source_to_group.setdefault(old_group, group) == group
        one_to_one &= group_to_source.setdefault(group, old_group) == old_group
        frozen_match = (
            one_to_one and source["specimen_key"] == new["specimen_key"] == key
            and old["allen_donor_id"] == old_group and new["donor_group"] == group
            and frozen["region"] == short["brain_region"] == old["region"] == new["region"]
            and all(frozen[field] == old[field] == new[field] for field in truth_fields)
            and frozen["Allen_gene_TPM_url"] == short["tpm_url"] == source["url"] == old_source["url"]
        )
        input_match = (
            source["status"] == "MATCH" and old_source["status"] == "PASS"
            and body_sha(old_source) == short["sha256"] == source["expected_sha256"] == body_sha(source)
        )
        logits_equal, logit_delta = compare_scores(old, new, logit_fields)
        scores_equal, score_delta = compare_scores(old, new, score_fields)
        ranking_match = old["ranked_class_ids"] == new["ranked_class_ids"] and all(
            old[field] == new[field] for field in prediction_fields
        )
        status_match = all(old[field] == new[field] for field in status_fields)
        checks.append({
            "specimen_key": key,
            "donor_group": group,
            "region": new["region"],
            "frozen_roster_match": frozen_match,
            "public_TPM_sha256_match": input_match,
            "logits_equal": logits_equal,
            "scores_equal": scores_equal,
            "logit_max_absolute_difference": logit_delta,
            "score_max_absolute_difference": score_delta,
            "rank_and_prediction_match": ranking_match,
            "status_and_failure_reason_match": status_match,
            "original_status": old["status"],
            "replay_status": new["status"],
        })

    old_donors = read_csv(original / "BY_DONOR.csv")
    expected_donors = {
        source_to_group[row["allen_donor_id"]]: {k: v for k, v in row.items() if k != "allen_donor_id"}
        for row in old_donors
    }
    new_donors = read_csv(output / "BY_DONOR.csv")
    replay_donors = {
        row["donor_group"]: {k: v for k, v in row.items() if k != "donor_group"}
        for row in new_donors
    }
    old_rate = donor_equal(original_rows.values(), "allen_donor_id")
    new_rate = donor_equal(new_rows, "donor_group")
    endpoint_fields = (
        "scored_specimens", "technical_failures", "model_misclassifications",
        "pooled_top1_hits", "pooled_top1", "primary_donor_equal_top1", "primary_exact_fraction", "by_region",
    )
    conditions = {
        "frozen_protocol": sha256(old_protocol_path) == config["source_protocol_sha256"],
        "original_roster_hash": sha256(old_roster_path) == old_protocol["frozen_artifacts"]["roster"]["sha256"],
        "minimal_roster_and_background_hashes": (
            sha256(minimal_path) == config["roster"]["sha256"] == new_summary["roster_sha256"]
            and sha256(background_path) == config["background"]["sha256"] == new_summary["background_sha256"]
            == old_protocol["frozen_artifacts"]["background"]["sha256"]
        ),
        "model_public_sha256": (
            sources["model"]["status"] == "MATCH"
            and body_sha(sources["model"]) == config["model"]["sha256"] == new_summary["model_sha256"]
            == old_protocol["frozen_artifacts"]["model"]["sha256"]
        ),
        "entry_and_config_hashes": (
            sha256(config_path) == new_summary["config_sha256"]
            and sha256(config_path.parent / "replay_s33.py") == new_summary["runner_sha256"]
        ),
        "same_21_specimens_and_11_donor_groups": (
            len(checks) == len(original_rows) == 21 and len(group_to_source) == len(source_to_group) == 11
            and len(old_donors) == len(new_donors) == 11
            and all(row["frozen_roster_match"] for row in checks)
        ),
        "all_TPM_bodies_match_original": all(row["public_TPM_sha256_match"] for row in checks),
        "all_logit_cells_match": all(row["logits_equal"] == len(classes) for row in checks),
        "all_softmax_score_cells_match": all(row["scores_equal"] == len(classes) for row in checks),
        "all_sort_orders_and_predictions_match": all(row["rank_and_prediction_match"] for row in checks),
        "all_statuses_and_failure_reasons_match": all(row["status_and_failure_reason_match"] for row in checks),
        "donor_rows_match_after_recoding": expected_donors == replay_donors,
        "reaggregated_donor_equal_top1_matches": (
            old_rate == new_rate
            and str(new_rate) == old_summary["primary_exact_fraction"] == new_summary["primary_exact_fraction"]
            and float(new_rate) == old_summary["primary_donor_equal_top1"] == new_summary["primary_donor_equal_top1"]
        ),
        "endpoints_match": all(old_summary[field] == new_summary[field] for field in endpoint_fields),
        "region_rows_and_confusion_counts_match": all(
            read_csv(original / name) == read_csv(output / name) for name in ("BY_REGION.csv", "CONFUSION.csv")
        ),
        "single_scoring_pass": (
            new_summary["inference_passes"] == 1
            and new_summary["classifier_fits"] == new_summary["tuning_runs"] == 0
        ),
    }
    with (output / "SPECIMEN_COMPARISON.csv").open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(checks[0]))
        writer.writeheader()
        writer.writerows(checks)
    files_to_bind = ("PREDICTIONS.csv", "BY_DONOR.csv", "BY_REGION.csv", "CONFUSION.csv", "RESULTS.json")
    receipt = {
        "status": "PASS" if all(conditions.values()) else "MISMATCH",
        "scope": "Local replay from the four-column ACT roster and the original fixed model",
        "original_directory": str(original),
        "config_sha256": sha256(config_path),
        "runtime": new_summary["runtime"],
        "verified": conditions,
        "comparison": {
            "specimens": len(checks), "donors": len(group_to_source),
            "same_TPM_sha256": sum(row["public_TPM_sha256_match"] for row in checks),
            "same_logit_values": sum(row["logits_equal"] for row in checks),
            "same_score_values": sum(row["scores_equal"] for row in checks),
            "maximum_absolute_logit_difference": max(row["logit_max_absolute_difference"] for row in checks),
            "maximum_absolute_score_difference": max(row["score_max_absolute_difference"] for row in checks),
            "same_rank_and_prediction_rows": sum(row["rank_and_prediction_match"] for row in checks),
            "same_status_and_failure_rows": sum(row["status_and_failure_reason_match"] for row in checks),
            "donor_equal_top1_original": str(old_rate),
            "donor_equal_top1_replay": str(new_rate),
            "endpoint_values_replay": {k: new_summary[k] for k in endpoint_fields},
        },
        "original_sha256": {name: sha256(original / name) for name in (
            "FROZEN_PROTOCOL.json", "FROZEN_ROSTER.csv", "TPM_SOURCE_RECEIPTS.json", *files_to_bind
        )},
        "replay_sha256": {name: sha256(output / name) for name in (
            "SOURCE_RECEIPTS.json", "SPECIMEN_COMPARISON.csv", *files_to_bind
        )},
    }
    with (output / "REPLAY_RECEIPT.json").open("x", encoding="utf-8", newline="\n") as stream:
        json.dump(receipt, stream, indent=2, ensure_ascii=False, allow_nan=False)
        stream.write("\n")
    print(receipt["status"], json.dumps(receipt["comparison"], ensure_ascii=False))
    return all(conditions.values())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True, help="Completed replay directory")
    parser.add_argument("--original", type=Path, required=True, help="Read-only original S33 directory")
    args = parser.parse_args()
    raise SystemExit(0 if verify(args.config.resolve(), args.output.resolve(), args.original.resolve()) else 1)
