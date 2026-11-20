"""Compare the fitted fold-local and full-source Group partitions (Methods S8)."""
import csv
import os
import json
from collections import defaultdict
from pathlib import Path
from statistics import fmean


OUTPUT = Path(__file__).resolve().parent
MODELS = OUTPUT.parent / "source_analysis" / "results"
HISTORICAL = Path(os.environ.get("BRAINTRACE_REVIEW_ROOT", OUTPUT.parent.parent)) / "review_round8_20261107" / "issue_resolution"
METRICS = ("exact", "jaccard", "restricted_exact", "restricted_jaccard")


def read_reference(path):
    obj = json.loads(path.read_text(encoding="utf-8"))
    parent_and_members = {
        (item["network_id"], item["region_id"]): {
            (item["network_id"], region) for region in item["group_members"]
        }
        for item in obj["annotations"].values()
    }
    return tuple(sorted(obj["beam"])), parent_and_members


def mean_metrics(rows):
    return {metric: fmean(row[metric] for row in rows) for metric in METRICS}


def write_csv(name, rows):
    with (OUTPUT / name).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0])
        writer.writeheader()
        writer.writerows(rows)


def old_keys_and_summary(folder, date):
    summary = json.loads((folder / f"RESULTS_{date}.json").read_text(encoding="utf-8-sig"))
    keys = set()
    with (folder / f"COMMON_SUPPORT_ANCHOR_RESULTS_{date}.jsonl").open(encoding="utf-8-sig") as file:
        for line in file:
            row = json.loads(line)
            keys.add((row["fold"], tuple(sorted(row["unordered_Networks"]))))
    return keys, summary["nine_donor_equal_weight_summary"]


def main():
    final = dict(read_reference(path) for path in (MODELS / "final" / "lower").glob("*/RULES_AND_MEMBERS.json"))
    with (MODELS / "source_predictions.csv").open(encoding="utf-8-sig", newline="") as file:
        donors = {int(row["fold"]): row["donor"] for row in csv.DictReader(file)}
    objects = []
    current_keys = set()
    by_donor = defaultdict(list)
    for path in sorted((MODELS / "folds").glob("*/lower/*/RULES_AND_MEMBERS.json")):
        fold = int(path.parents[2].name.split("_")[1])
        donor = donors[fold]
        beam, local = read_reference(path)
        current_keys.add((fold, beam))
        reference = final[beam]
        shared = set(local) & set(reference)
        scores = []
        for region in sorted(shared):
            a, b = local[region], reference[region]
            aa, bb = a & shared, b & shared
            scores.append(dict(exact=a == b, jaccard=len(a & b) / len(a | b), restricted_exact=aa == bb,
                               restricted_jaccard=len(aa & bb) / len(aa | bb)))
        row = dict(fold=fold, donor=donor, networks=json.dumps(beam), common_regions=len(shared),
                   partial_support=set(local) != set(reference), **mean_metrics(scores))
        objects.append(row)
        by_donor[donor].append(row)
    donor_rows = [dict(donor=donor, objects=len(rows), regions=sum(r["common_regions"] for r in rows),
                       **mean_metrics(rows)) for donor, rows in by_donor.items()]
    old_keys, old_summary = old_keys_and_summary(
        HISTORICAL / "P0-4_Group_member_stability_descriptive_results_20261111", "20261111"
    )
    secondary_keys, secondary_summary = old_keys_and_summary(
        HISTORICAL / "P0-4_original_route_group_stability_results_20261112", "20261112"
    )
    new_summary = dict(donors=len(donor_rows), objects=len(objects), support_anchors=sum(r["regions"] for r in donor_rows),
                       partial_support=sum(r["partial_support"] for r in objects),
                       possible=len(final) * len(donor_rows), unpaired=len(final) * len(donor_rows) - len(objects),
                       **mean_metrics(donor_rows))
    result = dict(fold_local_floor_clip=new_summary, original_common_background=old_summary,
                  old_only=sorted(old_keys - current_keys), new_only=sorted(current_keys - old_keys),
                  original_secondary=secondary_summary,
                  original_secondary_overlapping_new=sorted(secondary_keys & current_keys),
                  positions_without_revised_or_original_secondary_partition=(
                      new_summary["possible"] - len(current_keys | secondary_keys)
                  ))
    write_csv("S8_GROUP_PARTITIONS_PER_OBJECT.csv", objects)
    write_csv("S8_GROUP_PARTITIONS_PER_DONOR.csv", donor_rows)
    (OUTPUT / "S8_GROUP_STABILITY.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result))


if __name__ == "__main__":
    main()
