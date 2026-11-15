"""Compare donor-equal Network Pearson with archived floor_clip on the original AHBA task."""

from collections import defaultdict
import json
from statistics import mean
from time import monotonic

import numpy as np

from . import ahba_inputs as inputs
from ._config import arguments, location
from .shared_head_external import read_json, read_csv, write_csv, field, fail_unless, rank_scores

METRICS = ("top1_hit", "top3_hit", "MRR", "precision_at_3", "recall_at_3")


def reference(state, parameter_file, source_archive):
    with np.load(parameter_file, allow_pickle=False) as bundle, np.load(source_archive, allow_pickle=False) as source:
        rows = bundle["feature_indices"]
        panel = bundle["genes"].tolist()
        for source_key, bundle_key in (
            ("samples", "train_sample_ids"), ("labels", "train_labels"), ("donors", "train_donors")
        ):
            fail_unless(np.array_equal(source[source_key], bundle[bundle_key]),
                        f"Historical source specimen axis differs: {source_key}")
        fail_unless(source["genes"][rows].tolist() == panel, "Historical 5000 panel differs from source gene coordinates")
        _, train = inputs.transform(source["x"][rows].T.copy(), np.ones(len(panel), dtype=bool), state)
        classes = bundle["classes"].tolist()
        labels = bundle["train_labels"]
        donors = bundle["train_donors"]
        centroids = np.array([
            np.mean([train[(labels == label) & (donors == donor)].mean(axis=0, dtype=np.float64)
                     for donor in sorted(set(donors[labels == label]))], axis=0)
            for label in classes
        ])
        coef = bundle["coef"]
        intercept = bundle["intercept"]
    return panel, classes, centroids, coef, intercept


def score_queries(state, panel, classes, centroids, coef, intercept, cfg, base):
    frozen = read_json(location(cfg, base, "frozen_input"))
    counts = frozen["ahba_counts"]
    counts["metadata"]["path"] = str(location(cfg, base, "metadata"))
    for archive in counts["archives"]:
        archive["path"] = str(location(cfg["archives"], base, archive["donor"]))
    resolver = inputs.load_resolver(location(cfg, base, "author_map"), location(cfg, base, "orthology"))
    metadata = read_json(location(cfg, base, "feature_coverage"))["routes_by_donor"]
    archived = {r["profile_id"]: r for r in read_json(location(cfg, base, "old_predictions"))["routes"]["floor_clip"]["rows"]}
    result, donors, max_logit_error = {}, {}, 0.0
    for spec in inputs.prepare_ahba_specs(counts):
        axis, _ = inputs.read_axis(
            spec["path"], delimiter=",", zip_member="RNAseqCounts.csv", sample_ids=spec["sample_ids"]
        )
        donor = spec["donor"]
        mapping = resolver.platform(panel, axis)
        observed = np.array([not r["reason"] for r in mapping], dtype=bool)
        fail_unless(observed.tolist() == metadata["floor_clip|" + donor]["available"],
                    f"Historical feature-availability axis differs: {donor}")
        profiles, x, totals, reasons = inputs._stream_observed_counts(spec, mapping)
        usable = np.array([
            not reason and observed.sum() >= 2 and np.ptp(row[observed].astype(float)) > 1e-12
            for row, reason in zip(x, reasons, strict=True)
        ])
        _, prepared = inputs.transform(x[usable], observed, state)
        logits = prepared @ coef.T + intercept
        stored_logits = np.asarray([archived[profile]["decision"]
                                    for profile in np.asarray(profiles)[usable]])
        max_logit_error = max(max_logit_error, float(np.max(np.abs(logits - stored_logits))))
        fail_unless(np.allclose(logits, stored_logits, atol=1e-10, rtol=1e-12),
                    f"Restaged target logits differ from frozen floor_clip: {donor}")
        for i, profile in enumerate(profiles):
            fail_unless(float(totals[i]) == archived[profile]["full_library_total"],
                        f"Full-library count denominator differs: {profile}")
            result[profile] = dict(ranking=[], failure="" if usable[i] else (reasons[i] or "insufficient_observed_signal"),
                                   full_library_total=float(totals[i]))
        for profile, row in zip(np.asarray(profiles)[usable], prepared, strict=True):
            result[profile]["ranking"] = rank_scores(row.astype(np.float64), classes, centroids)
        donors[donor] = dict(n=len(profiles), observed_panel=int(observed.sum()),
                             missing_panel=int((~observed).sum()), query_failures=int((~usable).sum()))
    fail_unless(set(result) == set(archived), "Historical technical roster differs")
    return result, donors, max_logit_error


def metric_values(ranking, truth):
    if not truth:
        return {key: "" for key in METRICS}
    rank = next((i for i, label in enumerate(ranking, 1) if label in truth), 0)
    overlap = len(truth.intersection(ranking[:3]))
    return dict(top1_hit=int(rank == 1), top3_hit=int(0 < rank <= 3), MRR=1 / rank if rank else 0,
                precision_at_3=overlap / 3, recall_at_3=overlap / len(truth))


def rows_for_evaluation(predictions, endpoint_file):
    evaluation = []
    for r in read_csv(endpoint_file):
        if r["route_id"] != "floor_clip":
            continue
        profile = r["profile_id"]
        truth = set(json.loads(r["allowed_truth_json"]))
        frozen_ranking = json.loads(r["ranking_json"])
        new = predictions[profile]
        rank = [entry["label"] for entry in new["ranking"]]
        common = dict(profile_id=profile, donor_id=r["donor_id"], allowed_truth_json=r["allowed_truth_json"],
                      truth_cardinality=len(truth), full_library_total=new["full_library_total"])
        old_metrics = metric_values(frozen_ranking, truth)
        fail_unless(not truth or all(float(r[key]) == old_metrics[key] for key in
                    ("top1_hit", "top3_hit", "precision_at_3", "recall_at_3")),
                    f"Original evaluation rule mismatch: {profile}")
        evaluation.append(dict(
            **common, method="archived_floor_clip", ranking_json=field(frozen_ranking),
            scores_json="", failure=r["prediction_reason"], **old_metrics
        ))
        evaluation.append(dict(
            **common, method="common_input_pearson", ranking_json=field(rank),
            scores_json=field([entry["score"] for entry in new["ranking"]]), failure=new["failure"],
            **metric_values(rank, truth)
        ))
    return evaluation


def summary(rows):
    summary_rows, contrasts = [], []
    subsets = {
        "allowed": lambda n: n > 0, "singleton": lambda n: n == 1, "multi": lambda n: n > 1
    }
    methods = ("archived_floor_clip", "common_input_pearson")
    for subset, eligible in subsets.items():
        by_method = {
            method: [r for r in rows if r["method"] == method and eligible(r["truth_cardinality"])]
            for method in methods
        }
        for method, selected in by_method.items():
            by_donor = defaultdict(list)
            for row in selected:
                by_donor[row["donor_id"]].append(row)
            summary_rows.append(dict(
                method=method, subset=subset, donor="pooled", n=len(selected), failures=sum(bool(r["failure"]) for r in selected),
                top1_hits=sum(r["top1_hit"] for r in selected), top3_hits=sum(r["top3_hit"] for r in selected),
                **{key: mean(r[key] for r in selected) for key in METRICS}
            ))
            for donor, group in sorted(by_donor.items()):
                summary_rows.append(dict(
                    method=method, subset=subset, donor=donor, n=len(group), failures=sum(bool(r["failure"]) for r in group),
                    top1_hits=sum(r["top1_hit"] for r in group), top3_hits=sum(r["top3_hit"] for r in group),
                    **{key: mean(r[key] for r in group) for key in METRICS}
                ))
        old = {r["profile_id"]: r for r in by_method[methods[0]]}
        for metric in METRICS:
            diffs = [r[metric] - old[r["profile_id"]][metric] for r in by_method[methods[1]]]
            contrasts.append(dict(subset=subset, metric=metric, n=len(diffs), pooled_pearson_minus_floor_clip=mean(diffs),
                                  pearson_higher=sum(x > 0 for x in diffs), pearson_lower=sum(x < 0 for x in diffs),
                                  same=sum(x == 0 for x in diffs)))
    return summary_rows, contrasts


def main():
    start = monotonic()
    cfg, base, output = arguments(__doc__)
    source_archive = location(cfg, base, "source_archive")
    parameter_file = location(cfg, base, "parameters")
    state_file = location(cfg, base, "state")
    with np.load(state_file, allow_pickle=False) as state:
        panel, classes, centroids, coef, intercept = reference(state, parameter_file, source_archive)
        scores, donors, max_logit_error = score_queries(state, panel, classes, centroids, coef, intercept, cfg, base)
    rows = rows_for_evaluation(scores, location(cfg, base, "endpoint_rows"))
    summaries, contrasts = summary(rows)
    write_csv(output / "AHBA_PER_SPECIMEN.csv", rows)
    write_csv(output / "AHBA_POOLED_AND_DONOR_SUMMARY.csv", summaries)
    write_csv(output / "AHBA_PAIRED_DIFFERENCES.csv", contrasts)
    metadata = dict(
        route="AHBA count based final floor_clip development comparison",
        input_config=cfg, frozen_source=str(source_archive),
        train_profiles=819, train_donors=9, target_profiles=len(scores), target_donors=donors,
        reference_axis=len(panel), classes=classes,
        restaged_target_max_abs_floor_clip_logit_error=max_logit_error,
        elapsed_seconds=monotonic()-start
    )
    (output / "AHBA_RUN.json").write_text(json.dumps(metadata, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"AHBA: {len(scores)} technical profiles / {sum(r['truth_cardinality'] > 0 for r in rows) // 2} original truth-eligible")


if __name__ == "__main__":
    main()
