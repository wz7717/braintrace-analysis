"""Shared-head external lower-tier comparison using the original anatomy rosters."""

from collections import defaultdict
import csv
import itertools
import json
import platform
from statistics import mean
from time import monotonic

import numpy as np

from common.pearson_hierarchy import QueryFailure, build_reference, rank_scores

from ._config import arguments, location

METHODS = ("same_head_original_lower", "same_head_new_lower", "same_head_full_G_pearson")
COHORTS = ("GSE128537", "GSE179328")
ENDPOINTS = ("primary_unique", "secondary_ambiguous", "coarsened_parent_only")
METRICS = ("Exact_hit1", "Exact_hit3", "Exact_MRR", "Group_hit1", "Group_hit3", "Group_MRR",
           "group_top1_member_count", "group_top3_union_member_count")


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def read_csv(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def write_csv(path, rows):
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def field(value):
    return json.dumps(value, ensure_ascii=False, separators=(",", ":"))


def fail_unless(condition, message):
    if not condition:
        raise ValueError(message)


def source_reference(source_archive, membership_file, region_parent_file):
    membership = read_json(membership_file)
    parents = read_json(region_parent_file)
    with np.load(source_archive, allow_pickle=False) as source:
        for array_key, member_key in (
            ("samples", "samples"), ("donors", "donors"),
            ("regions", "regions"), ("labels", "networks"),
        ):
            fail_unless(source[array_key].tolist() == membership[member_key],
                        f"B5 source and historical training disagree: {array_key}")
        train = [
            dict(sample_id=sample, donor=donor, region=region, truth=network, column=col)
            for col, (sample, donor, region, network) in enumerate(zip(
                membership["samples"], membership["donors"], membership["regions"],
                membership["networks"], strict=True
            ))
        ]
        genes = source["genes"].tolist()
        reference = build_reference(source["x"], genes, train, parents)
    return reference, genes, parents, train


def frozen_beams(reference, parents, training, definitions):
    columns = defaultdict(list)
    for row in training:
        columns[row["region"]].append(row["column"])
    index = {label: i for i, label in enumerate(reference.fine_labels)}
    beams = {}
    for record, folder in definitions:
        definition = read_json(folder / "RULES_AND_MEMBERS.json")
        beam = tuple(sorted(definition["beam"]))
        candidates = definition["candidates"]
        expected = sorted(label for label in index if parents[label] in beam)
        fail_unless(candidates == expected and list(beam) == sorted(record["networks"]),
                    f"Source candidate mismatch: {record['id']}")
        fail_unless(definition["training_sample_indices"] == {label: columns[label] for label in candidates},
                    f"Beam training specimens differ from B5: {record['id']}")
        groups = defaultdict(list)
        for label in candidates:
            groups[definition["annotations"][label]["resolution_group"]].append(label)
        fail_unless(all(
            sorted(definition["annotations"][label]["group_members"]) == members
            for members in groups.values() for label in members
        ), f"Frozen Group partition mismatch: {record['id']}")
        beams[beam] = (record["id"], candidates, dict(groups), [index[label] for label in candidates])
    return beams


def query_scores(x, observed, samples, query_ids, heads, beams, reference):
    """Rank query IDs and frozen head routes before joining anatomy or donors."""
    columns = {sample: i for i, sample in enumerate(samples)}
    observed_g = observed[reference.gene_indices]
    results = {}
    for sample_id in query_ids:
        head = json.loads(heads[sample_id]["ranking_json"])
        beam_id, candidates, groups, reference_indices = beams[tuple(sorted(head[:3]))]
        values = x[reference.gene_indices, columns[sample_id]].astype(np.float64)
        measured = values[observed_g]
        result = dict(
            beam_id=beam_id, head=head, candidates=candidates, groups=groups,
            exact=[], group_labels=[], group_members=[], group_scores=[],
            Exact_failure="", Group_failure=""
        )
        try:
            if measured.size < 2:
                raise QueryFailure("QUERY_OBSERVED_COUNT_LT_2")
            if not np.isfinite(measured).all() or measured.max() <= measured.min():
                raise QueryFailure("QUERY_OBSERVED_NONFINITE_OR_CONSTANT")
            ranked = rank_scores(values, candidates, reference.fine_means[reference_indices])
        except QueryFailure as error:
            result["Exact_failure"] = result["Group_failure"] = str(error)
        else:
            by_label = {entry["label"]: entry["score"] for entry in ranked}
            group_scores = {label: max(by_label[member] for member in members)
                            for label, members in groups.items()}
            group_labels = sorted(groups, key=lambda label: (-group_scores[label], label.encode("utf-8")))
            result.update(
                exact=ranked, group_labels=group_labels,
                group_members=[groups[label] for label in group_labels],
                group_scores=[group_scores[label] for label in group_labels]
            )
        results[sample_id] = result
    return results


def first_hit(ranked_member_sets, allowed):
    return next((rank for rank, members in enumerate(ranked_member_sets, 1)
                 if allowed.intersection(members)), 0)


def endpoints(exact, group_members, allowed):
    out = {}
    for tier, ranked in (("Exact", [[label] for label in exact]), ("Group", group_members)):
        rank = first_hit(ranked, allowed)
        out.update({tier + "_allowed_rank": rank, tier + "_hit1": int(rank == 1),
                    tier + "_hit3": int(0 < rank <= 3), tier + "_MRR": 1 / rank if rank else 0.0})
    out["group_top1_member_count"] = len(group_members[0]) if group_members else 0
    out["group_top3_union_member_count"] = len(set().union(*group_members[:3]))
    return out


def comparable_rows(roster, historic, scores):
    rows = []
    for truth in roster:
        sample = truth["sample_id"]
        score = scores[sample]
        allowed = set(json.loads(truth["allowed_regions_json"]))
        common = {key: truth[key] for key in (
            "cohort", "sample_id", "GSM", "donor", "region_original", "fine_endpoint",
            "allowed_regions_json", "allowed_set_size"
        )}
        common.update(beam_id=score["beam_id"], network_top3_json=field(score["head"][:3]),
                      candidates_json=field(score["candidates"]), groups_json=field(score["groups"]),
                      allowed_in_pool=int(bool(allowed.intersection(score["candidates"]))))
        for method in METHODS[:2]:
            old = historic[(sample, method)]
            sample_msg = f"{truth['cohort']}:{sample}:{method}"
            fail_unless(all(old[key] == common[key] for key in (
                "GSM", "donor", "region_original", "fine_endpoint", "allowed_regions_json",
            )), f"Saved prediction truth binding differs: {sample_msg}")
            fail_unless(json.loads(old["Network_ranking_json"]) == score["head"] and
                        json.loads(old["candidate_regions_json"]) == score["candidates"],
                        f"Saved prediction candidates differ: {sample_msg}")
            exact = json.loads(old["Exact_ranking_json"]) if old["Exact_valid"] == "True" else []
            groups = json.loads(old["Group_members_json"]) if old["Group_valid"] == "True" else []
            group_labels = json.loads(old["Group_ranking_json"]) if groups else []
            fail_unless((not exact or sorted(exact) == score["candidates"]) and
                        (not groups or {key: sorted(vals) for key, vals in zip(group_labels, groups, strict=True)}
                         == score["groups"]), f"Saved ranking or Group partition differs: {sample_msg}")
            metrics = endpoints(exact, groups, allowed)
            fail_unless(all(float(old[key]) == metrics[key] for key in METRICS),
                        f"Saved metric and bound truth differ: {sample_msg}")
            rows.append(dict(
                **common, method=method, inference_status=old["inference_status"],
                Exact_failure=old["Exact_failure"], Group_failure=old["Group_failure"],
                exact_ranking_json=field(exact), exact_scores_json="",
                group_ranking_json=field(group_labels), group_ranking_members_json=field(groups),
                group_ranking_scores_json=old["Group_scores_json"], **metrics
            ))
        exact = [entry["label"] for entry in score["exact"]]
        rows.append(dict(
            **common, method=METHODS[2], inference_status="predicted",
            Exact_failure=score["Exact_failure"], Group_failure=score["Group_failure"],
            exact_ranking_json=field(exact), exact_scores_json=field([entry["score"] for entry in score["exact"]]),
            group_ranking_json=field(score["group_labels"]),
            group_ranking_members_json=field(score["group_members"]),
            group_ranking_scores_json=field(score["group_scores"]),
            **endpoints(exact, score["group_members"], allowed)
        ))
    return rows


def summaries(rows, cohorts=COHORTS, endpoints_to_compare=ENDPOINTS, methods=METHODS):
    summary, differences, donors_out, strata = [], [], [], []
    for cohort, endpoint in itertools.product(cohorts, endpoints_to_compare):
        part = [r for r in rows if r["cohort"] == cohort and r["fine_endpoint"] == endpoint]
        groups = {
            method: {r["sample_id"]: r for r in part if r["method"] == method}
            for method in methods
        }
        query = groups[methods[-1]]
        donors = sorted({r["donor"] for r in query.values()})
        for method, sample_rows in groups.items():
            values = list(sample_rows.values())
            region_means = []
            for region in sorted({row["region_original"] for row in values}):
                region_rows = [r for r in values if r["region_original"] == region]
                region_donors = sorted({r["donor"] for r in region_rows})
                result = {key: mean(mean(r[key] for r in region_rows if r["donor"] == donor)
                                    for donor in region_donors) for key in METRICS}
                region_means.append(result)
                strata.append(dict(cohort=cohort, endpoint=endpoint, method=method,
                                   region_original=region, n=len(region_rows), donors=len(region_donors),
                                   **result))
            donor_means = {}
            for donor in donors:
                donor_rows = [r for r in values if r["donor"] == donor]
                donor_means[donor] = {key: mean(r[key] for r in donor_rows) for key in METRICS}
                donors_out.append(dict(cohort=cohort, endpoint=endpoint, method=method,
                                       donor=donor, n=len(donor_rows), **donor_means[donor]))
            info = dict(
                cohort=cohort, endpoint=endpoint, method=method, n=len(values), donors=len(donors),
                original_anatomy_labels=len(region_means),
                n_allowed_in_pool=sum(r["allowed_in_pool"] for r in values),
                n_allowed_outside_pool=sum(1-r["allowed_in_pool"] for r in values),
                n_exact_failure=sum(bool(r["Exact_failure"]) for r in values),
                n_group_failure=sum(bool(r["Group_failure"]) for r in values),
                n_not_run=sum(r["inference_status"] == "not_run" for r in values),
            )
            for metric in METRICS:
                info["donor_equal_" + metric] = mean(row[metric] for row in donor_means.values())
                info["anatomy_label_equal_" + metric] = mean(row[metric] for row in region_means)
                info["pooled_" + metric] = mean(r[metric] for r in values)
                if metric.endswith(("hit1", "hit3")):
                    info["pooled_hits_" + metric] = sum(r[metric] for r in values)
            summary.append(info)
        for comparison in methods[:-1]:
            for metric in METRICS[:6]:
                pair = [(r[metric], groups[comparison][sample][metric], r["donor"])
                        for sample, r in query.items()]
                by_donor = [mean(a - b for a, b, d in pair if d == donor) for donor in donors]
                differences.append(dict(
                    cohort=cohort, endpoint=endpoint, comparator=comparison, metric=metric,
                    n=len(pair), donors=len(donors), pearson_minus_comparator_donor_equal=mean(by_donor),
                    samples_pearson_higher=sum(a > b for a, b, _ in pair),
                    samples_pearson_lower=sum(a < b for a, b, _ in pair),
                    samples_same=sum(a == b for a, b, _ in pair),
                    donors_pearson_higher=sum(v > 1e-12 for v in by_donor),
                    donors_pearson_lower=sum(v < -1e-12 for v in by_donor),
                    donors_same=sum(abs(v) <= 1e-12 for v in by_donor),
                ))
    return summary, differences, donors_out, strata


def main():
    start = monotonic()
    cfg, base, output = arguments(__doc__)
    source_archive = location(cfg, base, "source_archive")
    membership_file = location(cfg, base, "training_membership")
    region_parent_file = location(cfg, base, "region_parent")
    beam_root = location(cfg, base, "beam_rules_root")
    beam_manifest = location(cfg, base, "beam_manifest")
    reference, genes, parents, train = source_reference(source_archive, membership_file, region_parent_file)
    beams = frozen_beams(reference, parents, train, [
        (record, beam_root / record["id"]) for record in read_json(beam_manifest)
    ])
    admission = read_csv(location(cfg, base, "anatomy_admission"))
    all_heads = read_csv(location(cfg, base, "fixed_head"))
    all_rows, disposition, target_info = [], [], {}
    for cohort in COHORTS:
        roster = [r for r in admission if r["cohort"] == cohort and r["fine_endpoint"] in ENDPOINTS]
        heads = {r["sample_id"]: r for r in all_heads if r["cohort"] == cohort}
        target_config = cfg["cohorts"][cohort]
        prior = {(r["sample_id"], r["arm"]): r for r in read_csv(
            location(target_config, base, "old_predictions"))}
        query_ids = [r["sample_id"] for r in roster]
        with np.load(location(target_config, base, "target_input"), allow_pickle=False) as target:
            fail_unless(target["genes"].tolist() == genes, f"Frozen gene axis mismatch: {cohort}")
            observed = target["observed"]
            samples = target["samples"].tolist()
            scores = query_scores(target["x"], observed, samples, query_ids, heads, beams, reference)
        target_info[cohort] = dict(
            frozen_matrix_samples=len(samples), compared_samples=len(roster),
            reported_donors=len({r["donor"] for r in roster}),
            beams_used=len({r["beam_id"] for r in scores.values()}),
            measured_genes_in_G=int(observed[reference.gene_indices].sum()),
            imputed_genes_in_G=int((~observed[reference.gene_indices]).sum())
        )
        for row in (r for r in admission if r["cohort"] == cohort):
            disposition.append(dict(
                **row, round2_scored=row["sample_id"] in scores,
                round2_disposition=("evaluated_on_original_stratum" if row["sample_id"] in scores
                                    else "original_no_eligible_endpoint")
            ))
        all_rows.extend(comparable_rows(roster, prior, scores))
        print(f"{cohort}: {len(roster)} specimens, {len(heads)} frozen head rows", flush=True)
    summary, differences, donor_rows, anatomy_rows = summaries(all_rows)
    write_csv(output / "PER_SPECIMEN.csv", all_rows)
    write_csv(output / "DONOR_AND_POOLED_SUMMARY.csv", summary)
    write_csv(output / "PAIRED_DIFFERENCES.csv", differences)
    write_csv(output / "PER_DONOR.csv", donor_rows)
    write_csv(output / "PER_ANATOMY_LABEL.csv", anatomy_rows)
    write_csv(output / "ORIGINAL_DISPOSITION.csv", disposition)
    metadata = dict(
        method=METHODS[2], source_archive=str(source_archive), source_members_read=[
            "x", "genes", "samples", "regions", "labels", "donors"],
        source_specimens=len(train), source_donors=len({r["donor"] for r in train}),
        source_G=len(reference.gene_indices), source_networks=len(reference.network_labels),
        source_fine_regions=len(reference.fine_labels), source_beams_bound=len(beams),
        targets=target_info, total_original_disposition=len(disposition),
        evaluated_specimens=sum(info["compared_samples"] for info in target_info.values()),
        input_config=cfg, reused_pearson_source_code="common/pearson_hierarchy.py",
        exact_score="donor-equal regional centroid Pearson on fixed G, candidate-restricted",
        group_score="maximum member correlation using frozen beam-specific Group membership",
        ties="descending score and ascending UTF-8 label bytes at both tiers",
        python=platform.python_version(), numpy=np.__version__, elapsed_seconds=monotonic()-start
    )
    (output / "RUN.json").write_text(json.dumps(metadata, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
