"""Extend the shared-head lower comparison to the existing adult PsychENCODE strata."""

from collections import Counter
import json
from time import monotonic
 
import numpy as np

from ._config import arguments, location
from .shared_head_external import (
    METHODS, read_csv, read_json, write_csv, field, fail_unless,
    source_reference, frozen_beams, query_scores, endpoints, summaries
)

COHORT = "PsychENCODE_adult"
SETS = ("unique_operational", "ambiguous_set", "parent_compatible")
ARMS = (METHODS[0], METHODS[2])


def bind_beams(beam_root, beam_manifest):
    return [
        (record, beam_root / record["id"])
        for record in read_json(beam_manifest)
    ]


def paired_rows(roster, old, flat, scores):
    result = []
    for truth in roster:
        sid = truth["sample_id"]
        saved = old[sid]
        saved_hits = flat[sid]
        new = scores[sid]
        allowed = set(json.loads(truth["lower_allowed_json"]))
        fail_unless(new["head"] == saved["network_ranking"] and new["beam_id"] == saved["beam_id"]
                    and saved["candidate_regions"] == new["candidates"], f"Frozen candidates differ: {sid}")
        fail_unless(all(saved_hits[key] == truth[key] for key in
                        ("donor", "lower_set", "lower_allowed_json")), f"Frozen anatomy differs: {sid}")
        common = dict(
            cohort=COHORT, sample_id=sid, donor=truth["donor"], region_original=truth["region"],
            fine_endpoint=truth["lower_set"], allowed_regions_json=truth["lower_allowed_json"],
            allowed_set_size=len(allowed), beam_id=new["beam_id"], network_top3_json=field(new["head"][:3]),
            candidates_json=field(new["candidates"]), groups_json=field(new["groups"]),
            allowed_in_pool=int(bool(allowed.intersection(new["candidates"])))
        )
        old_exact = saved["Exact_ranking"] if saved["Exact_valid"] else []
        old_groups = saved["Group_members"] if saved["Group_valid"] else []
        old_group_labels = saved["Group_ranking"] if old_groups else []
        fail_unless((not old_exact or sorted(old_exact) == new["candidates"]) and
                    (not old_groups or {key: sorted(vals) for key, vals in
                     zip(old_group_labels, old_groups, strict=True)} == new["groups"]),
                    f"Saved ranking or Group partition differs: {sid}")
        metrics = endpoints(old_exact, old_groups, allowed)
        for layer in ("Exact", "Group"):
            for k in (1, 3):
                fail_unless(int(saved_hits[f"{layer}_Top{k}"]) == metrics[f"{layer}_hit{k}"],
                            f"Frozen hit calculation differs: {sid}:{layer}:Top{k}")
        result.append(dict(
            **common, method=ARMS[0], inference_status="predicted",
            Exact_failure=saved["Exact_failure"], Group_failure=saved["Group_failure"],
            exact_ranking_json=field(old_exact), exact_scores_json="",
            group_ranking_json=field(old_group_labels), group_ranking_members_json=field(old_groups),
            group_ranking_scores_json=field(saved["Group_scores"]), **metrics
        ))
        exact = [entry["label"] for entry in new["exact"]]
        result.append(dict(
            **common, method=ARMS[1], inference_status="predicted",
            Exact_failure=new["Exact_failure"], Group_failure=new["Group_failure"],
            exact_ranking_json=field(exact), exact_scores_json=field([entry["score"] for entry in new["exact"]]),
            group_ranking_json=field(new["group_labels"]), group_ranking_members_json=field(new["group_members"]),
            group_ranking_scores_json=field(new["group_scores"]), **endpoints(exact, new["group_members"], allowed)
        ))
    return result


def main():
    start = monotonic()
    cfg, base, output = arguments(__doc__)
    source_archive = location(cfg, base, "source_archive")
    membership_file = location(cfg, base, "training_membership")
    region_parent_file = location(cfg, base, "region_parent")
    beam_root = location(cfg, base, "beam_rules_root")
    beam_manifest = location(cfg, base, "beam_manifest")
    target_input = location(cfg, base, "target_input")
    original_labels = location(cfg, base, "original_labels")
    original_predictions = location(cfg, base, "old_predictions_json")
    original_rankings_csv = location(cfg, base, "old_predictions_csv")
    source_means = location(cfg, base, "source_raw_mean")
    reference, genes, parents, train = source_reference(source_archive, membership_file, region_parent_file)
    beams = frozen_beams(reference, parents, train, bind_beams(beam_root, beam_manifest))
    original_roster = read_csv(original_labels)
    roster = [r for r in original_roster if r["lower_set"] in SETS]
    old = {r["sample_id"]: r for r in read_json(original_predictions)}
    flat = {r["sample_id"]: r for r in read_csv(original_rankings_csv)}
    heads = {sid: {"ranking_json": field(r["network_ranking"])} for sid, r in old.items()}
    with np.load(target_input, allow_pickle=False) as input_file:
        fail_unless(input_file["genes"].tolist() == genes, "PsychENCODE gene axis differs from B5")
        samples = input_file["samples"].tolist()
        observed = input_file["observed"]
        x = input_file["x"].T.copy()
    mean = np.load(source_means, allow_pickle=False)
    x[~observed, :] = mean[~observed, None]
    scores = query_scores(x, observed, samples, [r["sample_id"] for r in roster], heads, beams, reference)
    rows = paired_rows(roster, old, flat, scores)
    summary, contrasts, donors, strata = summaries(rows, (COHORT,), SETS, ARMS)
    write_csv(output / "PSYCHENCODE_PER_SPECIMEN.csv", rows)
    write_csv(output / "PSYCHENCODE_DONOR_AND_POOLED_SUMMARY.csv", summary)
    write_csv(output / "PSYCHENCODE_PAIRED_DIFFERENCES.csv", contrasts)
    write_csv(output / "PSYCHENCODE_PER_DONOR.csv", donors)
    write_csv(output / "PSYCHENCODE_PER_ANATOMY_LABEL.csv", strata)
    write_csv(output / "PSYCHENCODE_ORIGINAL_DISPOSITION.csv", [
        dict(**row, round2_scored=row["sample_id"] in scores) for row in original_roster
    ])
    metadata = dict(
        source_specimens=len(train), source_donors=len({r["donor"] for r in train}),
        source_fine_regions=len(reference.fine_labels), source_G=len(reference.gene_indices),
        adult_prepared_input=str(target_input),
        original_predictions=str(original_predictions),
        original_labels=str(original_labels), frozen_lower_assets=str(beam_root.parent), input_config=cfg,
        comparison_rows_by_stratum=dict(Counter(r["lower_set"] for r in roster)),
        frozen_input_samples=len(samples), compared_samples=len(roster),
        reported_donors=len({r["donor"] for r in roster}),
        beams_used=len({r["beam_id"] for r in scores.values()}),
        measured_genes_in_G=int(observed[reference.gene_indices].sum()),
        imputed_genes_in_G=int((~observed[reference.gene_indices]).sum()),
        normalization="ln(1+CPM), using the 26,514 public common-annotation count rows as denominator",
        elapsed_seconds=monotonic()-start
    )
    (output / "PSYCHENCODE_RUN.json").write_text(
        json.dumps(metadata, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"Adult PsychENCODE: {len(roster)}/{len(samples)} specimens scored on original lower strata")


if __name__ == "__main__":
    main()
