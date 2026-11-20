"""Recalculate the source rows of S41B and the paired S38 donor diagnostics."""
import numpy as np
import pandas as pd

from source_tables import OUT, enrich, read_table, source_subsets


CATEGORIES = (
    "truth_absent_from_source_training_reference",
    "upstream_output_unavailable",
    "allowed_fine_truth_absent_from_saved_Network_pool",
    "lower_output_unavailable",
    "truth_in_pool_but_ranking_miss",
    "hit",
)


def failure_attribution(data, cohort="Source_OOF", sets=None):
    by_donor = []
    aggregates = []
    for endpoint, subset in (source_subsets(data) if sets is None else sets):
        for tier, hit_col, valid_col in (("Group", "Group", "group_valid"), ("Exact", "exact_observed", "exact_valid")):
            for k, group in subset.groupby("k"):
                category = np.select(
                    [~group.fine_supported, group.K.eq(0), ~group.in_pool, ~group[valid_col], group[hit_col].eq(0)],
                    CATEGORIES[:-1], default="hit",
                )
                counts = pd.crosstab(group.donor, category).reindex(columns=CATEGORIES, fill_value=0)
                n = counts.sum(axis=1)
                for label in CATEGORIES:
                    base = dict(cohort=cohort, evaluation_set=endpoint, tier=tier, k=k, category=label)
                    percentages = 100 * counts[label] / n
                    for donor in counts.index:
                        by_donor.append(dict(**base, donor=donor, samples=n[donor], count=counts.at[donor, label],
                                             donor_percent=percentages[donor]))
                    total = counts[label].sum()
                    aggregates.append(dict(**base, samples=len(group), donors=len(counts), count=total,
                                           pooled_percent=100*total/len(group), donor_equal_percent=percentages.mean()))
    return pd.DataFrame(by_donor), pd.DataFrame(aggregates)


def retrieval_contrasts(data, cohort="Source_OOF", sets=None):
    by_donor = []
    omissions = []
    for endpoint, subset in (source_subsets(data) if sets is None else sets):
        for k, group in subset.groupby("k"):
            for control in ("FlatPearson200", "Exact50_100"):
                base = dict(cohort=cohort, evaluation_set=endpoint, k_groups=k, control=control)
                donor_differences = {}
                for donor, paired in group.groupby("donor"):
                    group_hit = paired.Group.astype(bool)
                    control_hit = paired[control].astype(bool)
                    difference = 100 * (group_hit.mean() - control_hit.mean())
                    donor_differences[donor] = difference
                    by_donor.append(dict(**base, donor=donor, samples=len(paired), mean_candidate_budget=paired.m.mean(),
                                         group_hit_count=group_hit.sum(), control_hit_count=control_hit.sum(),
                                         both_hit=(group_hit & control_hit).sum(), group_only=(group_hit & ~control_hit).sum(),
                                         control_only=(~group_hit & control_hit).sum(), neither=(~group_hit & ~control_hit).sum(),
                                         group_hit_percent=100*group_hit.mean(), control_hit_percent=100*control_hit.mean(),
                                         group_minus_control_pp=difference))
                differences = pd.Series(donor_differences)
                for donor in differences.index:
                    retained = differences.drop(donor)
                    omissions.append(dict(**base, omitted_donor=donor, remaining_donors=len(retained),
                                          paired_difference_pp=retained.mean()))
    return pd.DataFrame(by_donor), pd.DataFrame(omissions)


def main():
    predictions = read_table(OUT.parent / "source_analysis" / "results" / "source_predictions.csv")
    data = enrich(predictions.set_index("sample_id"))
    category_donors, categories = failure_attribution(data)
    paired_donors, omissions = retrieval_contrasts(data)
    for name, frame in [
        ("S41B_SOURCE_BY_DONOR", category_donors), ("S41B_SOURCE_CATEGORY_SHARES", categories),
        ("S38_SOURCE_PAIRED_BY_DONOR", paired_donors), ("S38_SOURCE_DONOR_OMISSION", omissions),
    ]:
        frame.to_csv(OUT / (name + ".csv"), index=False)
    print(categories.pivot(index=["evaluation_set","tier","k"], columns="category", values="donor_equal_percent").round(3))


if __name__ == "__main__":
    main()
