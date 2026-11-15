"""SingleR on the frozen human RANK200 source and GTEx/HPA adapter vectors."""
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

import numpy as np
import pandas as pd
import scipy
from scipy.stats import rankdata

COHORTS = {
    "GTEx": ("percentile_ranks", "model_gene_ids"),
    "HPA": ("feature_ranks", "feature_gene_ids"),
}
INPUT_PATHS = (
    "source_profile_binding", "source_profiles", "training_axis", "matched_transform",
    "saved_comparator_predictions", "gtex_rank_input", "hpa_rank_input",
)
METHOD = "SingleR_shared_RANK200"
METRICS = ["Top1", "Top3", "MRR"]


def prepare(temp, originals, paths):
    source_meta = json.loads(paths["source_profile_binding"].read_text(encoding="utf-8"))
    reference_labels = pd.read_csv(paths["training_axis"], sep="\t", dtype=str, quoting=3)
    if reference_labels.n_specimens.ne("1").any() or len(reference_labels) != 1083:
        raise ValueError("Final PLM columns must be singleton donor-by-label profiles")
    with np.load(paths["matched_transform"], allow_pickle=False) as model:
        genes = model["selected_gene_ids"].astype(str)
        columns = model["selected_U_indices"]
        classes = sorted(model["class_ids"].astype(str).tolist())
        mean, scale = model["mean"], model["scale"]
    if not np.array_equal(np.asarray(source_meta["gene_ids"])[columns], genes):
        raise ValueError("Final PLM and the RANK200 source axis disagree")
    if reference_labels.ref_id.tolist() != source_meta["donor_label_ids"] or sorted(reference_labels.label.unique()) != classes:
        raise ValueError("Final PLM training profile identifiers or labels disagree")
    source = np.memmap(paths["source_profiles"], mode="r", dtype="<f8",
                       shape=tuple(source_meta["profile_shape"]))
    rank = (rankdata(source, axis=1, method="average")[:, columns] - 1) / (source.shape[1] - 1)
    standardized = (rank - mean) / scale
    pd.DataFrame(standardized.T, index=genes, columns=reference_labels.ref_id).to_csv(
        temp / "reference.tsv", sep="\t", index_label="gene_id")
    reference_labels[["ref_id", "label"]].to_csv(temp / "reference_labels.tsv", sep="\t", index=False)

    for cohort, (key, gene_key) in COHORTS.items():
        saved = originals[(originals.dataset == cohort) & (originals.method == "BrainTrace_RANK")]
        with np.load(paths[f"{cohort.lower()}_rank_input"], allow_pickle=False) as query:
            features = pd.DataFrame(query[key], index=query["sample_ids"].astype(str),
                                    columns=query[gene_key].astype(str))
        if set(features.index) != set(saved.sample_id) or set(features.columns) != set(genes):
            raise ValueError(f"{cohort}: query axes disagree with original endpoint")
        query = (features.loc[saved.sample_id, genes].to_numpy() - mean) / scale
        pd.DataFrame(query.T, index=genes, columns=saved.sample_id).to_csv(
            temp / f"{cohort}_query.tsv", sep="\t", index_label="gene_id")
    return classes


def compare(cohort, native, old, classes):
    saved = old[old.method == "BrainTrace_RANK"].set_index("sample_id")
    frame = saved[["dataset", "donor", "truth", "region", "subregion", "tier"]].join(
        native.set_index("sample_id"), how="left", validate="one_to_one")
    rows = []
    for sample, row in frame.iterrows():
        failure = row.reason if row.reason else ""
        scores = {k: row[k] for k in classes}
        ranking = sorted(classes, key=lambda c: (-scores[c], c)) if not failure else []
        rank = ranking.index(row.truth) + 1 if ranking else 0
        rows.append({"dataset": cohort, "sample_id": sample, "donor": row.donor,
                     "truth": row.truth, "region": row.region, "subregion": row.subregion,
                     "tier": row.tier, "method": METHOD, "origin": "new_comparator_prediction",
                     "evaluable": not failure, "failure_reason": failure,
                     "predicted_class": ranking[0] if ranking else "FAIL", "true_rank": rank,
                     "Top1": int(rank == 1), "Top3": int(0 < rank <= 3),
                     "MRR": 1 / rank if rank else 0, "ranking": "|".join(ranking),
                     "native_label": row.native_label,
                     **{f"score_{k}": scores[k] for k in classes}})
    return pd.DataFrame(rows)


def summary(cohort, data):
    tiers = [(stratum, data.tier == stratum) for stratum in data.tier.unique()]
    tiers.append(("all_six_tissues_descriptive" if cohort == "GTEx" else "combined_descriptive",
                  pd.Series(True, index=data.index)))
    if cohort == "GTEx":
        tiers.extend((f"tissue_{t}", data.region == t) for t in sorted(data.region.unique()))
    else:
        tiers.append(("provisional_frontal_parent", (data.tier == "supplementary_subregions")
                      & (data.truth == "4009")))
    result = []
    pairs = []
    for stratum, mask in tiers:
        strata = data[mask]
        for method, rows in strata.groupby("method", sort=False):
            pooled = rows[METRICS].mean()
            weighted = rows.groupby("donor")[METRICS].mean().mean() if cohort == "GTEx" else pooled
            result.append({"dataset": cohort, "stratum": stratum, "method": method,
                           "samples": len(rows), "donors": rows.donor.nunique() if cohort == "GTEx" else None,
                           "failures": int((~rows.evaluable).sum()), "Top1_hits": int(rows.Top1.sum()),
                           "Top3_hits": int(rows.Top3.sum()), **weighted.to_dict(),
                           **{f"pooled_{k}": v for k, v in pooled.items()}})
        base = strata[strata.method == "BrainTrace_RANK"].set_index("sample_id")
        sr = strata[strata.method == METHOD].set_index("sample_id").loc[base.index]
        for metric in METRICS:
            delta = sr[metric] - base[metric]
            row = {"dataset": cohort, "stratum": stratum, "metric": metric,
                   "samples": len(base), "donors": base.donor.nunique() if cohort == "GTEx" else None,
                   "SingleR_minus_BrainTrace": (delta.groupby(base.donor).mean().mean()
                                                if cohort == "GTEx" else delta.mean())}
            if metric != "MRR":
                row.update({"both_hit": int(((sr[metric] == 1) & (base[metric] == 1)).sum()),
                            "SingleR_only": int((delta == 1).sum()),
                            "BrainTrace_only": int((delta == -1).sum()),
                            "both_miss": int(((sr[metric] == 0) & (base[metric] == 0)).sum())})
            pairs.append(row)
    return result, pairs


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True, help="JSON file binding frozen inputs and R")
    parser.add_argument("--output", type=Path, required=True, help="New output directory")
    parser.add_argument("--work-dir", type=Path, help="Parent for temporary R workspace")
    args = parser.parse_args()

    config_file = args.config.resolve()
    config = json.loads(config_file.read_text(encoding="utf-8"))
    paths = {key: (config_file.parent / config[key]).resolve() for key in INPUT_PATHS}
    rscript = shutil.which(config["rscript"]) if Path(config["rscript"]).name == config["rscript"] else (
        config_file.parent / config["rscript"]).resolve()
    if rscript is None:
        parser.error(f"Rscript executable not on PATH: {config['rscript']}")
    out = args.output.resolve()
    out.mkdir(parents=True)

    originals = pd.read_csv(paths["saved_comparator_predictions"], keep_default_na=False,
                            dtype={"sample_id": str, "truth": str, "predicted_class": str})
    summaries, paired, predictions = [], [], []
    with tempfile.TemporaryDirectory(prefix="bt_s31_", dir=args.work_dir) as directory:
        temp = Path(directory).resolve()
        classes = prepare(temp, originals, paths)
        shutil.copyfile(Path(__file__).with_name("score_singler.R"), temp / "score_singler.R")
        env = os.environ.copy()
        for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
            env[name] = "1"
        for name in ("TMPDIR", "TEMP", "TMP"):
            env[name] = str(temp)
        if config.get("r_library"):
            env["R_LIBS_USER"] = str((config_file.parent / config["r_library"]).resolve())
        with (out / "NATIVE_RUN.log").open("w", encoding="utf-8") as log:
            subprocess.run([str(rscript), "--vanilla", "score_singler.R", "."], cwd=temp,
                           env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
        for name in ("marker_genes.txt", "R_SESSION.txt"):
            (out / name).write_bytes((temp / name).read_bytes())
        for cohort in COHORTS:
            file = f"{cohort}_native.tsv"
            native = pd.read_csv(temp / file, sep="\t", float_precision="round_trip",
                                 dtype={"sample_id": str, "native_label": str})
            native[["reason", "native_label"]] = native[["reason", "native_label"]].fillna("")
            old = originals[originals.dataset == cohort]
            if set(native.sample_id) != set(old.sample_id):
                raise ValueError(f"{cohort}: scoring lost original specimen IDs")
            new = compare(cohort, native, old, classes)
            predictions.append(new)
            combined = pd.concat([old, new], ignore_index=True)
            stats, contrasts = summary(cohort, combined)
            summaries.extend(stats)
            paired.extend(contrasts)
    pred = pd.concat(predictions, ignore_index=True)
    pred.to_csv(out / "SINGLER_PREDICTIONS.csv", index=False)
    pred[~pred.evaluable].to_csv(out / "FAILURES.csv", index=False)
    pd.DataFrame(summaries).to_csv(out / "ENDPOINTS.csv", index=False)
    pd.DataFrame(paired).to_csv(out / "PAIRED_WITH_BRAINTRACE.csv", index=False)
    manifest = {
        "input_config": str(config_file),
        "source_profile_binding": str(paths["source_profile_binding"]),
        "source_profiles": str(paths["source_profiles"]),
        "training_axis": str(paths["training_axis"]),
        "matched_transform": str(paths["matched_transform"]),
        "query_files": {c: str(paths[f"{c.lower()}_rank_input"]) for c in COHORTS},
        "saved_comparator_predictions": str(paths["saved_comparator_predictions"]),
        "python": sys.version.split()[0], "numpy": np.__version__, "pandas": pd.__version__,
        "scipy": scipy.__version__,
        "classes": classes, "source_profiles": 1083, "source_donors": 132,
        "new_SingleR_queries": len(pred), "new_SingleR_fits": 1,
        "query_failures": int((~pred.evaluable).sum()),
        "native_markers": len((out / "marker_genes.txt").read_text().splitlines()),
        "native_score_export_significant_digits": 17,
        "new_p_values_intervals": False,
    }
    (out / "RUN.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")


if __name__ == "__main__":
    main()
