"""Train the R9 macaque source Network, build its lower-tier references and evaluate donor holdouts."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from hashlib import file_digest, sha256
from pathlib import Path
import json
import time

import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits

from lower_tiers import FineReferences, write_json
from network import NetworkSplit, all_scores, choose_c
from source_data import read_source


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--science-root", type=Path, required=True,
                   help="Input directory containing external_data/Bo2023 with the original count matrix and specimen workbook.")
    p.add_argument("--output", type=Path, required=True, help="New directory for predictions, selection scores and reference assets.")
    p.add_argument("--variant", choices=("floor_clip","baseline","floor_only","clip_only"), default="floor_clip")
    p.add_argument("--workers", type=int, default=2)
    p.add_argument("--selection-scores", type=Path, help="Optional input-bound audit cache JSON for formal floor_clip.")
    return p.parse_args()


def cached_scores(path, repo, variant, data):
    cache = json.loads(path.read_text(encoding="utf8"))
    annotation = Path(__file__).with_name("mfas5_ensembl97_gene_map.csv")
    if cache.get("annotation_sha256") != sha256(annotation.read_bytes()).hexdigest():
        raise ValueError("Selection scores were not generated with this source annotation")
    if cache["variant"] != variant:
        raise ValueError("Selection scores belong to a different Network variant")
    for name, expected in cache["source_sha256"].items():
        with (repo / name).open("rb") as f:
            if file_digest(f, "sha256").hexdigest() != expected:
                raise ValueError(f"Selection cache differs from the input {name}")
    scores = pd.DataFrame(cache["scores"])
    donors = set(data.donors)
    expected = {(outer, validation, c) for outer in donors | {""} for validation in donors if validation != outer
                for c in (0.001,0.01,0.1)}
    actual = set(scores[["outer_holdout","validation_donor","C"]].itertuples(index=False,name=None))
    if actual != expected or len(scores) != len(expected):
        raise ValueError("Selection cache lacks the complete donor/C grid")
    return scores


def write_universe(folder, data, split):
    write_json(folder / "eligible_genes.json", data.symbols[split.eligible].tolist())


def evaluate_donor(data, donor, number, scores, variant, output):
    from core.scientific_validation import signal_reason

    folder = output / "folds" / f"outer_{number:02d}"
    folder.mkdir()
    split = NetworkSplit(data, (donor,), variant)
    c = choose_c(scores, donor)
    model = split.fit(c)
    prob, ranks = split.score(model)
    split.export(folder / "network", model)
    write_universe(folder, data, split)
    ref = FineReferences(data, split.train, split.eligible, folder / "lower")
    records = []
    for pos, idx in enumerate(split.query):
        network_truth = data.networks[idx]
        fine_truth = data.regions[idx]
        values = data.x[split.panel,idx]
        why = signal_reason(split.genes.tolist(), pd.Series(values,index=split.genes), values, prob[pos])
        if why:
            raise ValueError(f"Source specimen {data.samples[idx]} has no admissible Network signal: {why}")
        ranking = ranks[pos].tolist()
        gate = ranking[:3]
        rank = ranking.index(network_truth)+1
        row = dict(arm="floor_clip_hybrid" if variant=="floor_clip" else variant, fold=number, sample_id=data.samples[idx],
                   donor=donor, truth_network=network_truth, truth_region=fine_truth, fine_supported=fine_truth in ref.parents,
                   fine_support_reason="" if fine_truth in ref.parents else "truth_region_absent_from_training_donors",
                   Network_valid=True, Network_failure="", Network_ranking_json=json.dumps(ranking),
                   Network_hit1=int(rank==1), Network_hit3=int(rank<=3), Network_MRR=1/rank,
                   new_head_selected_C=c, network_origin="ensembl97_training_fold_detected_symbols",
                   beam_json=json.dumps(gate), truth_network_in_beam=network_truth in gate,
                   Network_scores_json=json.dumps(prob[pos].tolist()))
        row.update(ref.evaluate(idx,gate))
        records.append(row)
    ref.save_index()
    write_json(folder / "fit.json", {"C":c,"held_out_donor":donor,"train_sample_ids":data.samples[split.train].tolist(),
                                    "query_sample_ids":data.samples[split.query].tolist(),"eligible_symbol_count":len(split.eligible)})
    print(f"Outer fold {number}: donor {donor}, C={c:g}, {len(records)} source holdouts",flush=True)
    return records


def main():
    args = parse_args()
    repo = args.science_root.resolve()
    output = args.output.resolve()
    output.mkdir(parents=True,exist_ok=False)
    (output/"folds").mkdir()
    start = time.monotonic()
    data = read_source(repo)
    with threadpool_limits(limits=1):
        if args.selection_scores:
            scores = cached_scores(args.selection_scores,repo,args.variant,data)
        else:
            scores = all_scores(data,args.variant,args.workers)
        scores.to_csv(output/"selection_scores.csv",index=False)
        donors = sorted(set(data.donors))
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            pending = [pool.submit(evaluate_donor,data,d,i,scores,args.variant,output) for i,d in enumerate(donors,1)]
            rows = [row for work in pending for row in work.result()]
        predictions = pd.DataFrame(rows)
        predictions.to_csv(output/"source_predictions.csv",index=False)
        full = NetworkSplit(data, (), args.variant)
        final_c = choose_c(scores,"")
        final = output/"final"
        final.mkdir()
        final_model = full.fit(final_c)
        full.export(final/"network",final_model)
        write_universe(final,data,full)
        np.save(final/"SOURCE_RAW_MEAN.npy",data.x[full.eligible].mean(axis=1,dtype=np.float64).astype(np.float32))
        refs = FineReferences(data, full.train, full.eligible, final/"lower")
        refs.complete_final(final_model.classes_)
        refs.save_index()
    endpoints = [f"{tier}_hit{k}" for tier in ("Network","Group","Exact") for k in (1,3)]
    summary = dict(variant=args.variant, sample_count=len(rows), donors=len(donors),
                   selection="bound audit scores" if args.selection_scores else "fits from raw count input",
                   outer_C=predictions.groupby("donor").new_head_selected_C.first().to_dict(),final_C=final_c,
                   pooled_hits=predictions[endpoints].sum().to_dict(),
                   donor_equal=predictions.groupby("donor")[endpoints].mean().mean().to_dict(),
                   final_eligible_symbols=len(full.eligible), final_lower_combinations=len(refs.index),
                   seconds=round(time.monotonic()-start,1))
    write_json(output/"summary.json",summary)
    print(json.dumps(summary,ensure_ascii=False),flush=True)


if __name__ == "__main__":
    main()

