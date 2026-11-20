"""Reproject original complete libraries onto the Ensembl 97 macaque reference."""
from collections import defaultdict
import csv
import gzip
import importlib.util
import json
import os
from pathlib import Path
import sys

import joblib
import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits

HERE = Path(__file__).resolve().parent
PROJECT = Path(os.environ.get("BRAINTRACE_STUDY_ROOT", Path.cwd())).resolve()
REPO = PROJECT / "BrainTrace_algorithm_research_20261015"
RUNS = REPO / "research/runs"
PREP = RUNS / "macaque_external_preparation_20261031_01"
GEO = RUNS / "macaque_external_lower_validation_20261101_01"
PSYCH = PROJECT / "BrainTrace_PsychENCODE_Adult_Validation_20261101"
SOURCE = HERE.parent / "source_analysis"
FINAL = SOURCE / "results/final"
TOOL = HERE.parents[3] / "software/science"
OUT = HERE / "results"


def js(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def source_ids():
    table = pd.read_csv(SOURCE / "mfas5_ensembl97_gene_map.csv", dtype=str, keep_default_na=False)
    names = table.set_index("Gene.stable.ID")["Gene.name"]
    path = REPO / "external_data/Bo2023/mfas5_819samples_28415genes_featurecounts_counts.txt"
    by_symbol = defaultdict(set)
    with path.open(encoding="utf8") as stream:
        next(stream)
        for line in stream:
            sid = line.split("\t", 1)[0]
            by_symbol[names.at[sid]].add(sid)
    return by_symbol


def strict_map(features, ids, homology, index):
    users = defaultdict(set)
    mapping = []
    for gene in features:
        source = sorted(ids[gene])
        targets, issues = [], []
        for sid in source:
            edges = homology.get(sid, [])
            unique = {e["homology_gene_stable_id"] for e in edges}
            types = {e["homology_type"] for e in edges}
            if len(unique) != 1 or types != {"ortholog_one2one"}:
                issues.append("non_one2one_or_ambiguous" if edges else "no_homology")
            else:
                target = next(iter(unique))
                if target in index:
                    targets.append(target)
                else:
                    issues.append("ortholog_absent_from_matrix")
        if len(set(targets)) != len(targets):
            issues.append("within_feature_target_collision")
        if not source:
            raise ValueError(f"No source count IDs for {gene}")
        row = dict(feature=gene, source_ids=source, target_ids=sorted(set(targets)), issues=sorted(set(issues)), covered=not issues)
        if row["covered"]:
            for target in targets:
                users[target].add(gene)
        mapping.append(row)
    for row in mapping:
        if row["covered"] and any(len(users[t]) > 1 for t in row["target_ids"]):
            row["covered"] = False
            row["issues"].append("between_feature_target_collision")
    return mapping


def raw_input(cohort):
    if cohort == "PsychENCODE":
        path = PROJECT / "BrainTrace_Public_Macaque_Screen_20261101/sources/simFiltered.HPR.compGene.count.txt"
        with path.open(encoding="utf-8-sig") as f:
            header = next(f).rstrip().split("\t")
        columns = [i for i,s in enumerate(header) if i == 0 or s.startswith("RMB")]
        data = pd.read_csv(path, sep="\t", usecols=columns, index_col=0)
        data.index = data.index.str.split("|", regex=False).str[0]
        hpath = PSYCH / "ENSEMBL91_FASCICULARIS_HUMAN.json"
    else:
        if cohort == "GSE128537":
            path = RUNS / "macaque_external_discovery_20261031_01/downloads/GSE128537_Rhesus_RNA_gene_count.txt.gz"
            version = 91
        else:
            path = PREP / "downloads/GSE179328_est_counts.txt.gz"
            version = 101
        with gzip.open(path, "rt", encoding="utf8") as f:
            header = next(f).rstrip("\r\n").split("\t")[1:]
        data = pd.read_csv(path, sep="\t", header=None, skiprows=1, index_col=0, float_precision="round_trip")
        data.columns = header
        hpath = PREP / f"ENSEMBL_{version}_FASCICULARIS_MULATTA_HOMOLOGY.json"
    homology = defaultdict(list)
    for row in js(hpath):
        homology[row["gene_stable_id"]].append(row)
    return data, homology


def project(data, mapping):
    counts = np.ascontiguousarray(data.to_numpy(dtype=np.float64))
    library = counts.sum(axis=0)
    index = {gene:i for i,gene in enumerate(data.index)}
    x = np.full((len(data.columns),len(mapping)), np.nan, dtype=np.float32)
    mask = np.array([r["covered"] for r in mapping])
    for j in np.flatnonzero(mask):
        rows = [index[t] for t in mapping[j]["target_ids"]]
        x[:,j] = np.log1p(counts[rows].sum(axis=0) / library * 1e6)
    return x, mask


def ranks(kernel, model, state, x, mask):
    _, z, valid = kernel.transform_external(x, mask, state)
    probabilities, _, valid, reasons = kernel.score(model, z, valid)
    ordered = [kernel.endpoint_row(p, model.classes_, [], ok, reason)["ranking_json"] for p,ok,reason in zip(probabilities,valid,reasons)]
    return ordered, probabilities, valid, reasons


def load_lower(key, lookup, store):
    if key not in store:
        folder = FINAL / "lower" / lookup[key]
        rules = js(folder / "RULES_AND_MEMBERS.json")
        with np.load(folder / "INPUT.npz") as z:
            store[key] = dict(candidates=z["regions"].tolist(),order=z["gene_rows"],ref=z["reference"],ann=rules["annotations"],folder=folder.name)
    return store[key]


def score_cohort(cohort, ids, features, panel, kernel, model, state):
    data, homology = raw_input(cohort)
    mapping = strict_map(features, ids, homology, data.index)
    x, mask = project(data, mapping)
    (OUT / (cohort + "_mapping.json")).write_text(json.dumps(mapping, separators=(",",":")), "utf8")
    head = pd.Index(features).get_indexer(panel)
    assert min(head) >= 0
    if cohort == "PsychENCODE":
        hx, hm = x[:,head], mask[head]
        admit = data.columns
    else:
        headmap = strict_map(panel, ids, homology, data.index)
        hx, hm = project(data, headmap)
        if not np.array_equal(hm, mask[head]):
            (OUT / (cohort + "_head_mapping.json")).write_text(json.dumps(headmap), "utf8")
        original = pd.read_csv(GEO / "prepared/FIXED_HEAD.csv", dtype=str)
        admit = original.loc[original.cohort.eq(cohort), "sample_id"]
    sample_ids = list(admit)
    ix = data.columns.get_indexer(sample_ids)
    assert min(ix) >= 0 and len(set(sample_ids)) == len(ix)
    ranking, prob, valid, reasons = ranks(kernel, model, state, hx[ix].copy(), hm)
    mean = np.load(FINAL / "SOURCE_RAW_MEAN.npy")
    x[:,~mask] = mean[~mask]
    if cohort != "PsychENCODE":
        sample = "10_10" if cohort == "GSE128537" else "LID_100660"
        q = data.columns.get_loc(sample)
        mapping_id = cohort + ("_Ensembl91" if cohort == "GSE128537" else "_Ensembl101") + "_source97_20261120"
        assert np.array_equal(hm, mask[head]), "Runtime example requires the same head and lower-axis correspondence"
        totals = np.ascontiguousarray(data.to_numpy(dtype=np.float64)).sum(axis=0)
        np.savez_compressed(OUT / (cohort + "_example.npz"),x=x[q,:,None],genes=features,observed=mask,
                            samples=np.array([sample]),full_library_counts=totals[q:q+1],contract_version=np.array("braintrace.three_tier_input.v1"),
                            mapping_version=np.array(mapping_id),normalization=np.array("full_library_ln1p_CPM"),full_count_rows=np.int32(len(data)))
    lookup = {tuple(b["networks"]):b["id"] for b in js(FINAL / "lower/BEAMS.json")}
    scoring = module(TOOL / "lower/scoring.py", "lower_score")
    store, records = {}, []
    for sid,q,nr,values,ok,reason in zip(sample_ids,ix,ranking,prob,valid,reasons):
        row = dict(cohort=cohort,sample_id=sid,Network_valid=bool(ok),Network_failure=str(reason),
                   Network_ranking_json=nr,Network_scores_json=json.dumps(values.tolist()),Group_valid=False,Exact_valid=False,
                   Group_members_json="[]",Group_ranking_json="[]",Exact_ranking_json="[]",Group_failure="network_unavailable",Exact_failure="network_unavailable")
        if ok:
            ref = load_lower(tuple(sorted(json.loads(nr)[:3])),lookup,store)
            row.update(scoring.lower_predict(ref,x[q],features,"",observed_mask=mask))
            row["beam_id"] = ref["folder"]
        records.append(row)
    pd.DataFrame(records).to_csv(OUT / (cohort + "_predictions.csv"),index=False)
    return dict(cohort=cohort,count_rows=len(data),input_samples=data.shape[1],predicted=len(ix),head_observed=int(hm.sum()),lower_observed=int(mask.sum()),head_axis=len(panel),lower_axis=len(features))


def score_ahba(ids, panel, kernel, model, state):
    sys.path.insert(0,str(SOURCE))
    from research.source_final_v1 import coverage_audit as coverage
    from research.masked_inference_v1.external_eval import _stream_observed_counts
    root = RUNS / "c5000l_source_final_20261016T140049Z/coverage"
    roster = js(root / "COVERAGE_INPUTS_FROZEN.json")
    with Path(roster["orthology"]["path"]).open(encoding="utf-8-sig") as f:
        resolver = coverage.Resolver(ids,csv.DictReader(f,delimiter="\t"))
    rows, stats = [], []
    for spec in coverage.prepare_ahba_specs(roster["ahba_counts"]):
        axis,_ = coverage.read_axis(Path(spec["path"]),delimiter=",",zip_member="RNAseqCounts.csv",sample_ids=spec["sample_ids"])
        mapping = resolver.platform(panel,axis)
        mask = np.array([not m["reason"] for m in mapping])
        names,x,_,errors = _stream_observed_counts(coverage,spec,mapping)
        ranking,prob,valid,reasons = ranks(kernel,model,state,x,mask)
        for sid,nr,values,ok,reason,error in zip(names,ranking,prob,valid,reasons,errors):
            rows.append(dict(cohort="AHBA",sample_id=sid,donor=spec["donor"],Network_valid=bool(ok and not error),
                             Network_failure=error or str(reason),Network_ranking_json=nr if not error else "[]",Network_scores_json=json.dumps(values.tolist())))
        stats.append(dict(cohort="AHBA/"+spec["donor"],count_rows=len(axis),input_samples=len(spec["sample_ids"]),predicted=len(names),head_observed=int(mask.sum()),head_axis=len(panel)))
    pd.DataFrame(rows).to_csv(OUT/"AHBA_predictions.csv",index=False)
    return stats


def main():
    OUT.mkdir(exist_ok=True)
    features = np.array(js(FINAL / "eligible_genes.json"))
    panel = js(FINAL / "network/genes.json")
    ids = source_ids()
    kernel = module(TOOL / "macaque/inference.py","network_kernel")
    model = joblib.load(FINAL / "network/model.joblib")
    with np.load(FINAL / "network/state.npz") as state, threadpool_limits(limits=1):
        coverage = []
        for cohort in ("GSE128537","GSE179328","PsychENCODE"):
            item = score_cohort(cohort,ids,features,panel,kernel,model,state)
            coverage.append(item)
            print(item,flush=True)
        coverage.extend(score_ahba(ids,panel,kernel,model,state))
    pd.DataFrame(coverage).to_csv(OUT / "coverage.csv",index=False)


if __name__ == "__main__":
    main()
