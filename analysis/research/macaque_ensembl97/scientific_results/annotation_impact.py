"""Measure the impact against the immediately preceding fold-local source edition."""
from pathlib import Path
import json
import os

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
NEW = HERE.parent/"source_analysis/results"
OLD = Path(os.environ["BRAINTRACE_PREVIOUS_SOURCE_RESULTS"])


def read(path):
    return json.loads(path.read_text(encoding="utf8"))


def main():
    prior,new = [pd.read_csv(f/"source_predictions.csv", dtype={"sample_id":str,"donor":str}).set_index("sample_id").sort_index() for f in (OLD,NEW)]
    assert new.index.equals(prior.index)
    assert new[["donor","truth_network","truth_region"]].equals(prior[["donor","truth_network","truth_region"]])
    out = new[["donor","truth_region"]].copy()
    summary = {}
    for tier in ["Network","Group","Exact"]:
        col = "Group_members_json" if tier=="Group" else tier+"_ranking_json"
        def ordered(value):
            arr = json.loads(value)
            return [tuple(sorted(x)) for x in arr] if tier=="Group" else arr
        pairs = [(ordered(a),ordered(b)) for a,b in zip(prior[col],new[col])]
        for name,flags in [("all_ranks",[a!=b for a,b in pairs]),("top1",[a[:1]!=b[:1] for a,b in pairs]),("top3_set",[set(a[:3])!=set(b[:3]) for a,b in pairs])]:
            column = tier+"_"+name+"_changed"
            out[column] = flags
            summary[column] = sum(flags)
        for k in (1,3):
            col = tier+"_hit"+str(k)
            out[col+"_changed"] = prior[col] != new[col]
            summary[col+"_changed"] = int(out[col+"_changed"].sum())
    out.to_csv(HERE/"ANNOTATION_IMPACT_PER_SPECIMEN.csv")
    references = []
    for rel in ["final"]+[str(f.relative_to(NEW)) for f in sorted((NEW/"folds").iterdir())]:
        a,b = OLD/rel,NEW/rel
        prev,newp = [read(y/"network/genes.json") for y in (a,b)]
        oa,na = [set(read(y/"eligible_genes.json")) for y in (a,b)]
        oldbeams,newbeams = [{tuple(q["networks"]):y/"lower"/q["id"] for q in read(y/"lower/BEAMS.json")} for y in (a,b)]
        shared = oldbeams.keys() & newbeams.keys()
        panels,refs,partitions = 0,0,0
        for key in shared:
            aa,bb = oldbeams[key],newbeams[key]
            with np.load(aa/"INPUT.npz") as x, np.load(bb/"INPUT.npz") as y:
                panels += int(not np.array_equal(x["genes"],y["genes"]))
                refs += int(not np.array_equal(x["reference"],y["reference"]))
            ar,br = [read(v/"RULES_AND_MEMBERS.json")["annotations"] for v in (aa,bb)]
            partitions += int(any(ar[k]["resolution_group"]!=br[k]["resolution_group"] for k in ar))
        references.append(dict(fit=rel,old_eligible=len(oa),new_eligible=len(na),gained_keys=len(na-oa),lost_keys=len(oa-na),
                               head_set_replacements=len(set(newp)-set(prev)),head_positions_changed=sum(x!=y for x,y in zip(prev,newp)),
                               shared_gates=len(shared),old_only_gates=len(oldbeams)-len(shared),new_only_gates=len(newbeams)-len(shared),
                               changed_fine_panel_order=panels,changed_centroid_arrays=refs,changed_group_assignments=partitions))
    pd.DataFrame(references).to_csv(HERE/"ANNOTATION_IMPACT_REFERENCES.csv",index=False)
    (HERE/"ANNOTATION_IMPACT.json").write_text(json.dumps(summary,indent=2)+"\n",encoding="utf8")
    print(summary)


if __name__ == "__main__":
    main()
