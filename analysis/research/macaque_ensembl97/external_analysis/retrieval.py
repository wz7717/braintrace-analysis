"""Fine-region retrieval on projected external counts, using the original lists and budgets."""
from math import comb
import sys
import json
import numpy as np
import pandas as pd
from endpoints import load
from transfer import OUT, SOURCE

sys.path.insert(0,str(SOURCE.parent/"scientific_results"))
from source_tables import retrieval_tables, matched_sets
from records import SETS
from source_donor_details import failure_attribution, retrieval_contrasts


def chance(K, a, m):
    return 1 - (comb(K-a,m)/comb(K,m)) if K and a and m else 0.


def main():
    merged = [[],[],[],[]]
    diagnostics = [[],[],[],[]]
    for cohort in ("GSE128537","GSE179328","PsychENCODE"):
        new,_,cases = load(cohort)
        detail = []
        for sid,donor,subset,tier,allowed in cases:
            if tier != "Exact":
                continue
            r=new.loc[sid]
            pool=json.loads(r.candidate_regions_json)
            exact=json.loads(r.Exact_ranking_json)
            allowed=set(allowed)
            a=len(allowed.intersection(pool))
            K=len(pool)
            valid=bool(r.Exact_valid)
            group_valid=bool(r.Group_valid)
            rec=dict(candidates=pool,groups=json.loads(r.Group_members_json),exact_rank=exact,corr200=json.loads(r.corr200_candidate_order_json),
                     group_valid=group_valid,exact_valid=valid)
            for k in (1,3):
                lists,m,budget,ties=matched_sets(rec,k)
                values=dict(sample_id=sid,donor=donor,endpoint=subset,k=k,K=K,a=a,in_pool=bool(a),m=m,exact_valid=valid,group_valid=group_valid,fine_supported=True,
                            exact_observed=int(bool(allowed.intersection(exact[:k]))),exact_expected=chance(K,a,min(k,len(exact))),group_expected=chance(K,a,m),
                            truth_percentile=100*min(exact.index(x) for x in allowed.intersection(exact))/(K-1) if a and valid and K>1 else np.nan)
                for method,labels in lists.items():
                    values[method]=int(bool(allowed.intersection(labels)))
                detail.append(values)
        data=pd.DataFrame(detail)
        sets=[(i,endpoint,data[data.endpoint.eq(endpoint)]) for i,(c,endpoint) in enumerate(SETS) if c==cohort]
        grouped=[(name,g) for _,name,g in sets]
        detail_frames=(*failure_attribution(data,cohort,grouped),*retrieval_contrasts(data,cohort,grouped))
        for target,df in zip(diagnostics,detail_frames):
            target.append(df)
        for target,df in zip(merged,retrieval_tables(data,sets=sets)):
            target.append(df.assign(cohort=cohort))
    for name,frames in zip(("S41B_EXTERNAL_BY_DONOR","S41B_EXTERNAL_CATEGORY_SHARES","S38_EXTERNAL_PAIRED_BY_DONOR","S38_EXTERNAL_DONOR_OMISSION"),diagnostics):
        pd.concat(frames,ignore_index=True).to_csv(OUT/(name+".csv"),index=False)
    for name,frames in zip(("MATCHED","ATTRIBUTION","CHANCE","RANKS"),merged):
        pd.concat(frames,ignore_index=True).to_csv(OUT/("RETRIEVAL_"+name+".csv"),index=False)


if __name__=="__main__":
    main()
