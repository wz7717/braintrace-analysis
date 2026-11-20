"""BrainTrace numerical helpers used by the Ensembl 97 reproduction route."""
from __future__ import annotations
import json
import numpy as np
import pandas as pd
from core.models import trace_corr
from core.scientific_validation import signal_reason
from core.signal_evidence import ranking_is_undifferentiated
from scripts.run_bo2023_leave_one_monkey_out_validation import build_region_reference, zscore
from scripts.run_bo2023_resolution_tier_validation import build_resolution_groups

def rank_metrics(ranking,truth):
    pos=ranking.index(truth)+1 if truth in ranking else 0
    return int(pos==1),int(0<pos<=3),1/pos if pos else 0.


def lower_predict(b,xq,genes,truth,observed_mask=None):
    cand=b['candidates'];order=b['order'];ref=b['ref'];ann=b['ann']
    q=xq[order];observed=pd.Series(xq,index=genes) if observed_mask is None else pd.Series(xq[observed_mask],index=genes[observed_mask])
    a=trace_corr(ref[:50],q[:50]);c=trace_corr(ref[:100],q[:100]);g=trace_corr(ref,q)
    fused=.25*zscore(a)+.75*zscore(c)
    ereason=signal_reason(genes[order[:50]].tolist(),observed,q[:50],a)
    if not ereason:ereason=signal_reason(genes[order[:100]].tolist(),observed,q[:100],c)
    if not ereason and ranking_is_undifferentiated(fused):ereason='exact_fusion_undifferentiated'
    greason=signal_reason(genes[order].tolist(),observed,q,g)
    # Preserve production ordering: rows first in reverse argsort of exact fusion,
    # then stable descending Group Top200 score and first occurrence of each group.
    initial=np.argsort(fused)[::-1].tolist() if not ereason else list(range(len(cand)))
    exact=[cand[j] for j in initial] if not ereason else []
    groups=[];members=[];gscores=[]
    if not greason:
        for j in sorted(initial,key=lambda j:float(g[j]),reverse=True):
            name=ann[cand[j]]['resolution_group']
            if name not in groups:
                groups.append(name);members.append(ann[cand[j]]['group_members']);gscores.append(float(g[j]))
        if len(groups)>1 and ranking_is_undifferentiated(np.array(gscores)):
            greason='group_aggregated_scores_undifferentiated';groups=[];members=[];gscores=[]
    truth_group=next((groups[j] for j,m in enumerate(members) if truth in m),'')
    eh1,eh3,emrr=rank_metrics(exact,truth);gh1,gh3,gmrr=rank_metrics(groups,truth_group)
    return {'candidate_regions_json':json.dumps(cand),'n_candidate_regions':len(cand),
            'n_candidate_groups':len({ann[r]['resolution_group'] for r in cand}),
            'Group_valid':not bool(greason),'Group_failure':greason or '',
            'Exact_valid':not bool(ereason),'Exact_failure':ereason or '',
            'Group_ranking_json':json.dumps(groups),'Group_members_json':json.dumps(members),'Group_scores_json':json.dumps(gscores),
            'Exact_ranking_json':json.dumps(exact),'exact_fusion_candidate_order_json':json.dumps(fused.tolist()),
            'corr50_candidate_order_json':json.dumps(a.tolist()),'corr100_candidate_order_json':json.dumps(c.tolist()),
            'corr200_candidate_order_json':json.dumps(g.tolist()),
            'Group_hit1':gh1,'Group_hit3':gh3,'Group_MRR':gmrr,'Exact_hit1':eh1,'Exact_hit3':eh3,'Exact_MRR':emrr,
            'group_top1_member_count':len(members[0]) if members else 0,
            'group_top3_union_member_count':len(set().union(*map(set,members[:3]))) if members else 0,
            'true_group_members_json':json.dumps(ann[truth]['group_members']) if truth in ann else '[]',
            'exact_score_tie_count':len(fused)-len(set(fused)),
            'group_region_score_tie_count':len(g)-len(set(g)),
            'beam_artifact':b['folder']}
