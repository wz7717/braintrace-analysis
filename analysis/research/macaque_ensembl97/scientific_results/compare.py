"""Pure set/rank comparison of saved scores; never import model/inference code."""
import math

METHODS=('Group','FlatPearson200','Exact50_100')
CONTROLS=METHODS[1:]

def ordered_union(groups):
    return list(dict.fromkeys(r for group in groups for r in group))

def flat_ranking(record):
    candidates=record['candidates'];scores=record['corr200']
    assert len(candidates)==len(scores) and all(math.isfinite(v) for v in scores)
    score_by_region=dict(zip(candidates,scores))
    initial=record['exact_rank'] if record['exact_valid'] else candidates
    assert len(initial)==len(candidates) and set(initial)==set(candidates)
    return sorted(initial,key=score_by_region.__getitem__,reverse=True)

def matched_sets(record,k):
    assert k in (1,3)
    if not record['group_valid']:
        return {m:[] for m in METHODS},0,False,0
    group=ordered_union(record['groups'][:k]);m=len(group)
    assert m>0 and set(group).issubset(record['candidates'])
    flat=flat_ranking(record)
    exact=record['exact_rank'][:m] if record['exact_valid'] else []
    assert len(flat)>=m and (not record['exact_valid'] or len(exact)==m)
    tied=0
    if m<len(flat):
        scores=dict(zip(record['candidates'],record['corr200']))
        if scores[flat[m-1]]==scores[flat[m]]:
            tied=sum(v==scores[flat[m-1]] for v in scores.values())
    outputs={'Group':group,'FlatPearson200':flat[:m],'Exact50_100':exact}
    assert len(outputs['Group'])==len(outputs['FlatPearson200'])==m
    assert all(len(v)==len(set(v)) for v in outputs.values())
    return outputs,m,True,tied

def hit(values,allowed,evaluable=True):
    return int(bool(set(values)&set(allowed))) if evaluable else None
