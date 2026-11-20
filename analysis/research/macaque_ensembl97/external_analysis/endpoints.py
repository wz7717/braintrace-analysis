"""Compare complete predictions and frozen evaluation sets at both reference versions."""
import json
import pandas as pd
from transfer import GEO, OUT, PROJECT, PSYCH, RUNS, js


def groups(row, tier):
    key = "Group_members_json" if tier == "Group" else tier + "_ranking_json"
    values = json.loads(row[key])
    return values if tier == "Group" else [[x] for x in values]


def load(cohort):
    new = pd.read_csv(OUT / (cohort+"_predictions.csv"), dtype={"sample_id":str}, keep_default_na=False).set_index("sample_id")
    if cohort == "AHBA":
        f = RUNS / "floor_clip_ahba_development_20261017T124707Z/external/AHBA_ENDPOINT_ROWS.csv"
        old = pd.read_csv(f).query("route_id=='floor_clip'").set_index("profile_id")
        old = old.rename(columns={"ranking_json":"Network_ranking_json"})
        truth = pd.read_csv(RUNS/"c5000l_source_final_20261016T140049Z/coverage/AHBA_TRUTH_ONLY_LEDGER.csv").set_index("profile_id")
        cases = []
        for sid,r in truth.iterrows():
            allowed = json.loads(r.allowed_truth_json)
            if allowed:
                for key in ["allowed", "singleton" if len(allowed)==1 else "multiple"]:
                    cases.append((sid,r.donor_id,key,"Network",allowed))
        return new,old,cases
    if cohort == "PsychENCODE":
        old = pd.DataFrame(js(PROJECT/"BrainTrace_PsychENCODE_Adult_Lower_Completion_20261101/ALL_77_PREDICTIONS.json")).set_index("sample_id")
        for tier,orig in [("Network_ranking_json","network_ranking"),("Group_members_json","Group_members"),("Exact_ranking_json","Exact_ranking")]:
            old[tier] = old[orig].map(json.dumps)
        truth = pd.read_csv(PSYCH/"FROZEN_SAMPLE_LABELS.csv",keep_default_na=False)
        cases = []
        for r in truth.itertuples():
            for subset,ok in [("primary_coarse_34",r.network_primary),("extended_coarse_68",r.network_extended)]:
                if ok:
                    cases.append((r.sample_id,r.donor,subset,"Network",json.loads(r.network_allowed_json)))
            if r.lower_set != "not_scored":
                for tier in ["Group","Exact"]:
                    cases.append((r.sample_id,r.donor,r.lower_set,tier,json.loads(r.lower_allowed_json)))
        return new,old,cases
    head = pd.read_csv(GEO/"prepared/FIXED_HEAD.csv",dtype={"sample_id":str})
    head = head[head.cohort.eq(cohort)].set_index("sample_id")
    prior = pd.read_csv(GEO/"external/ALL_PREDICTIONS.csv",dtype={"sample_id":str},keep_default_na=False)
    old = prior[prior.cohort.eq(cohort) & prior.arm.eq("same_head_original_lower")].set_index("sample_id")
    old = old.loc[head.index].copy()
    old["Network_ranking_json"] = head.ranking_json
    cases = [(sid,r.donor,r.endpoint,"Network",json.loads(r.allowed_networks_json)) for sid,r in head.iterrows()]
    admit = pd.read_csv(GEO/"prepared/ANALYSIS_ADMISSION.csv",dtype={"sample_id":str})
    for r in admit[admit.cohort.eq(cohort)].itertuples():
        if r.fine_endpoint in ("primary_unique","secondary_ambiguous","coarsened_parent_only"):
            for tier in ("Group","Exact"):
                cases.append((r.sample_id,r.donor,r.fine_endpoint,tier,json.loads(r.allowed_regions_json)))
    return new,old,cases


def main():
    points, summary, changes, sample_changes = [], [], [], []
    for cohort in ("GSE128537","GSE179328","PsychENCODE","AHBA"):
        new,old,cases = load(cohort)
        assert set(new.index) == set(old.index)
        tiers = ["Network"] if cohort == "AHBA" else ["Network","Group","Exact"]
        for tier in tiers:
            ranks = {}
            for sid,r in new.iterrows():
                a,b = groups(old.loc[sid],tier),groups(r,tier)
                if tier=="Group":
                    a,b = [sorted(g) for g in a], [sorted(g) for g in b]
                ranks[sid] = a,b
                sample_changes.append(dict(cohort=cohort,sample_id=sid,tier=tier,ordering_changed=a!=b,top1_changed=a[:1]!=b[:1],top3_set_changed={tuple(g) for g in a[:3]}!={tuple(g) for g in b[:3]}))
            c = pd.DataFrame(sample_changes[-len(ranks):])
            changes.append(dict(cohort=cohort,tier=tier,profiles=len(ranks),**{k:int(c[k].sum()) for k in ("ordering_changed","top1_changed","top3_set_changed")}))
            for sid,donor,subset,layer,allowed in cases:
                if layer != tier:
                    continue
                p = dict(cohort=cohort,subset=subset,tier=tier,sample_id=sid,donor=donor)
                for version,ranking in zip(("original","corrected"),ranks.get(sid,([],[]))):
                    rank = next((i for i,g in enumerate(ranking,1) if set(g).intersection(allowed)),0)
                    p[version+"_available"] = int(bool(ranking))
                    p[version+"_first_size"] = len(ranking[0]) if ranking else 0
                    p[version+"_top3_size"] = len({v for group in ranking[:3] for v in group})
                    p[version+"_hit1"] = int(rank==1)
                    p[version+"_hit3"] = int(0<rank<=3)
                points.append(p)
    point = pd.DataFrame(points)
    for (cohort,subset,tier),g in point.groupby(["cohort","subset","tier"],sort=False):
        row = dict(cohort=cohort,subset=subset,tier=tier,specimens=len(g),donors=g.donor.nunique())
        for version in ("original","corrected"):
            row[version+"_available"] = int(g[version+"_available"].sum())
            for size in ("first_size","top3_size"):
                row[version+"_"+size] = g.groupby("donor")[version+"_"+size].mean().mean()
            for k in (1,3):
                metric = version+"_hit"+str(k)
                row[metric] = int(g[metric].sum())
                row[metric+"_donor_pct"] = 100*g.groupby("donor")[metric].mean().mean()
        summary.append(row)
    for name,data in [("FROZEN_ENDPOINTS",pd.DataFrame(summary)),("RANK_CHANGES",pd.DataFrame(changes)),("ENDPOINTS_PER_SPECIMEN",point),("RANK_CHANGES_PER_SPECIMEN",pd.DataFrame(sample_changes))]:
        data.to_csv(OUT/(name+".csv"),index=False)
    print(pd.DataFrame(changes).to_string(index=False))


if __name__ == "__main__":
    main()
