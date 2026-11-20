"""Summarize the revised 819 source holdouts against the named historical endpoints.

Run from the BrainTrace project root after source_analysis has exported source_predictions.csv.
Pass --historical-check to exercise the same calculations on the original source instead.
"""
from pathlib import Path
import argparse
import os
import io
import json
import subprocess
import sys
import tempfile
import zipfile

import numpy as np
import pandas as pd

PROJECT = Path(os.environ.get("BRAINTRACE_STUDY_ROOT", Path.cwd())).resolve()
OUT = Path(__file__).resolve().parent
RUNS = PROJECT / 'BrainTrace_algorithm_research_20261015/research/runs'
BRIDGE = RUNS / 'hierarchical_bridge_20261101_01'
BUDGET = RUNS / 'group_equal_candidate_budget_20261102_01'
REGIONAL = RUNS / 'regional_balance_psychmask_20261102_01/stage1_region_balance'
COORDINATES = ('Network','Group','Exact')
from compare import matched_sets

def mean_ci(by_donor, draws):
    values=by_donor.to_numpy(dtype=float)
    return [float(values.mean()), *np.quantile(values[draws].mean(axis=1), [.025,.975]).tolist()]

def read_table(path):
    return pd.read_csv(path, dtype={'sample_id':str,'donor':str}, keep_default_na=False, float_precision='round_trip')

def summaries(new, historical, draw):
    rows=[]
    for label,group in [('corrected',new),('archived_vsd',historical['original_three_tier']),('archived_common_universe',historical['floor_clip_hybrid'])]:
        for tier in COORDINATES:
            row={'version':label,'tier':tier,'samples':len(group),'donors':group.donor.nunique()}
            for metric,field in [('Top1','hit1'),('Top3','hit3'),('MRR','MRR')]:
                key=f'{tier}_{field}'
                vals=group.groupby('donor')[key].mean()
                row[metric],row[metric+'_CI95_low'],row[metric+'_CI95_high']=mean_ci(vals,draw)
                if field!='MRR':
                    row[metric+'_hits']=int(group[key].sum())
            rows.append(row)
    pairs=[]
    for label,model in [('corrected',new),('archived_common_universe',historical['floor_clip_hybrid'])]:
        vsd=historical['original_three_tier'].loc[model.index]
        for tier in COORDINATES:
            for metric,field in [('Top1','hit1'),('Top3','hit3'),('MRR','MRR')]:
                key=f'{tier}_{field}'
                delta=(model[key]-vsd[key]).groupby(model.donor).mean()*100
                mean,lo,hi=mean_ci(delta,draw) if draw is not None else (float(delta.mean()),np.nan,np.nan)
                pairs.append(dict(version=label,comparator='archived_vsd',tier=tier,metric=metric,difference_pp=mean,CI95_low_pp=lo,CI95_high_pp=hi))
    return pd.DataFrame(rows),pd.DataFrame(pairs)

def enrich(new):
    enriched=[]
    for sid,r in new.iterrows():
        truth=r.truth_region
        pool=json.loads(r.candidate_regions_json)
        group_lists=json.loads(r.Group_members_json)
        exact=json.loads(r.Exact_ranking_json)
        in_pool=truth in pool
        supported=bool(r.fine_supported)
        rec=dict(candidates=pool, groups=group_lists, exact_rank=exact, corr200=json.loads(r.corr200_candidate_order_json),
                 group_valid=bool(r.Group_valid),exact_valid=bool(r.Exact_valid))
        for k in [1,3]:
            lists,m,budget_ok,ties=matched_sets(rec,k)
            e=dict(sample_id=sid,donor=r.donor,truth=truth,fine_supported=supported,in_pool=in_pool,
                K=len(pool),a=int(in_pool),k=k,m=m,budget_ok=budget_ok,flat_ties=ties,
                exact_valid=bool(r.Exact_valid),group_valid=bool(r.Group_valid))
            for meth,vals in lists.items():
                e[meth]=int(truth in vals)
            e['exact_observed']=int(truth in exact[:k])
            e['group_expected']=(m/len(pool)) if in_pool else 0.0
            e['exact_expected']=(min(k,len(exact))/len(pool)) if in_pool else 0.0
            e['truth_percentile']=100*exact.index(truth)/(len(pool)-1) if in_pool and len(pool)>1 else np.nan
            enriched.append(e)
    return pd.DataFrame(enriched)

def source_subsets(data):
    return [('all_samples',data),('supported_812',data[data.fine_supported])]

def retrieval_tables(data, sets=None):
    matched=[];attribution=[];random=[];rank_stats=[]
    if sets is None:
        sets = [(i, key, frame) for i,(key,frame) in enumerate(source_subsets(data))]
    for idx,endpoint,sub in sets:
        draw=np.load(BUDGET/f'analysis_attempt01/DONOR_DRAWS_{idx:02d}.npy') if idx < 8 else None
        for k,g in sub.groupby('k'):
            by=g.groupby('donor')
            donor_values=by[['Group','FlatPearson200','Exact50_100']].mean()
            row=dict(endpoint=endpoint,k=k,n=len(g),donors=g.donor.nunique(),mean_m=float(by.m.mean().mean()),
                     min_m=g.m.min(),max_m=g.m.max())
            for method in donor_values:
                row[method+'_pct']=100*donor_values[method].mean()
            for control in ['FlatPearson200','Exact50_100']:
                delta=(donor_values.Group-donor_values[control])*100
                row[control+'_delta_pp'],row[control+'_low_pp'],row[control+'_high_pp']=mean_ci(delta,draw) if draw is not None else (float(delta.mean()),np.nan,np.nan)
                row['Group_only_vs_'+control]=int(((g.Group==1)&(g[control]==0)).sum())
                row[control+'_only_vs_Group']=int(((g.Group==0)&(g[control]==1)).sum())
            matched.append(row)
            val=dict(endpoint=endpoint,k=k,n=len(g),donors=g.donor.nunique(),K_median=g.K.median(),K_min=g.K.min(),K_max=g.K.max(),
                     a_min=g.a.min(),a_max=g.a.max(),in_pool_count=g.in_pool.sum(),in_pool_pct=100*by.in_pool.mean().mean(),mean_m=row['mean_m'])
            for tier,ob,ex in [('Exact','exact_observed','exact_expected'),('Group','Group','group_expected')]:
                val[tier+'_observed_pct']=100*by[ob].mean().mean()
                val[tier+'_random_pct']=100*by[ex].mean().mean()
                val[tier+'_difference_pp']=val[tier+'_observed_pct']-val[tier+'_random_pct']
                eligible=g[g.in_pool]
                is_valid=eligible.exact_valid if tier=='Exact' else eligible.group_valid
                hits=int(g[ob].sum())
                attribution.append(dict(endpoint=endpoint,tier=tier,k=k,n=len(g),donors=g.donor.nunique(),
                     unsupported=int((~g.fine_supported).sum()),supported_truth_outside_pool=int((g.fine_supported&(~g.in_pool)).sum()),
                     unavailable=int((~is_valid).sum()),in_pool_rank_miss=int(((eligible[ob]==0)&is_valid).sum()),
                     hits=hits,in_pool_n=len(eligible),in_pool_donors=eligible.donor.nunique(),
                     in_pool_donor_hit_percent=100*eligible.groupby('donor')[ob].mean().mean()))
            random.append(val)
            if k==1:
                ranks=g.truth_percentile.dropna()
                rank_stats.append(dict(endpoint=endpoint,percentile_n=len(ranks),mean_percentile=ranks.mean(),median_percentile=ranks.median(),
                    q25=ranks.quantile(.25),q75=ranks.quantile(.75),eligible_donor_count=g.loc[g.truth_percentile.notna(),'donor'].nunique()))
    return [pd.DataFrame(a) for a in (matched,attribution,random,rank_stats)]

def regional_tables(new):
    saved=read_table(REGIONAL/'ALL_WEIGHTED_RECORDS.csv')
    metrics=('Exact_Top1','Exact_Top3','Group_Top1','Group_Top3')
    result=[]
    for cohort,cohort_table in saved.groupby('comparison',sort=False):
        source=cohort_table[cohort_table.condition=='source_complete'].copy()
        external=cohort_table[cohort_table.condition=='external']
        for metric in metrics:
            corrected=new.loc[source.sample_id,metric.replace('Top','hit')].to_numpy()
            source[metric+'_corrected']=corrected
            measure=100*float(corrected@source.weight)
            old=100*float(source[metric]@source.weight)
            ext=100*float(external[metric]@external.weight)
            deletions=[]
            for donor in source.donor.unique():
                omit=source[source.donor!=donor]
                if omit.label.nunique()==source.label.nunique():
                    rates=omit.groupby(['donor','label'])[metric+'_corrected'].mean().groupby('label').mean()
                    deletions.append(100*rates.mean())
            result.append(dict(comparison=cohort,metric=metric,labels=source.label.nunique(),samples=len(source),donors=source.donor.nunique(),
                   revised_source_pct=measure,archived_unmasked_source_pct=old,external_pct=ext,
                   external_minus_revised_source_pp=ext-measure,omissions=len(deletions),omission_min_pct=min(deletions),omission_max_pct=max(deletions)))
    return pd.DataFrame(result)

def read_archive_script(z,suffix):
    key=next(s for s in z.namelist() if s.endswith('/'+suffix))
    return z.read(key)

def figures(summary,matched,dest):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from records import SETS
    dest.mkdir(exist_ok=True,parents=True)
    archive=next(PROJECT.glob('BrainTrace_R8_*/*R9_Archive*/published_files/*.zip'))
    with zipfile.ZipFile(archive) as z:
        source=pd.read_csv(io.BytesIO(read_archive_script(z,'SOURCE_THREE_TIER.csv')))
        for i,row in source.iterrows():
            if row.arm=='floor_clip_hybrid':
                s=summary[(summary.version=='corrected')&(summary.tier==row.tier)].iloc[0]
                for key in ('Top1','Top3','Top1_CI95_low','Top1_CI95_high','Top3_CI95_low','Top3_CI95_high'):
                    source.at[i,key]=s[key]
        source.to_csv(dest/'SOURCE_THREE_TIER.csv',index=False)
        external=pd.read_csv(io.BytesIO(read_archive_script(z,'EXTERNAL_THREE_TIER.csv')))
        samples=pd.read_csv(OUT.parent/'external_analysis/results/ENDPOINTS_PER_SPECIMEN.csv')
        draws_dir=RUNS/'macaque_external_comparators_20261101_01/analysis'
        order=json.loads((draws_dir/'BOOTSTRAP_DRAW_IDENTITIES.json').read_text(encoding='utf8'))
        for i,row in external.iterrows():
            records=samples[(samples.cohort==row.cohort)&(samples.subset==row.evaluation_set)&(samples.tier==row.tier)]
            donor_list=order[row.cohort]['donor_order']
            draws=np.load(draws_dir/(row.cohort+'_DONOR_DRAWS.npy'))
            for k in (1,3):
                values=records.groupby('donor')['corrected_hit'+str(k)].mean().loc[donor_list].to_numpy()
                if not np.isclose(row['Top'+str(k)],values.mean()) or (records['original_hit'+str(k)]!=records['corrected_hit'+str(k)]).any():
                    low,high=np.quantile(values[draws].mean(axis=1),[.025,.975])
                    external.loc[i,['Top'+str(k),'Top'+str(k)+'_CI95_low','Top'+str(k)+'_CI95_high']]=[values.mean(),low,high]
        external.to_csv(dest/'EXTERNAL_THREE_TIER.csv',index=False)
        with tempfile.TemporaryDirectory() as tmp:
            tool=Path(tmp)
            for name in ('plot_common.py','plot_figure_7.py'):
                (tool/name).write_bytes(read_archive_script(z,name))
            subprocess.run([sys.executable,str(tool/'plot_figure_7.py'),'--inputs',str(dest),'--output',str(dest/'Figure_7.png')],check=True)
    pairs=pd.read_csv(BUDGET/'analysis_attempt01/PAIRED_COMPARISONS.csv')
    ext=pd.read_csv(OUT.parent/'external_analysis/results/RETRIEVAL_MATCHED.csv')
    plot_rows=pd.concat([matched.assign(cohort='Source_OOF'),ext],ignore_index=True)
    for row in plot_rows.itertuples():
        for control in ('FlatPearson200','Exact50_100'):
            flag=(pairs.cohort==row.cohort)&(pairs.evaluation_set==row.endpoint)&(pairs.k_groups==row.k)&(pairs.control==control)
            assert flag.sum()==1
            for source_name,dest_name in [('delta_pp','Group_minus_control_pp'),('low_pp','ci95_low_pp'),('high_pp','ci95_high_pp')]:
                pairs.loc[flag,dest_name]=getattr(row,control+'_'+source_name)
            pairs.loc[flag,'Group_hit_percent']=row.Group_pct
            pairs.loc[flag,'control_hit_percent']=getattr(row,control+'_pct')
            pairs.loc[flag,'pooled_Group_only_hits']=getattr(row,'Group_only_vs_'+control)
            pairs.loc[flag,'pooled_control_only_hits']=getattr(row,control+'_only_vs_Group')
    keys=['cohort','evaluation_set','k_groups','control','samples','donors','Group_minus_control_pp','ci95_low_pp','ci95_high_pp']
    pairs[keys].to_csv(dest/'FIGURE_S10_POINTS.csv',index=False)
    # Run just the original display block against the revised coordinates, without its report side effects.
    text=(BUDGET/'report.py').read_text(encoding='utf8')
    block=text[text.index('plt.rcParams.update'):text.index('for ext in (')]
    scope=dict(plt=plt,pd=pd,np=np,pairs=pairs,SETS=SETS,Line2D=Line2D)
    exec(compile(block,'archived_S10_plot_block','exec'),scope)
    fig=scope['fig']
    for extension in ('png','pdf'):
        fig.savefig(dest/('Figure_S10.'+extension),dpi=190,facecolor='white')
    plt.close(fig)

def run():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--historical-check',action='store_true')
    parser.add_argument('--figures', action='store_true', help='Also recreate figure input tables using the original study archive and resampling draws.')
    args=parser.parse_args()
    orig=read_table(BRIDGE/'execution/PREDICTIONS.csv')
    old={s:g.set_index('sample_id').sort_index() for s,g in orig.groupby('arm')}
    if args.historical_check:
        new=old['floor_clip_hybrid']
    else:
        new=read_table(OUT.parent/'source_analysis/results/source_predictions.csv').set_index('sample_id').sort_index()
    assert len(new)==819 and new.index.equals(old['floor_clip_hybrid'].index)
    # Make axis identity explicit before any paired calculation.
    assert (new[['donor','truth_region','truth_network']]==old['floor_clip_hybrid'][['donor','truth_region','truth_network']]).all().all()
    draw=np.load(BRIDGE/'analysis/DONOR_DRAWS.npy')
    summary,pairs=summaries(new,old,draw)
    enriched=enrich(new)
    matched,attribution,random,rank_stats=retrieval_tables(enriched)
    regional=regional_tables(new)
    if args.historical_check:
        print(summary[['version','tier','Top1','Top3']].to_string(index=False))
        print(matched.iloc[0].to_dict())
        print(random.iloc[0].to_dict())
        return
    metrics=['Network_hit1','Network_hit3','Group_hit1','Group_hit3','Exact_hit1','Exact_hit3','Group_MRR','Exact_MRR']
    before=old['floor_clip_hybrid'].loc[new.index]
    donors=new[['donor','truth_region']+metrics].copy()
    for m in metrics:
        donors['change_'+m]=new[m]-before[m]
    donors.groupby('donor')[metrics+['change_'+m for m in metrics]].mean().reset_index().to_csv(OUT/'SOURCE_BY_DONOR.csv',index=False)
    points=new[['donor','truth_region']].copy()
    for m in metrics:
        points['old_'+m]=before[m]
        points['new_'+m]=new[m]
    changed=np.any(new[metrics].to_numpy()!=before[metrics].to_numpy(),axis=1)
    points.loc[changed].reset_index().to_csv(OUT/'SOURCE_CHANGED_ENDPOINTS.csv',index=False)
    for name,df in [('SOURCE_TIER_RATES',summary),('VSD_PAIRED',pairs),('TABLE_7_S38',matched),('TABLE_S41',attribution),
                    ('TABLE_9_S45',random),('WITHIN_POOL_PERCENTILE',rank_stats),('TABLE_S42_UNMASKED',regional)]:
        df.to_csv(OUT/(name+'.csv'),index=False)
    if args.figures:
        figures(summary,matched,OUT/'figures')
    print(summary[summary.version=='corrected'][['tier','Top1','Top3']].to_string(index=False))
    print('Source summaries and both figures saved.')

if __name__=='__main__':
    run()
