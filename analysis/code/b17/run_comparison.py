"""Budgeted final-binding comparison; real mode requires the delegated review lock."""
import argparse
import csv
import json
import os
from pathlib import Path
import subprocess
import sys
import time
import zipfile

from budget_control import ByteLedger, PRODUCTION_CAPS
from process_resources import process_tree_rss, terminate_tree

METHODS = ('PEARSON','SINGLER')
SCIENCE_PYTHON = sys.executable
THREAD_ENV = {name:'1' for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS','VECLIB_MAXIMUM_THREADS')}


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def run(args):
    if args.mode not in ('synthetic', 'real'):
        print(json.dumps({'status':'EXECUTION_MODE_DISABLED','real_data_execution_enabled':False})); return 3
    is_real = args.mode == 'real'
    prefix = 'REAL' if is_real else 'SYNTHETIC'
    execution_lock = None
    if is_real:
        execution_lock = read(args.execution_lock)
        if (execution_lock.get('supervision_checkpoint') != 'PASSED'
                or execution_lock.get('real_execution_authorized') is not True
                or execution_lock.get('real_data_execution_enabled') is not True):
            print(json.dumps({'status':'AWAITING_DELEGATED_SUPERVISION_BEFORE_NUMERIC_READ'})); return 3
        if args.synthetic_overrides:
            raise ValueError('Real execution uses frozen production caps without overrides')
    manifest = read(args.manifest)
    if manifest.get('synthetic') is not (not is_real):
        print(json.dumps({'status':'MANIFEST_DATA_MODE_MISMATCH_BEFORE_NUMERIC_READ'})); return 3
    if is_real and (manifest.get('real_execution_authorized') is not True
                    or manifest.get('real_data_execution_enabled') is not True):
        print(json.dumps({'status':'REAL_MANIFEST_EXECUTION_DISABLED_BEFORE_NUMERIC_READ'})); return 3
    root = Path(args.output).absolute()
    root.mkdir(parents=True, exist_ok=False)
    for name in ('cache','queries','jobs'):
        (root/name).mkdir()
    records_file = Path(args.manifest).parent / manifest['queries']['records']
    rows = list(csv.DictReader(records_file.open(encoding='utf-8-sig',newline='')))
    roster = root/'ROSTER.csv'
    with roster.open('w',encoding='utf-8',newline='') as stream:
        writer = csv.DictWriter(stream,fieldnames=['candidate_specimen_key','reported_donor_case_id','model_truth_id'])
        writer.writeheader()
        for row in rows:
            writer.writerow({'candidate_specimen_key':row.get('candidate_specimen_key') or row['sample_id'],
                             'reported_donor_case_id':row['donor_id'],'model_truth_id':row['truth_id']})
    caps = dict(PRODUCTION_CAPS)
    overrides = read(args.synthetic_overrides) if args.synthetic_overrides else {}
    for key,value in overrides.items():
        if key in caps:
            if value <= 0 or value > PRODUCTION_CAPS[key]:
                raise ValueError('Synthetic overrides may reduce, never raise frozen caps')
            caps[key] = value
    for mapping,limit in (('method_build_seconds',caps['build_seconds']),('method_query_seconds',caps['query_seconds'])):
        for method,value in overrides.get(mapping,{}).items():
            if method not in METHODS or value <= 0 or value > limit:
                raise ValueError('Invalid reduced synthetic method cap')
    ledger = ByteLedger(root, args.preparation_root if is_real else None)
    provenance = manifest.get('binding_provenance', {})
    payload = {'input_data_kind':'REAL_GSE25219_FROZEN_FINAL_PLM' if is_real else 'SYNTHETIC_B17_END_TO_END',
               'mode':args.mode,'synthetic':not is_real,'methods':list(METHODS),
               'source_binding_id':execution_lock.get('source_binding_id') if is_real else None,
               'query_binding_id':execution_lock.get('query_binding_id') if is_real else None,
               'source_binding_state':'B15_ACCEPTED_METADATA_LINKED_CURRENT_NUMERIC_FORMAL_ATTEMPT' if is_real else 'SYNTHETIC_GENERATED',
               'source_input_contract_status':'NOT_RUN','binding_provenance':provenance,
               'shared_input_contract_status':'NOT_RUN',
               'method_reference_states':{m:{'status':'NOT_RUN'} for m in METHODS},
               'records':[{'candidate_specimen_key':r.get('candidate_specimen_key') or r['sample_id'],
                           'sample_id':r['sample_id'],'GSM_sample_id':r['sample_id'],
                           'original_roster_index_1based':r.get('original_roster_index_1based',i+1),
                           'reported_donor_case_id':r['donor_id'],'original_truth_id':r['truth_id'],
                           'methods':{m:{'status':'NOT_RUN'} for m in METHODS}} for i,r in enumerate(rows)]}
    state = {'scope':'REAL_FINAL_EXTERNAL' if is_real else 'SYNTHETIC_ONLY','mode':args.mode,
             'real_data_execution_enabled':is_real,
             'real_G_reference_marker_or_prediction_authorized':is_real,'synthetic_execution_authorized':not is_real,
             'execution_lock':str(Path(args.execution_lock).absolute()) if is_real else None,
             'production_caps_unchanged':dict(PRODUCTION_CAPS),'effective_caps':caps,'synthetic_overrides':overrides,
             'manifest':str(Path(args.manifest).absolute()),'original_query_N':len(rows),
             'original_query_D':len({r['donor_id'] for r in rows}),
             'first_source_G_reference_or_prediction_monotonic':time.monotonic(),
             'first_source_G_reference_or_prediction_epoch':time.time(), 'peak_process_tree_RSS_bytes':0,
             'global_stop':None,'method_stops':{},'jobs':[],'status':'RUNNING_'+prefix}
    started = state['first_source_G_reference_or_prediction_monotonic']
    shared_cache, pearson_cache, singler_cache = root/'cache/shared', root/'cache/pearson', root/'cache/singler'
    stopped = set()

    def save():
        state['elapsed_seconds'] = time.monotonic()-started
        state['byte_accounting_before_state_write'] = ledger.sample()
        ledger.write_json(root/'RUN_STATE.json',state)
        ledger.write_json(root/'OBSERVATIONS.json',payload)

    def global_reason():
        usage = ledger.sample()
        rss = process_tree_rss(os.getpid())
        state['peak_process_tree_RSS_bytes'] = max(state['peak_process_tree_RSS_bytes'],rss)
        if time.monotonic()-started >= caps['whole_seconds']: return 'WHOLE_WALL_CAP'
        if rss >= caps['rss_bytes']: return 'PROCESS_TREE_RSS_CAP'
        if usage['cache'] >= caps['cache_bytes']: return 'CUMULATIVE_CACHE_CAP'
        if usage['output'] >= caps['output_bytes']: return 'CUMULATIVE_OUTPUT_CAP'
        return None

    def job(action, local_deadline, **extra):
        reason = global_reason()
        if reason:
            state['global_stop'] = reason; return {'kind':'GLOBAL_BUDGET_STOP','reason':reason}
        number = len(state['jobs'])+1
        folder = root/'jobs'/f'{number:03d}_{action}'
        folder.mkdir()
        result_path = folder/'RESULT.json'
        spec = {'mode':args.mode,'action':action,'manifest':str(Path(args.manifest).absolute()),
                'shared_cache':str(shared_cache),'pearson_cache':str(pearson_cache),'singler_cache':str(singler_cache),
                'result':str(result_path),**extra}
        ledger.write_json(folder/'JOB.json',spec)
        record = {'index':number,'action':action,'started_epoch':time.time(),'started_monotonic':time.monotonic()}
        state['jobs'].append(record)
        with (folder/'stdout.txt').open('w',encoding='utf-8') as stdout,(folder/'stderr.txt').open('w',encoding='utf-8') as stderr:
            process = subprocess.Popen([SCIENCE_PYTHON,str(Path(__file__).with_name('comparison_worker.py')),
                '--mode',args.mode,'--job',str(folder/'JOB.json')],
                env={**os.environ,**THREAD_ENV,'PYTHONDONTWRITEBYTECODE':'1'},stdout=stdout,stderr=stderr)
            record['pid'] = process.pid
            killed = None
            while process.poll() is None:
                reason = global_reason()
                if reason:
                    state['global_stop'] = reason; killed = 'GLOBAL_BUDGET_STOP'
                elif time.monotonic() >= local_deadline:
                    reason = 'LOCAL_STAGE_WALL_CAP'; killed = 'LOCAL_BUDGET_STOP'
                if killed:
                    record['termination'] = terminate_tree(process.pid)
                    process.wait(); break
                time.sleep(0.05)
        record.update(returncode=process.returncode,finished_epoch=time.time(),elapsed_seconds=time.monotonic()-record['started_monotonic'])
        ledger.sample()
        if killed:
            result = {'kind':killed,'reason':reason}
        elif result_path.exists():
            result = read(result_path)
        else:
            result = {'kind':'METHOD_REFERENCE_FAILURE' if 'score_' not in action else 'QUERY_FAILURE',
                      'error':'Worker exited without result','returncode':process.returncode}
        record['result_kind'] = result['kind']
        # Capture byte/RSS caps even when a fast child completed between polls.
        reason = global_reason()
        if reason:
            state['global_stop'] = reason
        save()
        return result

    def shared_failed(result):
        if result['kind'] == 'SHARED_CONTRACT_FAILURE':
            payload['shared_input_contract_status'] = 'SHARED_CONTRACT_FAILURE'
            payload['source_input_contract_status'] = 'SHARED_CONTRACT_FAILURE'
            state['global_stop'] = 'SHARED_INPUT_OR_AXIS_CONTRACT_FAILURE'
            return True
        return False

    save()
    first_build = time.monotonic()
    result = job('prepare_pearson',first_build+overrides.get('method_build_seconds',{}).get('PEARSON',caps['build_seconds']))
    shared_ready = shared_cache/'SHARED_PREPARED_READY.json'
    if shared_ready.exists():
        payload['shared_input_contract_status'] = 'PASS'
        payload['source_input_contract_status'] = 'PASS'
    if not shared_failed(result):
        if result['kind'] == 'DONE':
            payload['method_reference_states']['PEARSON'] = {'status':'DONE'}
        else:
            stopped.add('PEARSON'); state['method_stops']['PEARSON'] = result
            payload['method_reference_states']['PEARSON'] = {'status':'BUDGET_STOP_NOT_RUN' if 'BUDGET' in result['kind'] else 'METHOD_CONTRACT_FAILURE','reason':result}
        if not shared_ready.exists() and not state['global_stop']:
            state['global_stop'] = 'SHARED_PREPARATION_NOT_COMPLETED'
    if not state['global_stop']:
        result = job('build_singler',time.monotonic()+overrides.get('method_build_seconds',{}).get('SINGLER',caps['build_seconds']))
        if result['kind'] == 'DONE':
            payload['method_reference_states']['SINGLER'] = {'status':'DONE','native_mocked':False}
        else:
            shared_failed(result)
            stopped.add('SINGLER'); state['method_stops']['SINGLER'] = result
            payload['method_reference_states']['SINGLER'] = {'status':'BUDGET_STOP_NOT_RUN' if 'BUDGET' in result['kind'] else 'METHOD_CONTRACT_FAILURE','reason':result}
    save()
    for index,row in enumerate(rows):
        if state['global_stop'] or len(stopped) == len(METHODS): break
        query_folder = root/'queries'/row['sample_id']; query_folder.mkdir()
        common = None
        for owner in METHODS:
            if owner in stopped or state['global_stop']: continue
            query_started = time.monotonic()
            # A stopped method never retries; another method owns a distinct
            # common-read attempt/cache under its own remaining query budget.
            query_cache = query_folder/owner/'query_G.f64le'
            common = job('prepare_query',query_started+overrides.get('method_query_seconds',{}).get(owner,caps['query_seconds']),
                         query_id=row['sample_id'],query_cache=str(query_cache))
            if shared_failed(common) or state['global_stop']: break
            if common['kind'] == 'DONE' or common['kind'] == 'QUERY_FAILURE': break
            stopped.add(owner); state['method_stops'][owner] = common
        if state['global_stop']: break
        if common is None: break
        if common['kind'] != 'DONE':
            if common['kind'] == 'QUERY_FAILURE':
                for method in METHODS:
                    if method not in stopped:
                        payload['records'][index]['methods'][method] = {'status':'TECHNICAL_FAILURE_ZERO','failure_reason':common,
                                                                       'reference_status_at_scoring':'DONE'}
            save(); continue
        for method in METHODS:
            if method in stopped or state['global_stop']: continue
            start = query_started if method == owner else time.monotonic()
            result = job('score_pearson' if method == 'PEARSON' else 'score_singler',
                         start+overrides.get('method_query_seconds',{}).get(method,caps['query_seconds']),
                         query_id=row['sample_id'],query_cache=str(query_cache),native_output=str(query_folder/'native_singler'))
            if shared_failed(result): break
            if result['kind'] == 'DONE':
                outcome = result['prediction']
                payload['records'][index]['methods'][method] = {**outcome,'reference_status_at_scoring':'DONE'}
            elif result['kind'] == 'QUERY_FAILURE':
                payload['records'][index]['methods'][method] = {'status':'TECHNICAL_FAILURE_ZERO','failure_reason':result,
                                                               'reference_status_at_scoring':'DONE'}
            elif result['kind'] == 'METHOD_REFERENCE_FAILURE':
                stopped.add(method); state['method_stops'][method] = result
                payload['method_reference_states'][method] = {'status':'METHOD_CONTRACT_FAILURE','reason':result}
            else:
                stopped.add(method); state['method_stops'][method] = result
        save()
        print(json.dumps({'query':row['sample_id'],'methods':{m:payload['records'][index]['methods'][m]['status'] for m in METHODS},
                          'global_stop':state['global_stop']}),flush=True)
    state['status'] = prefix+'_STOPPED' if state['global_stop'] or stopped else prefix+'_FINISHING'
    state['synthetic_execution_enabled'] = False
    state['real_data_execution_enabled'] = False
    save()
    state['method_completion'] = {method:'NOT_EXPORTED_BUDGET_OR_CONTRACT_STOP' for method in METHODS}
    # Export is another monitored child, sharing the original global clock/RSS/bytes.
    if not state['global_stop']:
        result = job('export',started+caps['whole_seconds'],roster=str(roster),
                     observations=str(root/'OBSERVATIONS.json'),export_output=str(root/'results'))
        if result['kind'] == 'DONE':
            state['method_completion'] = result['method_completion']
        elif not state['global_stop']:
            state['global_stop'] = 'EXPORT_CONTRACT_FAILURE'
            state['export_failure'] = result
    reason = global_reason()
    if reason: state['global_stop'] = reason
    if not state['global_stop']:
        archive = root/'RUN_BUNDLE.zip'
        state['archive_role'] = 'Stage snapshot; current terminal RUN_STATE and FINAL_BUDGET outside this archive are authoritative'
        save()
        before_archive = [p for p in root.rglob('*') if p.is_file()]
        with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED) as package:
            for path in before_archive:
                reason = global_reason()
                if reason:
                    state['global_stop'] = reason; break
                with path.open('rb') as source,package.open(path.relative_to(root).as_posix(),'w') as target:
                    while True:
                        reason = global_reason()
                        if reason:
                            state['global_stop'] = reason; break
                        chunk = source.read(65536)
                        if not chunk: break
                        target.write(chunk)
                        package.fp.flush()
                if state['global_stop']: break
        ledger.sample()
        reason = global_reason()
        if reason: state['global_stop'] = reason
    state['status'] = prefix+'_STOPPED' if state['global_stop'] or stopped else prefix+'_COMPLETE'
    if is_real and args.delivery_root and not state['global_stop']:
        from delivery_batch17 import package_delivery
        delivery = package_delivery(root, Path(args.preparation_root), Path(args.delivery_root),
                                    Path(args.receipt), state, ledger, global_reason)
        state['delivery'] = delivery
        if delivery.get('stop_reason'):
            state['global_stop'] = delivery['stop_reason']
            state['status'] = prefix+'_STOPPED'
    save()
    final_budget = {'elapsed_seconds_including_aggregation_and_run_archive':time.monotonic()-started,
                    'production_caps_unchanged':PRODUCTION_CAPS,'effective_caps':caps,
                    'global_stop':state['global_stop'],'method_stops':state['method_stops'],
                    'peak_process_tree_RSS_bytes':state['peak_process_tree_RSS_bytes'],
                    'ledger_before_own_final_record':ledger.snapshot(),
                    'real_data_execution_enabled':False,'numerical_threads':1,
                    'monitor_interval_seconds':0.05,
                    'new_actual_real_predictions':'SEE_FORMAL_PER_METHOD_RECORDS' if is_real else 'NOT_RUN'}
    ledger.write_json(root/'FINAL_BUDGET.json',final_budget)
    if is_real:
        terminal_lock = {
            'batch':17,'status':state['status'],'real_execution_authorized':True,
            'real_data_execution_enabled':False,'supervision_checkpoint':'PASSED',
            'global_stop':state['global_stop'],'method_stops':state['method_stops'],
            'original_N':len(rows),'original_D_reported':len({r['donor_id'] for r in rows}),
            'method_completion':state['method_completion'],
            'formal_first_attempt_retained':True,'automatic_retry_or_follow_on':False,
            'run_state':'RUN_STATE.json','budget':'FINAL_BUDGET.json'}
        ledger.write_json(root/'TERMINAL_EXECUTION_LOCK_20261108.json',terminal_lock)
    # Final accounting/status writes may themselves reach a tiny test cap.
    # Preserve a minimal stop record rather than leaving COMPLETE after a cap.
    reason = global_reason()
    if reason and not state['global_stop']:
        state['global_stop'] = reason
        state['status'] = prefix+'_STOPPED'
        final_budget['global_stop'] = reason
        final_budget['terminal_preservation_only_after_cap'] = True
        ledger.write_json(root/'FINAL_BUDGET.json',final_budget)
        save()
        if is_real:
            terminal_lock.update(status=state['status'], global_stop=reason,
                                 real_data_execution_enabled=False)
            ledger.write_json(root/'TERMINAL_EXECUTION_LOCK_20261108.json',terminal_lock)
    print(json.dumps({'status':state['status'],'method_completion':state['method_completion'],
                      'seconds':final_budget['elapsed_seconds_including_aggregation_and_run_archive'],
                      'bytes':ledger.sample(),'global_stop':state['global_stop']},ensure_ascii=True),flush=True)
    return 0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode',choices=('disabled','synthetic','real'),default='disabled')
    parser.add_argument('--manifest')
    parser.add_argument('--output')
    parser.add_argument('--synthetic-overrides')
    parser.add_argument('--execution-lock')
    parser.add_argument('--preparation-root')
    parser.add_argument('--delivery-root')
    parser.add_argument('--receipt')
    args = parser.parse_args()
    return run(args)


if __name__ == '__main__':
    raise SystemExit(main())
