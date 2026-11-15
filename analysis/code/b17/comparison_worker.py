"""One monitored build/query job for an explicitly selected data mode."""
import argparse
import json
import os
from pathlib import Path
import time
import traceback

for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[key] = '1'


def save(path, payload):
    path = Path(path)
    data = json.dumps(payload, ensure_ascii=False, indent=2, allow_nan=False)
    # Each job/ready receipt is written once at a unique path. The controller
    # consumes RESULT only after the child exits, so no rewrite/temp alias is needed.
    path.write_text(data, encoding='utf-8')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--mode', default='disabled')
    parser.add_argument('--job')
    args = parser.parse_args()
    if args.mode not in ('synthetic', 'real'):
        print('EXECUTION_MODE_DISABLED'); return 3
    import numpy as np
    import shared_core as core
    import pearson_method as pearson
    import singler_method as singler
    job = json.loads(Path(args.job).read_text(encoding='utf-8'))
    if job['mode'] != args.mode:
        return 3
    mode = args.mode
    result = {'action': job['action'], 'mode': mode,
              'input_data_kind': 'SYNTHETIC' if mode == 'synthetic' else 'REAL_GSE25219_FROZEN_FINAL_PLM',
              'started_epoch': time.time(), 'kind': 'NOT_RUN'}
    try:
        if job['action'] == 'prepare_pearson':
            data = core.load_input(job['manifest'], mode=mode)
            prepared = core.prepare_source(data, job['shared_cache'])
            save(Path(job['shared_cache']) / 'SHARED_PREPARED_READY.json',
                 {'synthetic': mode == 'synthetic', 'mode': mode, 'G': len(prepared.gene_ids), 'source_only': True})
            result['shared_prepared'] = True
            result['reference'] = pearson.build_reference(prepared, job['pearson_cache'], mode=mode)
        elif job['action'] == 'build_singler':
            result['reference'] = singler.build_reference(job['shared_cache'], job['singler_cache'], mode=mode)
        elif job['action'] == 'prepare_query':
            data = core.load_query_metadata(job['manifest'], mode=mode)
            axis = core.load_query_axis(job['shared_cache'], mode=mode)
            record = next(row for row in data.queries if row['sample_id'] == job['query_id'])
            values = core.read_query_on_G(data, record, axis)
            path = Path(job['query_cache']); path.parent.mkdir(parents=True, exist_ok=True)
            np.asarray(values, dtype='<f8').tofile(path)
            result.update(query_id=record['sample_id'], gene_ids=list(axis.gene_ids), feature_count=len(values))
        elif job['action'] in ('score_pearson', 'score_singler'):
            axis = core.load_query_axis(job['shared_cache'], mode=mode)
            values = tuple(float(x) for x in np.fromfile(job['query_cache'], dtype='<f8'))
            if len(values) != len(axis.gene_ids):
                raise core.SharedContractError('QUERY_G_CACHE_AXIS_MISMATCH')
            if job['action'] == 'score_pearson':
                reference = pearson.load_reference(Path(job['pearson_cache']) / 'PEARSON_REFERENCE.json', mode=mode)
                result['prediction'] = {'status': 'DONE', **pearson.score_query(reference, values, mode=mode)}
            else:
                result['prediction'] = singler.score_query(job['singler_cache'], values, axis.gene_ids,
                    job['query_id'], job['native_output'], mode=mode)
        elif job['action'] == 'export':
            from export_comparison import export
            summary,_,_ = export(job['roster'],job['observations'],job['export_output'],mode=mode)
            result['method_completion'] = {method:entry['completion_status'] for method,entry in summary['methods'].items()}
        else:
            raise core.SharedContractError('UNKNOWN_ACTION')
        result['kind'] = 'DONE'
    except core.SharedContractError as error:
        result.update(kind='SHARED_CONTRACT_FAILURE', error=str(error), error_code=error.code)
    except (core.ReferenceContractError, singler.MethodReferenceFailure) as error:
        result.update(kind='METHOD_REFERENCE_FAILURE', error=str(error), error_code=getattr(error, 'code', type(error).__name__))
    except core.QueryFailure as error:
        result.update(kind='QUERY_FAILURE', error=str(error), error_code=error.code)
    except Exception as error:
        kind = 'QUERY_FAILURE' if job['action'].startswith('score_') else 'METHOD_REFERENCE_FAILURE'
        result.update(kind=kind, error=str(error), error_code=type(error).__name__, traceback=traceback.format_exc())
    result['finished_epoch'] = time.time()
    save(job['result'], result)
    print(json.dumps({'action': job['action'], 'kind': result['kind']}), flush=True)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
