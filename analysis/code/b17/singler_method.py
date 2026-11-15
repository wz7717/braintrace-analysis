"""Pinned SingleR adapter with explicit real/synthetic provenance and disabled default."""
from pathlib import Path
import csv
import json
import math
import os
import shutil
import subprocess

METHOD_ID = 'SingleR2.14.2_donor-balanced_score-ranking_adaptation'
CLASS_IDS = ('4009','4132','4180','4219','4275','4391','4696','9001','9512')
THREADS = ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS')
NATIVE_ENTRY = Path(__file__).with_name('native_singler_external.R')


class MethodReferenceFailure(RuntimeError):
    pass


def write_json(path, value):
    Path(path).write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False),encoding='utf-8')


def read_tsv(path):
    with Path(path).open(encoding='utf-8-sig',newline='') as stream:
        return list(csv.DictReader(stream,delimiter='\t'))


def invoke(action, outdir, arguments, mode='disabled'):
    if mode not in ('synthetic', 'real'):
        raise PermissionError('EXPLICIT_EXECUTION_MODE_REQUIRED before native input reads')
    rscript = Path(os.environ['BRAINTRACE_B17_RSCRIPT']).resolve(strict=True)
    library = Path(os.environ['BRAINTRACE_B17_R_LIBRARY']).resolve(strict=True)
    out = Path(outdir).absolute()
    out.mkdir(parents=True,exist_ok=False)
    for name in ('user','tmp'):
        (out/name).mkdir()
    script = out/'native_singler_external.R'
    shutil.copyfile(NATIVE_ENTRY,script)
    env = {**os.environ, **{key:'1' for key in THREADS},
        'R_LIBS':str(library),'R_LIBS_USER':str(library),'R_LIBS_SITE':str(library),
        'BRAINTRACE_B17_R_LIBRARY':str(library),'R_USER':str(out/'user'),
        'TMPDIR':str(out/'tmp'),'TMP':str(out/'tmp'),'TEMP':str(out/'tmp'),
        'LANG':'C','LC_ALL':'C','R_DEFAULT_PACKAGES':'NULL','R_ENVIRON_USER':'NUL','R_PROFILE_USER':'NUL'}
    env.pop('R_HOME', None)
    command = [str(rscript),'--vanilla','--default-packages=NULL',str(script),
        '--mode',mode,'--action',action,'--out',str(out),*arguments]
    with (out/'stdout.txt').open('w',encoding='utf-8') as stdout,(out/'stderr.txt').open('w',encoding='utf-8') as stderr:
        process = subprocess.run(command,env=env,cwd=out,stdout=stdout,stderr=stderr)
    write_json(out/'command_receipt.json',{'mode':mode,'action':action,'argv':command,'exit_code':process.returncode,
        'native_mocked':False,'synthetic':mode=='synthetic',
        'source_scope':'REAL_FINAL_EXTERNAL' if mode=='real' else 'SYNTHETIC_ONLY',
        'real_data_execution_enabled':mode=='real','thread_environment':{key:env[key] for key in THREADS}})
    contract = out/'contract_status.tsv'
    if not contract.is_file():
        raise MethodReferenceFailure(f'NATIVE_RUNTIME_FAILED:{out}')
    status = read_tsv(contract)[0]
    if status['status'] != 'COMPLETED' or process.returncode:
        raise MethodReferenceFailure(f"{status['code']}:{status.get('detail','')}:{out}")
    return out


def build_reference(source_cache_dir, native_cache_dir, mode='disabled'):
    if mode not in ('synthetic', 'real'):
        raise PermissionError('EXPLICIT_EXECUTION_MODE_REQUIRED')
    source = Path(source_cache_dir).absolute()
    provenance = json.loads((source/'SOURCE_PREPARED.json').read_text(encoding='utf-8-sig'))
    expected_synthetic = mode == 'synthetic'
    scope = 'SYNTHETIC_ONLY' if expected_synthetic else 'REAL_FINAL_EXTERNAL'
    if provenance.get('synthetic') is not expected_synthetic:
        raise PermissionError('Prepared-source synthetic flag does not match explicit mode')
    if mode == 'real' and (provenance.get('mode') != 'real' or provenance.get('source_scope') != scope):
        raise PermissionError('Real prepared-source mode and REAL_FINAL_EXTERNAL provenance required')
    out = invoke('build',native_cache_dir,['--reference-dir',str(source)],mode)
    markers = [row['gene_id'] for row in read_tsv(out/'marker_genes.tsv')]
    info = {'method':METHOD_ID,'status':'DONE','mode':mode,'synthetic':expected_synthetic,
        'source_scope':scope,'binding_provenance':provenance.get('binding_provenance'), 'native_mocked':False,
        'native_cache_dir':str(out),'marker_genes':markers,'marker_count':len(markers),
        'source_prepared_cache_dir':str(source),'real_data_execution_enabled':mode=='real'}
    write_json(out/'SINGLER_REFERENCE.json',info)
    return info


def write_query(path, gene_ids, values_on_G, query_id):
    def scalar(value):
        value = float(value)
        return 'NaN' if math.isnan(value) else 'Inf' if value==math.inf else '-Inf' if value==-math.inf else repr(value)
    with Path(path).open('w',encoding='utf-8',newline='') as stream:
        stream.write('gene_id\t'+query_id+'\n')
        for gene,value in zip(gene_ids,values_on_G):
            stream.write(str(gene)+'\t'+scalar(value)+'\n')


def score_query(native_cache_dir, values_on_G, gene_ids, query_id, outdir, mode='disabled'):
    if mode not in ('synthetic', 'real'):
        raise PermissionError('EXPLICIT_EXECUTION_MODE_REQUIRED')
    if mode == 'real':
        provenance = json.loads((Path(native_cache_dir)/'SINGLER_REFERENCE.json').read_text(encoding='utf-8-sig'))
        if (provenance.get('synthetic') is not False or provenance.get('mode') != 'real'
                or provenance.get('source_scope') != 'REAL_FINAL_EXTERNAL'):
            raise PermissionError('Real native-reference provenance required before query materialization')
    out = Path(outdir).absolute()
    query = out.with_name(out.name+'_query.tsv')
    query.parent.mkdir(parents=True,exist_ok=True)
    write_query(query,gene_ids,values_on_G,query_id)
    invoke('score',out,['--cache-dir',str(Path(native_cache_dir).absolute()),'--query',str(query)],mode)
    row = read_tsv(out/'native_query_status.tsv')[0]
    if row['status'] != 'OK':
        result = {'method':METHOD_ID,'status':'TECHNICAL_FAILURE_ZERO','failure_reason':row['reason'],
                  'native_label':None,'native_class_ids':None,'native_score9':None,'class_ids':list(CLASS_IDS),'score9':None}
    else:
        scores = read_tsv(out/'native_scores.tsv')[0]
        axis = [key for key in scores if key!='record_id']
        values = [float(scores[label]) for label in axis]
        result = {'method':METHOD_ID,'status':'DONE','native_label':row['native_label'],
                  'native_class_ids':axis,'native_score9':values,'class_ids':list(CLASS_IDS),
                  'score9':[scores_value for label in CLASS_IDS for axis_label,scores_value in zip(axis,values) if label==axis_label]}
    result.update(mode=mode, synthetic=mode=='synthetic',
                  source_scope='REAL_FINAL_EXTERNAL' if mode=='real' else 'SYNTHETIC_ONLY',
                  real_data_execution_enabled=mode=='real')
    write_json(out/'SINGLER_QUERY.json',result)
    return result
