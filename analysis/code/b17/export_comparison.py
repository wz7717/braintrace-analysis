"""Export supplied multi-method observations; no feature/model or native calls.

The Batch17 CLI requires explicit --mode synthetic or --mode real. Real supplied
observations require the frozen-final-PLM data kind and mode='real' in the API.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
import json
import math
from pathlib import Path

CLASS_IDS = ('4009', '4132', '4180', '4219', '4275', '4391', '4696', '9001', '9512')
PRIMARY_TRUTH_IDS = ('4219', '4696')
DEFAULT_METHODS = ('PEARSON', 'SINGLER')
RESOLVED = ('DONE', 'TECHNICAL_FAILURE_ZERO')
METRICS = ('top1', 'top3', 'mrr')
PLAN_SHA256 = 'af571e775d42164bee031749856c4f6377f9f4792a9d383491190e11d0dac315'
REAL_DATA_KIND = 'REAL_GSE25219_FROZEN_FINAL_PLM'


def mean(values):
    return sum(values) / len(values) if values else None


def null_metrics():
    return {name: None for name in METRICS}


def metric_means(rows):
    return {name: mean([row[name] for row in rows]) for name in METRICS}


def reference_status(value):
    return value.get('status', 'NOT_RUN') if isinstance(value, dict) else value


def ranking9(entry):
    if entry['class_ids'] != list(CLASS_IDS):
        raise ValueError('DONE requires the frozen canonical nine-class score axis')
    scores = entry['score9']
    if len(scores) != 9 or not all(math.isfinite(value) for value in scores):
        raise ValueError('DONE requires nine finite original-precision method scores')
    indices = sorted(range(9), key=lambda index: -scores[index])
    # Stable sorting keeps the ascending canonical axis for exact ties only.
    return [CLASS_IDS[index] for index in indices], scores


def normalized_rows(roster, payload, methods, truth_ids):
    observations = payload.get('records', [])
    by_key = {row['candidate_specimen_key']: row for row in observations}
    keys = [row['candidate_specimen_key'] for row in roster]
    if len(by_key) != len(observations) or len(set(keys)) != len(keys) or set(by_key) - set(keys):
        raise ValueError('Observation keys must uniquely belong to the original unique roster')
    global_shared = payload.get('shared_input_contract_status', 'NOT_RUN')
    references = payload.get('method_reference_states', {})
    rows, discrepancies = [], []
    for index, original in enumerate(roster, 1):
        key = original['candidate_specimen_key']
        truth = original.get('model_truth_id') or None
        predeclared = str(original.get('primary_included', 'true')).lower() == 'true'
        known = predeclared and truth in truth_ids
        if predeclared and truth and truth not in truth_ids:
            discrepancies.append({'candidate_specimen_key': key, 'truth_id': truth,
                                  'issue': 'UNEXPECTED_TRUTH_IN_PREDECLARED_SCOPE_STOP_REVIEW'})
        record = by_key.get(key, {})
        shared = record.get('shared_input_contract_status', global_shared)
        for method in methods:
            entry = record.get('methods', {}).get(method, {})
            status = entry.get('status', 'NOT_RUN')
            if status not in (*RESOLVED, 'NOT_RUN'):
                raise ValueError('Unknown query-method status: ' + str(status))
            ref = reference_status(references.get(method, 'NOT_RUN'))
            ref_at_score = entry.get('reference_status_at_scoring', entry.get('method_reference_status'))
            reason = entry.get('failure_reason', entry.get('stop_reason', entry.get('reason', '')))
            supplied_status = status
            # A failed/unbuilt reference is not an attempted query failure. Earlier
            # DONE outputs remain saved; historical query fail0 needs a valid-ref record.
            if status == 'TECHNICAL_FAILURE_ZERO' and ref != 'DONE' and ref_at_score != 'DONE':
                status, reason = 'NOT_RUN', 'REFERENCE_UNRESOLVED_NOT_QUERY_FAILURE'
            if status == 'TECHNICAL_FAILURE_ZERO' and shared != 'PASS':
                status, reason = 'NOT_RUN', 'SHARED_CONTRACT_UNRESOLVED_NOT_QUERY_FAILURE'
            row = {
                'input_data_kind': payload['input_data_kind'],
                'original_roster_index_1based': index,
                'candidate_specimen_key': key,
                'sample_id': original.get('sample_id', original.get('GSM_sample_id', record.get('sample_id', ''))),
                'reported_donor_case_id': original['reported_donor_case_id'],
                'model_truth_id': truth,
                'primary_predeclared': predeclared,
                'primary_included': known,
                'metadata_scope': 'PRIMARY_KNOWN_TRUTH' if known else ('UNKNOWN_TRUTH_NA' if truth is None else 'OUTSIDE_PRIMARY'),
                'method': method,
                'shared_input_contract_status': shared,
                'method_reference_status': ref,
                'reference_status_at_scoring': ref_at_score,
                'method_status': status,
                'supplied_method_status': supplied_status,
                'failure_or_stop_reason': reason,
                'source_binding_id': entry.get('source_binding_id', record.get('source_binding_id', payload.get('source_binding_id'))),
                'query_binding_id': entry.get('query_binding_id', record.get('query_binding_id', payload.get('query_binding_id'))),
                'source_binding_state': entry.get('source_binding_state', record.get('source_binding_state', payload.get('source_binding_state'))),
                'source_input_contract_status': entry.get('source_input_contract_status', record.get('source_input_contract_status', payload.get('source_input_contract_status'))),
                'class_ids': list(CLASS_IDS), 'score9': None, 'ranking9': None,
                'prediction': None, 'truth_rank': None, 'tied_scores': None,
                'uninformative_tie': None,
                'native_label': entry.get('native_label'),
                'native_class_ids': entry.get('native_class_ids'),
                'native_score9': entry.get('native_score9'),
                'native_outputs': {name: value for name, value in entry.items() if name.startswith('native_')},
                'method_observation': entry,
                **null_metrics(),
            }
            if status == 'DONE':
                ranking, scores = ranking9(entry)
                row.update(score9=scores, ranking9=ranking, prediction=ranking[0],
                           tied_scores=len(set(scores)) < 9, uninformative_tie=len(set(scores)) == 1)
                if known:
                    rank = ranking.index(truth) + 1
                    row.update(truth_rank=rank, top1=int(rank == 1), top3=int(rank <= 3), mrr=1 / rank)
            elif status == 'TECHNICAL_FAILURE_ZERO' and known:
                row.update(top1=0, top3=0, mrr=0.0)
            rows.append(row)
    return rows, discrepancies


def cell_summary(rows):
    counts = Counter(row['method_status'] for row in rows)
    n = len(rows)
    complete = n > 0 and counts['NOT_RUN'] == 0
    return {'N_specimens': n, 'N_DONE': counts['DONE'],
            'N_TECHNICAL_FAILURE_ZERO': counts['TECHNICAL_FAILURE_ZERO'], 'N_NOT_RUN': counts['NOT_RUN'],
            'status': 'NA_ABSENT_TRUTH' if not n else ('COMPLETE_CELL' if complete else 'NOT_ESTIMATED_INCOMPLETE'),
            **(metric_means(rows) if complete else null_metrics())}


def summarize(rows, discrepancies, payload, methods, truth_ids):
    primary = [row for row in rows if row['primary_included']]
    donors = sorted({row['reported_donor_case_id'] for row in primary})
    donor_rows, flat, method_summaries = [], [], {}
    shared_bad = payload.get('shared_input_contract_status', 'NOT_RUN') != 'PASS' or any(
        row['shared_input_contract_status'] == 'SHARED_CONTRACT_FAILURE' for row in rows)
    refs = payload.get('method_reference_states', {})
    for method in methods:
        selected = [row for row in primary if row['method'] == method]
        ref = reference_status(refs.get(method, 'NOT_RUN'))
        complete = bool(selected) and all(row['method_status'] in RESOLVED for row in selected)
        if discrepancies:
            completion = 'STOP_REVIEW_TRUTH_CONTRACT'
        elif shared_bad:
            completion = 'NOT_ESTIMATED_SHARED_CONTRACT_UNRESOLVED'
        elif ref != 'DONE':
            completion = 'NOT_ESTIMATED_REFERENCE_UNRESOLVED'
        else:
            completion = 'ESTIMATED_COMPLETE' if complete else 'NOT_ESTIMATED_INCOMPLETE'
        reportable = completion == 'ESTIMATED_COMPLETE'
        cells = {}
        for donor in donors:
            for scope in ('ALL', *truth_ids):
                group = [row for row in selected if row['reported_donor_case_id'] == donor
                         and (scope == 'ALL' or row['model_truth_id'] == scope)]
                cell = cell_summary(group)
                if group and (discrepancies or shared_bad or ref != 'DONE'):
                    cell.update(status=completion, **null_metrics())
                cells[(donor, scope)] = cell
                donor_rows.append({'method': method, 'reported_donor_case_id': donor,
                                   'truth_scope': scope, 'method_reference_status': ref, **cell})
        scopes = {}
        for scope in ('ALL', *truth_ids):
            group = [row for row in selected if scope == 'ALL' or row['model_truth_id'] == scope]
            contributing = [donor for donor in donors if cells[(donor, scope)]['N_specimens'] > 0]
            counts = Counter(row['method_status'] for row in group)
            pooled = metric_means(group) if reportable else null_metrics()
            equal = {metric: mean([cells[(donor, scope)][metric] for donor in contributing])
                     for metric in METRICS} if reportable else null_metrics()
            scope_status = 'NA_ABSENT_TRUTH' if not group else completion
            scopes[scope] = {'status': scope_status, 'N_specimens_original': len(group),
                             'D_contributing_original': len(contributing), 'N_DONE': counts['DONE'],
                             'N_TECHNICAL_FAILURE_ZERO': counts['TECHNICAL_FAILURE_ZERO'],
                             'N_NOT_RUN': counts['NOT_RUN'], 'pooled': pooled, 'contributing_donor_equal': equal}
            flat.append({'method': method, 'truth_scope': scope, 'method_reference_status': ref,
                         **{key: value for key, value in scopes[scope].items() if key not in ('pooled', 'contributing_donor_equal')},
                         **{'pooled_' + key: value for key, value in pooled.items()},
                         **{'donor_equal_' + key: value for key, value in equal.items()}})
        method_summaries[method] = {'role': 'FIXED_COMPARATOR_DESCRIPTIVE', 'completion_status': completion,
                                    'method_reference_status': ref, 'scopes': scopes}
    method_count = len(methods)
    predeclared = [row for row in rows if row['primary_predeclared']]
    summary = {
        'schema': 'braintrace.external_fixed_comparison.export.v1',
        'B15_frozen_plan_sha256': PLAN_SHA256,
        'canonical_class_ids': list(CLASS_IDS), 'primary_truth_ids': list(truth_ids),
        'N_metadata_scope': len(rows) // method_count,
        'N_predeclared_scope': len(predeclared) // method_count,
        'D_predeclared_scope': len({row['reported_donor_case_id'] for row in predeclared}),
        'N_primary_original': len(primary) // method_count, 'D_primary_original': len(donors),
        'N_primary_original_definition': 'Known declared truth accuracy denominator; original declared scope retained separately',
        'N_unknown_truth_metadata_scope': sum(row['metadata_scope'] == 'UNKNOWN_TRUTH_NA' for row in rows) // method_count,
        'primary_descriptive_endpoint_per_method': 'ALL contributing_donor_equal top1',
        'shared_input_contract_status': payload.get('shared_input_contract_status', 'NOT_RUN'),
        'method_reference_states': refs, 'methods': method_summaries,
        'truth_contract_discrepancies': discrepancies,
        'output_competitors_without_truth_observations': [label for label in CLASS_IDS if label not in {row['model_truth_id'] for row in primary}],
        'intervals_tests_resampling': 'NOT_IMPLEMENTED_NOT_REQUESTED',
        'BrainTrace_B13_baseline_recomputed': False,
        'native_or_model_or_feature_execution_by_this_exporter': False,
    }
    return summary, donor_rows, flat


def write_json(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + '\n', encoding='utf-8')


def write_csv(path, rows, fields):
    with path.open('w', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator='\n')
        writer.writeheader()
        for row in rows:
            writer.writerow({name: json.dumps(row.get(name), ensure_ascii=False, allow_nan=False)
                             if isinstance(row.get(name), (dict, list)) else row.get(name) for name in fields})


def export(roster_path, observations_path, output, methods=None, mode=None):
    payload = json.loads(Path(observations_path).read_text(encoding='utf-8-sig'))
    kind = payload.get('input_data_kind', '')
    if str(kind).startswith('SYNTHETIC'):
        if mode not in (None, 'synthetic'):
            raise PermissionError('Synthetic observations require synthetic mode')
        export_mode = 'synthetic'
    elif kind == REAL_DATA_KIND:
        if mode != 'real':
            raise PermissionError('Real supplied observations require explicit mode=real')
        export_mode = 'real'
    else:
        raise ValueError('Unsupported input_data_kind for the Batch17 exporter')
    with Path(roster_path).open(encoding='utf-8-sig', newline='') as stream:
        roster = list(csv.DictReader(stream))
    methods = tuple(methods or payload.get('methods', DEFAULT_METHODS))
    if not methods or len(set(methods)) != len(methods):
        raise ValueError('Explicit unique method IDs are required')
    truth_ids = tuple(payload.get('primary_truth_ids', PRIMARY_TRUTH_IDS))
    rows, discrepancies = normalized_rows(roster, payload, methods, truth_ids)
    summary, donor_rows, flat = summarize(rows, discrepancies, payload, methods, truth_ids)
    summary.update(input_data_kind=kind, export_mode=export_mode,
                   source_binding_id=payload.get('source_binding_id'),
                   query_binding_id=payload.get('query_binding_id'),
                   source_binding_state=payload.get('source_binding_state'),
                   source_input_contract_status=payload.get('source_input_contract_status'),
                   source_roster=str(Path(roster_path).absolute()),
                   source_observations=str(Path(observations_path).absolute()))
    for row in (*donor_rows, *flat):
        row['input_data_kind'] = kind
    output = Path(output); output.mkdir(parents=True, exist_ok=True)
    write_json(output / 'PER_QUERY_METHOD.json', rows)
    fields = list(rows[0]) if rows else ['input_data_kind', 'candidate_specimen_key', 'method', 'method_status']
    write_csv(output / 'PER_QUERY_METHOD.csv', rows, fields)
    ranking_rows, native_rows = [], []
    for row in rows:
        if row['method_status'] == 'DONE':
            scores = dict(zip(CLASS_IDS, row['score9']))
            for rank, label in enumerate(row['ranking9'], 1):
                ranking_rows.append({'input_data_kind': kind,
                                     'candidate_specimen_key': row['candidate_specimen_key'], 'method': row['method'],
                                     'rank': rank, 'class_id': label, 'score': scores[label],
                                     'is_known_primary_truth': label == row['model_truth_id'] if row['primary_included'] else None})
        if row['native_outputs']:
            native_rows.append({key: row[key] for key in ('input_data_kind', 'candidate_specimen_key', 'method', 'method_status',
                                                         'native_label', 'native_class_ids', 'native_score9', 'native_outputs')})
    write_csv(output / 'FULL_RANKINGS.csv', ranking_rows, ['input_data_kind', 'candidate_specimen_key', 'method', 'rank', 'class_id', 'score', 'is_known_primary_truth'])
    write_json(output / 'NATIVE_OUTPUTS.json', native_rows)
    write_csv(output / 'NATIVE_OUTPUTS.csv', native_rows, ['input_data_kind', 'candidate_specimen_key', 'method', 'method_status', 'native_label', 'native_class_ids', 'native_score9', 'native_outputs'])
    write_csv(output / 'PER_DONOR.csv', donor_rows, ['input_data_kind', 'method', 'reported_donor_case_id', 'truth_scope', 'method_reference_status',
                                                  'N_specimens', 'N_DONE', 'N_TECHNICAL_FAILURE_ZERO', 'N_NOT_RUN', 'status', *METRICS])
    write_json(output / 'SUMMARY.json', summary)
    write_csv(output / 'SUMMARY.csv', flat, list(flat[0]) if flat else ['method', 'truth_scope'])
    return summary, rows, donor_rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode', choices=('disabled', 'synthetic', 'real'), default='disabled')
    parser.add_argument('--roster', type=Path, required=True)
    parser.add_argument('--observations', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--methods', nargs='+')
    args = parser.parse_args()
    if args.mode == 'disabled':
        parser.error('Batch17 requires explicit --mode synthetic or --mode real')
    summary, _, _ = export(args.roster, args.observations, args.output, args.methods, mode=args.mode)
    print(json.dumps({'input_data_kind': summary['input_data_kind'],
                      'methods': {key: value['completion_status'] for key, value in summary['methods'].items()}}))


if __name__ == '__main__':
    main()
