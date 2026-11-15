"""Thin B15 metadata-to-reader bridge; never opens or stats numeric inputs."""
from pathlib import Path
import argparse
import csv
import json
import os

CLASS_IDS = ('4009', '4132', '4180', '4219', '4275', '4391', '4696', '9001', '9512')
SCALE = 'FROZEN_TRAINING_PLM_FEATURE_SCALE'


def _read(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def _absolute(base, value):
    # Lexical paths only: resolve()/stat() would touch the existing input target.
    path = Path(value)
    return Path(os.path.abspath(path if path.is_absolute() else Path(base) / path))


def _write_csv(path, rows, fields):
    with path.open('w', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator='\n')
        writer.writeheader(); writer.writerows(rows)


def build_manifest(source_handoff_path, source_rows_path, query_binding_path, output_dir, mode='disabled'):
    """Convert saved binding fields only; no expression targets are accessed.

    mode='synthetic' requires both input bindings to declare synthetic=true and
    a common owned fixture root containing every declared numeric/axis path.
    Explicit real mode emits a real-marked authorization metadata manifest; the
    parent controller keeps its preparation snapshot disabled until review.
    Default/disabled modes continue to emit final_external_disabled.
    """
    source_handoff_path, source_rows_path, query_binding_path = map(Path, (source_handoff_path, source_rows_path, query_binding_path))
    source, query = _read(source_handoff_path), _read(query_binding_path)
    with source_rows_path.open(encoding='utf-8-sig', newline='') as stream:
        original_source = list(csv.DictReader(stream))
    records = query['records']
    synthetic = mode == 'synthetic'
    real = mode == 'real'
    if real and (source.get('synthetic') is True or query.get('synthetic') is True):
        raise ValueError('REAL_MODE_REQUIRES_REAL_BINDING_METADATA')
    output = _absolute(Path.cwd(), output_dir)
    source_file = _absolute(source_handoff_path.parent, source['selected_feature_file'])
    source_gene = _absolute(source_handoff_path.parent, source.get('gene_axis_path_existing', str(Path(source['selected_root']) / 'inputs/GENE_IDS.json')))
    query_axes = {_absolute(query_binding_path.parent, row['gene_axis_path_existing']) for row in records}
    if len(query_axes) != 1:
        raise ValueError('QUERY_BINDING_REQUIRES_ONE_DECLARED_QUERY_AXIS')
    query_gene = next(iter(query_axes))
    query_paths = [_absolute(query_binding_path.parent, row['feature_path']) for row in records]
    measurable_value = source.get('measurable_gene_axis_path')
    measurable = _absolute(source_handoff_path.parent, measurable_value) if measurable_value else None
    if synthetic:
        if source.get('synthetic') is not True or query.get('synthetic') is not True:
            raise ValueError('BOTH_BINDINGS_MUST_BE_SYNTHETIC; real B15 bindings cannot be relabelled synthetic')
        if not source.get('synthetic_fixture_root'):
            raise ValueError('SYNTHETIC_FIXTURE_ROOT_REQUIRED')
        fixture_root = _absolute(source_handoff_path.parent, source['synthetic_fixture_root'])
        if output != fixture_root:
            raise ValueError('BRIDGED_SYNTHETIC_MANIFEST_MUST_SHARE_NUMERIC_FIXTURE_ROOT')
        for path in [source_file, source_gene, query_gene, *query_paths, *([measurable] if measurable else [])]:
            try:
                path.relative_to(fixture_root)
            except ValueError as error:
                raise ValueError('SYNTHETIC_BINDING_POINTS_OUTSIDE_OWNED_FIXTURE_ROOT') from error
    output.mkdir(parents=True, exist_ok=True)
    source_rows = [{'sample_id': row['sample_id'], 'donor_id': row.get('reported_donor', row.get('donor_id')),
                    'truth_id': row.get('truth_class_id', row.get('truth_id'))} for row in original_source]
    query_rows, roster = [], []
    for i, (row, path) in enumerate(zip(records, query_paths), 1):
        key = row.get('candidate_specimen_key', row['sample_id'])
        donor = row['reported_donor_case_id']
        truth = row.get('model_truth_id') or ''
        original_index = row.get('original_roster_index_1based', i)
        query_rows.append({'sample_id': row['sample_id'], 'donor_id': donor, 'truth_id': truth, 'feature_file': str(path),
                           'candidate_specimen_key': key, 'original_roster_index_1based': original_index})
        roster.append({'candidate_specimen_key': key, 'reported_donor_case_id': donor, 'model_truth_id': truth,
                       'sample_id': row['sample_id'], 'original_roster_index_1based': original_index,
                       'primary_included': str(row.get('primary_included', True)).lower()})
    _write_csv(output / 'SOURCE_RECORDS.csv', source_rows, ['sample_id', 'donor_id', 'truth_id'])
    _write_csv(output / 'QUERY_RECORDS.csv', query_rows,
               ['sample_id', 'donor_id', 'truth_id', 'feature_file', 'candidate_specimen_key', 'original_roster_index_1based'])
    _write_csv(output / 'ROSTER.csv', roster,
               ['candidate_specimen_key', 'reported_donor_case_id', 'model_truth_id', 'sample_id', 'original_roster_index_1based', 'primary_included'])
    declared = source['declared_layout']
    manifest = {
        'schema': 'braintrace.batch17.binding_bridge_manifest.v1',
        'synthetic': synthetic, 'mode': 'synthetic' if synthetic else ('real' if real else 'final_external_disabled'),
        'source_scope': 'SYNTHETIC_ONLY' if synthetic else 'REAL_FINAL_EXTERNAL',
        'real_execution_authorized': real,
        'real_data_execution_enabled': real, 'synthetic_execution_authorized': synthetic,
        'scale': SCALE, 'class_ids': list(query.get('canonical_nine_class_ids', CLASS_IDS)),
        'gene_axis': str(source_gene),
        'source': {'features': str(source_file), 'records': 'SOURCE_RECORDS.csv',
                   'shape': declared['shape'], 'dtype': declared.get('numpy_dtype', '<f8'), 'layout': declared.get('layout', 'sample-major'),
                   'gene_axis': str(source_gene)},
        'queries': {'records': 'QUERY_RECORDS.csv', 'dtype': '<f8', 'gene_axis': str(query_gene)},
        'roster': 'ROSTER.csv',
        'binding_provenance': {'source_handoff': str(_absolute(Path.cwd(), source_handoff_path)),
                               'source_row_identity': str(_absolute(Path.cwd(), source_rows_path)),
                               'query_binding': str(_absolute(Path.cwd(), query_binding_path)),
                               'field_conversion_only': True, 'numeric_files_opened_statted_or_rehashed': False,
                               'source_identity_declarations_reused': source.get('final_reference_chain', {}),
                               'query_identity_declarations_reused': [{'sample_id': row['sample_id'],
                                   'feature_SHA': row.get('feature_sha256_reused_from_B13_QC'),
                                   'gene_axis_SHA': row.get('gene_axis_sha256_reused')} for row in records]}
    }
    if measurable:
        manifest['measurable_gene_axis'] = str(measurable)
    target = output / 'BRIDGED_MANIFEST.json'
    target.write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding='utf-8')
    return target


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-handoff', required=True)
    parser.add_argument('--source-rows', required=True)
    parser.add_argument('--query-binding', required=True)
    parser.add_argument('--output-dir', required=True)
    parser.add_argument('--mode', choices=['disabled', 'final_external_disabled', 'synthetic', 'real'], default='disabled')
    args = parser.parse_args()
    result = build_manifest(args.source_handoff, args.source_rows, args.query_binding, args.output_dir, args.mode)
    print(json.dumps({'manifest': str(result), 'mode': args.mode, 'numeric_input_reads_or_stats': 0}, ensure_ascii=True))


if __name__ == '__main__':
    main()
