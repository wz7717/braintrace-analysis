"""Generate fresh nine-class synthetic files; never opens real bindings or values."""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import shutil

import numpy as np

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output-root', required=True)
parser.add_argument('--index-root', required=True)
args = parser.parse_args()
ROOT = Path(args.output_root).absolute()
INDEX = Path(args.index_root).absolute()
CLASS_IDS = ('4009', '4132', '4180', '4219', '4275', '4391', '4696', '9001', '9512')
GENES = [f'g{number:02d}' for number in range(36)] + ['g_CONST', 'g_NONFINITE', 'g_TINY']
QUERY_GENES = GENES[11:] + GENES[:11]
ROOT.mkdir(parents=True, exist_ok=True)
INDEX.mkdir(parents=True, exist_ok=True)
(ROOT / 'query_features').mkdir(exist_ok=True)


def write_json(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False), encoding='utf-8')


def digest(data):
    return hashlib.sha256(data).hexdigest()


def write_csv(path, rows):
    with path.open('w', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator='\n')
        writer.writeheader(); writer.writerows(rows)


def prototype(label):
    # A simple data-generating process, not a rule embedded in either scorer.
    values = np.full(36, 6.0, dtype=np.float64)
    start = CLASS_IDS.index(label) * 4
    values[start:start + 4] += np.array([5.0, 5.3, 5.1, 5.2])
    return values


coordinates = np.arange(36, dtype=np.float64)
donor_offsets = {
    0: .03 * np.sin((coordinates + 1) * .37),
    1: -.02 * np.cos((coordinates + 1) * .51),
    2: .01 * ((coordinates % 4) - 1.5),
}
source_records, source_rows = [], []
label_order = list(reversed(CLASS_IDS))
for label in label_order:
    for donor in (2, 0, 1):
        replicates = (1, 0) if donor == 0 else (0,)
        for replicate in replicates:
            source_id = f'SYN_SOURCE_D{donor}_{label}_R{replicate}'
            rep_offset = (.015 if replicate else -.015) * ((coordinates % 3) - 1) if donor == 0 else np.zeros(36)
            numeric = np.concatenate([prototype(label) + donor_offsets[donor] + rep_offset,
                                      np.array([7.0, 4.0, 2.0], dtype=np.float64)])
            source_records.append({'sample_id': source_id, 'donor_id': f'SYN_SOURCE_D{donor}', 'truth_id': label})
            source_rows.append(numeric)
source_values = np.asarray(source_rows, dtype='<f8')
source_values[0, GENES.index('g_NONFINITE')] = np.nan
tiny_row = next(index for index, row in enumerate(source_records) if row['sample_id'].startswith('SYN_SOURCE_D0_') and row['sample_id'].endswith('_R1'))
source_values[tiny_row, GENES.index('g_TINY')] = np.nextafter(np.float64(2.0), np.float64(3.0))
source_bytes = source_values.tobytes(order='C')
(ROOT / 'SOURCE.sample_major.f64le').write_bytes(source_bytes)
write_csv(ROOT / 'SOURCE_RECORDS.csv', source_records)
write_json(ROOT / 'GENE_IDS.json', GENES)
write_json(ROOT / 'QUERY_GENE_IDS.json', QUERY_GENES)
write_json(ROOT / 'MEASURABLE_GENES.json', GENES)

source_identity = []
for index, row in enumerate(source_records):
    source_identity.append({'row_index_0based': index, 'row_index_1based': index + 1,
                            'sample_id': row['sample_id'], 'reported_donor': row['donor_id'],
                            'truth_class_id': row['truth_id'], 'historical_outer_fold_trace_only': '',
                            'bank': 'SYNTHETIC_CONSTRUCTOR_NO_HUMAN_SAMPLES'})
write_csv(ROOT / 'FINAL_SOURCE_ROW_IDENTITY.csv', source_identity)

source_gene_hash = digest((ROOT / 'GENE_IDS.json').read_bytes())
query_gene_hash = digest((ROOT / 'QUERY_GENE_IDS.json').read_bytes())
reference_hash = digest(b'SYNTHETIC_FINAL_REFERENCE_IDENTITY_NO_REAL_REFERENCE')
annotation_hash = digest(b'SYNTHETIC_STATIC_ANNOTATION_IDENTITY_NO_REAL_ANNOTATION')
protocol_hash = digest(b'SYNTHETIC_SHARED_SCALE_PROTOCOL_IDENTITY_NO_REAL_PROTOCOL')

query_specifications = [
    ('SYN_Q_A1', 'SYN_QUERY_A', '4696', '4696'),
    ('SYN_Q_A2', 'SYN_QUERY_A', '4696', '4009'),
    ('SYN_Q_A3', 'SYN_QUERY_A', '4219', '4219'),
    ('SYN_Q_B1', 'SYN_QUERY_B', '4219', '4219'),
    ('SYN_Q_B2', 'SYN_QUERY_B', '4219', '4219'),
    ('SYN_Q_C1', 'SYN_QUERY_C', '4696', '4696'),
]
query_records, query_bindings = [], []
for index, (sample_id, donor_id, truth, generating_label) in enumerate(query_specifications):
    noise = .01 * np.sin((coordinates + 1) * (index + 2) * .13)
    original_order = np.concatenate([prototype(generating_label) + noise, np.array([7.0, np.nan, 2.0])])
    rotated = np.array([original_order[GENES.index(gene)] for gene in QUERY_GENES], dtype='<f8')
    feature_file = 'query_features/' + sample_id + '.f64le'
    raw = rotated.tobytes(order='C')
    (ROOT / feature_file).write_bytes(raw)
    query_records.append({'sample_id': sample_id, 'donor_id': donor_id, 'truth_id': truth,
                          'feature_file': feature_file, 'candidate_specimen_key': sample_id,
                          'reported_donor_case_id': donor_id, 'model_truth_id': truth})
    query_bindings.append({
        'original_roster_index_1based': index + 1, 'candidate_specimen_key': sample_id,
        'sample_id': sample_id, 'reported_donor_case_id': donor_id,
        'model_truth_id': truth, 'model_truth_label': 'cerebellar cortex' if truth == '4696' else 'hippocampus',
        'reported_identity_scope': 'Synthetic named groups, no real people',
        'feature_path': str(ROOT / feature_file), 'feature_file_exists': True,
        'feature_storage': 'Synthetic headerless little-endian float64 single row on declared rotated query gene axis',
        'logical_B13_python_shape': f'1 x {len(GENES)}', 'expected_gene_n': len(GENES),
        'expected_feature_bytes': len(raw), 'observed_feature_bytes_stat_only': len(raw),
        'feature_sha256_reused_from_B13_QC': digest(raw),
        'feature_sha256_recomputed_in_B15': False,
        'synthetic_feature_SHA_computed_on_fresh_generated_bytes': True,
        'B13_RAW_feature_sha256_reused': digest(raw), 'B13_RANK_feature_sha256_reused': digest(raw),
        'B13_RAW_RANK_feature_identity_same': True,
        'B13_shared_PLM_status': 'SYNTHETIC_GENERATOR_NO_PLM',
        'B13_RAW_status': 'NOT_RUN', 'B13_RANK_status': 'NOT_RUN',
        'frozen_reference_sha256_reused': reference_hash,
        'annotation_sha256_reused': annotation_hash,
        'protocol_sha256_reused': protocol_hash,
        'gene_axis_sha256_reused': query_gene_hash,
        'gene_axis_path_existing': str(ROOT / 'QUERY_GENE_IDS.json'),
        'binding_status': 'SYNTHETIC_GENERATED_SAME_SCHEMA_NOT_REAL_B13',
        'binding_issues': '', 'synthetic': True,
        'B15_Pearson_status': 'NOT_RUN', 'B15_SingleR_status': 'NOT_RUN',
        'B15_real_comparator_execution_authorized': False, 'B15_real_data_execution_enabled': False,
    })
write_csv(ROOT / 'QUERY_RECORDS.csv', query_records)
write_csv(ROOT / 'ROSTER.csv', query_records)
generic_manifest = {
    'synthetic': True, 'input_data_kind': 'SYNTHETIC_FRESH_CONSTRUCTOR',
    'scale': 'FROZEN_TRAINING_PLM_FEATURE_SCALE', 'class_ids': list(CLASS_IDS),
    'gene_axis': 'GENE_IDS.json', 'measurable_gene_axis': 'MEASURABLE_GENES.json',
    'source': {'features': 'SOURCE.sample_major.f64le', 'records': 'SOURCE_RECORDS.csv',
               'shape': list(source_values.shape), 'dtype': '<f8', 'layout': 'sample-major'},
    'queries': {'records': 'QUERY_RECORDS.csv', 'dtype': '<f8', 'gene_axis': 'QUERY_GENE_IDS.json'},
}
write_json(ROOT / 'GENERIC_MANIFEST.json', generic_manifest)
handoff = {
    'schema': 'braintrace.batch15.final_source_training_input_handoff.v1',
    'date': '2026-11-08', 'status': 'SYNTHETIC_SAME_SCHEMA_GENERATED_NOT_REAL_SOURCE',
    'synthetic': True, 'synthetic_fixture_root': str(ROOT),
    'selected_feature_file': str(ROOT / 'SOURCE.sample_major.f64le'), 'selected_root': str(ROOT),
    'gene_axis_path_existing': str(ROOT / 'GENE_IDS.json'),
    'measurable_gene_axis_path': str(ROOT / 'MEASURABLE_GENES.json'),
    'declared_layout': {'headerless': True, 'shape': list(source_values.shape), 'dtype': 'float64',
                        'endian': 'little', 'numpy_dtype': '<f8', 'layout': 'sample-major',
                        'extent_bytes': len(source_bytes), 'row_byte_stride': len(GENES) * 8},
    'identity_summary': {'specimens': len(source_records), 'unique_sample_ids': len(source_records),
                         'reported_donors': 3, 'truth_classes': list(CLASS_IDS), 'ordered_genes': len(GENES),
                         'all_source_records_are_final_training_rows': True, 'query_sample_ids_in_final_export': 0,
                         'donor_identity_limit': 'Entirely synthetic IDs; no real-source identity assertion'},
    'final_reference_chain': {'feature_sha256_inherited': digest(source_bytes),
                              'reference_sha256_inherited': reference_hash,
                              'annotation_sha256_inherited': annotation_hash,
                              'protocol_sha256_inherited': protocol_hash,
                              'gene_ids_sha256_inherited': source_gene_hash,
                              'new_hashes_computed': True,
                              'representation': 'Synthetic values generated on the declared shared scale; no real PLM/reference loaded'},
    'files': {'row_identity': 'FINAL_SOURCE_ROW_IDENTITY.csv'},
    'real_input_values_read': False, 'G_reference_marker_or_predictions_computed_by_generator': False,
}
write_json(ROOT / 'SOURCE_BINDING_HANDOFF.json', handoff)
write_json(ROOT / 'QUERY_FEATURE_BINDINGS.json', {
    'schema': 'braintrace.batch15.query_feature_metadata_binding.v1', 'date': '2026-11-08',
    'synthetic': True, 'synthetic_fixture_root': str(ROOT),
    'scope': 'Fresh synthetic analogue of B15 metadata schema; no actual B13 vectors or values',
    'canonical_nine_class_ids': list(CLASS_IDS), 'records': query_bindings,
})
provenance = {
    'input_data_kind': 'SYNTHETIC_FRESH_CONSTRUCTOR', 'synthetic': True,
    'generator': str(Path(__file__).absolute()), 'fixture_root': str(ROOT),
    'source_shape': list(source_values.shape), 'query_files': len(query_records),
    'source_donors': 3, 'query_donors': 3, 'class_axis': list(CLASS_IDS),
    'source_gene_order': GENES, 'query_gene_order_rotated_by': 11,
    'source_record_order': 'Unsorted label reverse, donor2/0/1, donor0 replicate1/0; binary rows match CSV exactly',
    'source_replicates': 'donor0 each label has2 specimens; donor1/donor2 each label have1',
    'data_formula': 'Nine prototypes with four peak genes/class, weak deterministic donor/replicate offsets, and small deterministic query sine noise',
    'G_expected_from_source': {'include': [f'g{x:02d}' for x in range(36)] + ['g_TINY'],
                             'exclude_constant': ['g_CONST'], 'exclude_nonfinite': ['g_NONFINITE'],
                             'g_TINY': '2 everywhere except one fresh source cell nextafter(2,3), exact nonconstant before donor averaging',
                             'tiny_source_row_index0': tiny_row, 'tiny_source_sample_id': source_records[tiny_row]['sample_id']},
    'query_semantics': [{'sample_id': sample_id, 'declared_truth': truth, 'generating_prototype': prototype_label}
                        for sample_id, _, truth, prototype_label in query_specifications],
    'explicit_synthetic_mislabel': 'SYN_Q_A2 truthCBC but values generated from frontal4009 prototype; scorers contain no fixture-ID rule',
    'missing_truth_cell': 'SYN_QUERY_C hasCBC only, HIP isNA',
    'query_nonfinite_outside_expected_G': 'g_NONFINITE isNaN in allqueries; eligibility remains source-only',
    'real_source_query_models_or_CELs_opened': False,
    'G_reference_marker_or_method_scores_computed': False,
    'numeric_files': [{'file': 'SOURCE.sample_major.f64le', 'bytes': len(source_bytes), 'sha256': digest(source_bytes)}]
                     + [{'file': row['feature_file'], 'bytes': (ROOT / row['feature_file']).stat().st_size,
                         'sha256': binding['feature_sha256_reused_from_B13_QC']} for row, binding in zip(query_records, query_bindings)],
}
write_json(ROOT / 'FIXTURE_PROVENANCE.json', provenance)
for name in ('GENERIC_MANIFEST.json', 'SOURCE_BINDING_HANDOFF.json', 'FINAL_SOURCE_ROW_IDENTITY.csv',
             'QUERY_FEATURE_BINDINGS.json', 'GENE_IDS.json', 'QUERY_GENE_IDS.json', 'MEASURABLE_GENES.json',
             'SOURCE_RECORDS.csv', 'QUERY_RECORDS.csv', 'ROSTER.csv', 'FIXTURE_PROVENANCE.json'):
    shutil.copyfile(ROOT / name, INDEX / name)
write_json(INDEX / 'ACTUAL_SYNTHETIC_FILE_LOCATIONS.json', {'fixture_root': str(ROOT),
           'source_binding': str(ROOT / 'SOURCE_BINDING_HANDOFF.json'),
           'source_rows': str(ROOT / 'FINAL_SOURCE_ROW_IDENTITY.csv'),
           'query_binding': str(ROOT / 'QUERY_FEATURE_BINDINGS.json'),
           'generic_manifest': str(ROOT / 'GENERIC_MANIFEST.json'), 'roster': str(ROOT / 'ROSTER.csv'),
           'numeric_copies_in_work_index': False})
print(json.dumps({'fixture_root': str(ROOT), 'source_binding': str(ROOT / 'SOURCE_BINDING_HANDOFF.json'),
                  'source_rows': str(ROOT / 'FINAL_SOURCE_ROW_IDENTITY.csv'),
                  'query_binding': str(ROOT / 'QUERY_FEATURE_BINDINGS.json'),
                  'generic_manifest': str(ROOT / 'GENERIC_MANIFEST.json'),
                  'roster': str(ROOT / 'ROSTER.csv'), 'source_shape': list(source_values.shape),
                  'query_n': len(query_records), 'G_not_built': True}, indent=2))
