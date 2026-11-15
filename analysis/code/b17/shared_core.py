"""Generic final-feature input with explicit synthetic or authorized real mode."""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
import csv
import json
import math
import os

import numpy as np

CLASS_IDS = ('4009', '4132', '4180', '4219', '4275', '4391', '4696', '9001', '9512')
SCALE = 'FROZEN_TRAINING_PLM_FEATURE_SCALE'


class SharedContractError(RuntimeError):
    def __init__(self, code, detail=''):
        self.code, self.detail = code, detail
        super().__init__(code + (': ' + detail if detail else ''))


class ReferenceContractError(RuntimeError):
    def __init__(self, code, detail=''):
        self.code, self.detail = code, detail
        super().__init__(code + (': ' + detail if detail else ''))


class QueryFailure(RuntimeError):
    def __init__(self, code, detail=''):
        self.code, self.detail = code, detail
        super().__init__(code + (': ' + detail if detail else ''))


def canonical(values):
    return tuple(sorted(values, key=lambda value: value.encode('utf-8')))


def _json(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def _ids(values, kind):
    if not isinstance(values, list) or not values or any(not isinstance(x, str) or not x or any(c in x for c in '\t\r\n') for x in values):
        raise SharedContractError('INVALID_' + kind.upper() + '_IDS')
    if len(set(values)) != len(values):
        raise SharedContractError('DUPLICATE_' + kind.upper() + '_IDS')
    return tuple(values)


def _fixture_file(root, location):
    path = (root / location).resolve()
    if root != path and root not in path.parents:
        raise SharedContractError('SYNTHETIC_FILE_OUTSIDE_FIXTURE_ROOT', str(path))
    return path


def _input_file(root, location, mode):
    if mode == 'synthetic':
        return _fixture_file(root, location)
    path = Path(location)
    return Path(os.path.abspath(path if path.is_absolute() else root / path))


def _manifest_mode(path, mode):
    # Default/disabled exits before reading even metadata; real flags are checked
    # before any gene-axis or expression target is opened or statted.
    if mode not in ('synthetic', 'real'):
        raise SharedContractError('REAL_EXECUTION_DISABLED')
    manifest = _json(path)
    if mode == 'synthetic':
        if manifest.get('synthetic') is not True:
            raise SharedContractError('SYNTHETIC_MANIFEST_REQUIRED')
    elif (manifest.get('synthetic') is not False or manifest.get('mode') != 'real'
          or manifest.get('real_execution_authorized') is not True
          or manifest.get('real_data_execution_enabled') is not True):
        raise SharedContractError('REAL_MANIFEST_AUTHORIZATION_REQUIRED')
    return manifest


def _cache_mode(meta, mode):
    if mode not in ('synthetic', 'real') or meta.get('synthetic') is not (mode == 'synthetic'):
        raise SharedContractError('CACHE_DATA_MODE_MISMATCH')
    # B16 synthetic caches without an explicit mode remain readable as synthetic.
    if meta.get('mode', 'synthetic') != mode:
        raise SharedContractError('CACHE_DATA_MODE_MISMATCH')


def _source_scope(mode):
    return 'REAL_FINAL_EXTERNAL' if mode == 'real' else 'SYNTHETIC_ONLY'


def _records(path, role):
    with path.open(encoding='utf-8-sig', newline='') as stream:
        rows = list(csv.DictReader(stream))
    _ids([row.get('sample_id') for row in rows], role + '_sample')
    for row in rows:
        if not isinstance(row.get('donor_id'), str) or not row['donor_id'] or any(c in row['donor_id'] for c in '\t\r\n'):
            raise SharedContractError('MISSING_DONOR_ID', row['sample_id'])
        row['truth_id'] = row.get('truth_id') or None
        if row['truth_id'] not in CLASS_IDS and (role == 'source' or row['truth_id'] is not None):
            raise SharedContractError('TRUTH_OUTSIDE_FIXED_NINE_CLASSES', row['sample_id'])
    return rows


@dataclass
class InputData:
    manifest_path: Path
    manifest: dict
    fixture_root: Path
    source_records: list
    source_values: np.ndarray
    gene_ids: tuple
    measurable_gene_ids: tuple
    query_gene_ids: tuple
    queries: list
    mode: str = 'synthetic'


@dataclass
class QueryInputData:
    manifest_path: Path
    manifest: dict
    fixture_root: Path
    query_gene_ids: tuple
    queries: list
    mode: str = 'synthetic'


@dataclass
class QueryAxis:
    gene_ids: tuple
    query_gene_indices: tuple
    synthetic: bool = True
    mode: str = 'synthetic'
    source_scope: str = 'SYNTHETIC_ONLY'
    binding_provenance: dict = field(default_factory=dict)


@dataclass
class PreparedSource:
    gene_ids: tuple
    excluded_genes: tuple
    train_record_ids: tuple
    training_donors: tuple
    supported_labels: tuple
    donor_label_pairs: tuple
    donor_label_ids: tuple
    donor_label_counts: tuple
    donor_label_profiles: tuple
    label_means: tuple
    label_contributing_donors: tuple
    query_gene_indices: tuple
    cache_dir: Path
    synthetic: bool = True
    mode: str = 'synthetic'
    source_scope: str = 'SYNTHETIC_ONLY'
    binding_provenance: dict = field(default_factory=dict)


def load_input(manifest_path, mode='disabled'):
    """Read a fixture or explicitly authorized final input; disabled is default."""
    if mode not in ('synthetic', 'real'):
        raise SharedContractError('REAL_EXECUTION_DISABLED')
    manifest_path = Path(manifest_path).resolve()
    manifest = _manifest_mode(manifest_path, mode)
    root = manifest_path.parent
    if manifest.get('scale') != SCALE:
        raise SharedContractError('SHARED_SCALE_MISMATCH')
    if tuple(manifest.get('class_ids', CLASS_IDS)) != CLASS_IDS:
        raise SharedContractError('CANONICAL_CLASS_AXIS_MISMATCH')
    source, query = manifest['source'], manifest['queries']
    gene_path = source.get('gene_axis', manifest.get('gene_axis'))
    query_gene_path = query.get('gene_axis', manifest.get('gene_axis'))
    genes = _ids(_json(_input_file(root, gene_path, mode)), 'gene')
    query_genes = _ids(_json(_input_file(root, query_gene_path, mode)), 'query_gene')
    measurable = _ids(_json(_input_file(root, manifest['measurable_gene_axis'], mode)), 'measurable_gene') if manifest.get('measurable_gene_axis') else genes
    if not set(measurable).issubset(genes):
        raise SharedContractError('MEASURABLE_GENE_OUTSIDE_SOURCE_AXIS')
    source_records = _records(_input_file(root, source['records'], mode), 'source')
    queries = _records(_input_file(root, query['records'], mode), 'query')
    if {r['truth_id'] for r in source_records} != set(CLASS_IDS):
        raise SharedContractError('ALL_NINE_SOURCE_CLASSES_REQUIRED')
    if set(r['donor_id'] for r in source_records) & set(r['donor_id'] for r in queries):
        raise SharedContractError('LITERAL_SOURCE_QUERY_DONOR_OVERLAP')
    if set(r['sample_id'] for r in source_records) & set(r['sample_id'] for r in queries):
        raise SharedContractError('SOURCE_QUERY_SAMPLE_ID_OVERLAP')
    if source.get('dtype', '<f8') != '<f8' or source.get('layout', 'sample-major') != 'sample-major' or query.get('dtype', '<f8') != '<f8':
        raise SharedContractError('FEATURE_LAYOUT_OR_DTYPE_MISMATCH')
    shape = (len(source_records), len(genes))
    if tuple(source.get('shape', shape)) != shape:
        raise SharedContractError('SOURCE_DECLARED_SHAPE_MISMATCH')
    source_file = _input_file(root, source['features'], mode)
    if source_file.stat().st_size != 8 * shape[0] * shape[1]:
        raise SharedContractError('SOURCE_BINARY_EXTENT_MISMATCH')
    # First expression read occurs here; caller starts its shared clock before load.
    source_values = np.fromfile(source_file, dtype='<f8').reshape(shape)
    for row in queries:
        row['feature_path'] = str(_input_file(root, row['feature_file'], mode))
    return InputData(manifest_path, manifest, root, source_records, source_values, genes, measurable, query_genes, queries, mode=mode)


def load_query_metadata(manifest_path, mode='disabled'):
    """Read guarded query metadata; never reopen source values."""
    if mode not in ('synthetic', 'real'):
        raise SharedContractError('REAL_EXECUTION_DISABLED')
    path = Path(manifest_path).resolve()
    manifest = _manifest_mode(path, mode)
    if manifest.get('scale') != SCALE or tuple(manifest.get('class_ids', CLASS_IDS)) != CLASS_IDS:
        raise SharedContractError('SHARED_SCALE_OR_CLASS_AXIS_MISMATCH')
    query = manifest['queries']
    if query.get('dtype', '<f8') != '<f8':
        raise SharedContractError('FEATURE_LAYOUT_OR_DTYPE_MISMATCH')
    genes = _ids(_json(_input_file(path.parent, query.get('gene_axis', manifest.get('gene_axis')), mode)), 'query_gene')
    rows = _records(_input_file(path.parent, query['records'], mode), 'query')
    for row in rows:
        row['feature_path'] = str(_input_file(path.parent, row['feature_file'], mode))
    return QueryInputData(path, manifest, path.parent, genes, rows, mode=mode)


def _ordered_mean(values):
    try:
        value = math.fsum(values) / len(values)
    except (OverflowError, ValueError) as error:
        raise SharedContractError('INVALID_TRAINING_REDUCTION', str(error)) from error
    if not math.isfinite(value):
        raise SharedContractError('INVALID_TRAINING_REDUCTION')
    return value


def prepare_source(data, cache_dir):
    """Select G on original training specimens, then form donor-label means."""
    _cache_mode(data.manifest, data.mode)
    if data.mode == 'real' and (data.manifest.get('real_execution_authorized') is not True or data.manifest.get('real_data_execution_enabled') is not True):
        raise SharedContractError('REAL_MANIFEST_AUTHORIZATION_REQUIRED')
    genes = canonical(data.gene_ids)
    index = {gene: i for i, gene in enumerate(data.gene_ids)}
    measurable = set(data.measurable_gene_ids)
    order = sorted(range(len(data.source_records)), key=lambda i: data.source_records[i]['sample_id'].encode('utf-8'))
    eligible, excluded = [], []
    for gene in genes:
        if gene not in measurable:
            excluded.append((gene, 'OUTSIDE_FIXED_M')); continue
        column = [float(data.source_values[i, index[gene]]) for i in order]
        if not all(math.isfinite(value) for value in column):
            excluded.append((gene, 'TRAIN_NONFINITE')); continue
        if max(column) <= min(column):
            excluded.append((gene, 'TRAIN_CONSTANT')); continue
        eligible.append(gene)
    if len(eligible) < 2:
        raise SharedContractError('INSUFFICIENT_G', str(len(eligible)))
    donors = canonical(set(r['donor_id'] for r in data.source_records))
    pairs, profiles, counts = [], [], []
    for donor in donors:
        for label in CLASS_IDS:
            cell = [i for i in order if data.source_records[i]['donor_id'] == donor and data.source_records[i]['truth_id'] == label]
            if not cell:
                continue
            pairs.append((donor, label)); counts.append(len(cell))
            profiles.append(tuple(_ordered_mean([float(data.source_values[i, index[g]]) for i in cell]) for g in eligible))
    means, contributions = [], []
    for label in CLASS_IDS:
        group = [profile for pair, profile in zip(pairs, profiles) if pair[1] == label]
        if not group:
            raise SharedContractError('MISSING_SOURCE_LABEL', label)
        contributions.append(len(group))
        means.append(tuple(_ordered_mean([profile[j] for profile in group]) for j in range(len(eligible))))
    query_index = {gene: i for i, gene in enumerate(data.query_gene_ids)}
    if not set(eligible).issubset(query_index):
        raise SharedContractError('QUERY_GENE_AXIS_MISSING_FIXED_G')
    ids = tuple(json.dumps([donor, label], ensure_ascii=False, separators=(',', ':')) for donor, label in pairs)
    prepared = PreparedSource(tuple(eligible), tuple(excluded), tuple(data.source_records[i]['sample_id'] for i in order), donors,
        CLASS_IDS, tuple(pairs), ids, tuple(counts), tuple(profiles), tuple(means), tuple(contributions),
        tuple(query_index[g] for g in eligible), Path(cache_dir), synthetic=data.mode == 'synthetic', mode=data.mode,
        source_scope=_source_scope(data.mode), binding_provenance=data.manifest.get('binding_provenance', {}))
    write_prepared(prepared)
    return prepared


def write_prepared(prepared):
    path = prepared.cache_dir
    path.mkdir(parents=True, exist_ok=True)
    np.asarray(prepared.donor_label_profiles, dtype='<f8').tofile(path / 'source_profiles.f64le')
    (path / 'genes.txt').write_text('\n'.join(prepared.gene_ids) + '\n', encoding='utf-8')
    with (path / 'reference.tsv').open('w', encoding='utf-8', newline='') as stream:
        stream.write('\t'.join(['gene_id', *prepared.donor_label_ids]) + '\n')
        for j, gene in enumerate(prepared.gene_ids):
            stream.write('\t'.join([gene, *(format(profile[j], '.17g') for profile in prepared.donor_label_profiles)]) + '\n')
    with (path / 'reference_meta.tsv').open('w', encoding='utf-8', newline='') as stream:
        stream.write('ref_id\tdonor\tlabel\tn_specimens\n')
        for ref_id, pair, count in zip(prepared.donor_label_ids, prepared.donor_label_pairs, prepared.donor_label_counts):
            stream.write('\t'.join([ref_id, *pair, str(count)]) + '\n')
    payload = {key: getattr(prepared, key) for key in ('gene_ids', 'excluded_genes', 'train_record_ids', 'training_donors', 'supported_labels',
        'donor_label_pairs', 'donor_label_ids', 'donor_label_counts', 'label_means', 'label_contributing_donors', 'query_gene_indices')}
    payload.update(synthetic=prepared.synthetic, mode=prepared.mode, source_scope=prepared.source_scope,
        binding_provenance=prepared.binding_provenance,
        profile_shape=[len(prepared.donor_label_profiles), len(prepared.gene_ids)],
        profile_binary='source_profiles.f64le', dtype='<f8', layout='profile-major',
        source_only_G_selected_before_aggregation=True, native_reference_orientation='reference.tsv:genes_by_profiles')
    (path / 'SOURCE_PREPARED.json').write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding='utf-8')
    (path / 'QUERY_AXIS.json').write_text(json.dumps({'synthetic': prepared.synthetic, 'mode': prepared.mode,
        'source_scope': prepared.source_scope, 'binding_provenance': prepared.binding_provenance, 'gene_ids': prepared.gene_ids,
        'query_gene_indices': prepared.query_gene_indices}, ensure_ascii=False), encoding='utf-8')


def load_query_axis(cache_dir, mode='synthetic'):
    meta = _json(Path(cache_dir) / 'QUERY_AXIS.json')
    _cache_mode(meta, mode)
    return QueryAxis(tuple(meta['gene_ids']), tuple(meta['query_gene_indices']), synthetic=meta['synthetic'], mode=mode,
                     source_scope=meta.get('source_scope', _source_scope(mode)), binding_provenance=meta.get('binding_provenance', {}))


def load_prepared(cache_dir, mode='synthetic'):
    path = Path(cache_dir)
    meta = _json(path / 'SOURCE_PREPARED.json')
    _cache_mode(meta, mode)
    profiles = np.fromfile(path / meta['profile_binary'], dtype='<f8').reshape(meta['profile_shape'])
    values = {key: tuple(tuple(x) if isinstance(x, list) else x for x in meta[key]) for key in ('gene_ids', 'excluded_genes', 'train_record_ids',
        'training_donors', 'supported_labels', 'donor_label_pairs', 'donor_label_ids', 'donor_label_counts', 'label_means',
        'label_contributing_donors', 'query_gene_indices')}
    return PreparedSource(**values, donor_label_profiles=tuple(tuple(float(x) for x in row) for row in profiles), cache_dir=path,
        synthetic=meta['synthetic'], mode=mode, source_scope=meta.get('source_scope', _source_scope(mode)), binding_provenance=meta.get('binding_provenance', {}))


def read_query(data, record):
    _cache_mode(data.manifest, data.mode)
    if data.mode == 'real' and (data.manifest.get('real_execution_authorized') is not True or data.manifest.get('real_data_execution_enabled') is not True):
        raise SharedContractError('REAL_MANIFEST_AUTHORIZATION_REQUIRED')
    path = Path(record['feature_path'])
    if not path.is_file():
        raise SharedContractError('QUERY_FEATURE_FILE_MISSING', record['sample_id'])
    if path.stat().st_size != 8 * len(data.query_gene_ids):
        raise SharedContractError('QUERY_BINARY_EXTENT_MISMATCH', record['sample_id'])
    try:
        return np.fromfile(path, dtype='<f8')
    except OSError as error:
        raise QueryFailure('QUERY_FEATURE_READ_FAILURE', str(error)) from error


def read_query_on_G(data, record, prepared):
    values = read_query(data, record)
    return tuple(float(values[i]) for i in prepared.query_gene_indices)
