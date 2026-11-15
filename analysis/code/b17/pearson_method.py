"""Frozen B8 Pearson arithmetic, separated from the old PhaseA wrapper."""
import json
import math
from pathlib import Path

from shared_core import CLASS_IDS, QueryFailure, ReferenceContractError

METHOD_ID = 'all-eligible-gene_Pearson_donor-label_means'


def _centered(values):
    if not all(math.isfinite(v) for v in values):
        return None, None, 'NONFINITE_VALUES'
    try:
        mean = math.fsum(values) / len(values)
        centered = tuple(v - mean for v in values)
        norm = math.fsum(v * v for v in centered)
    except (OverflowError, ValueError):
        return None, None, 'REDUCTION_OVERFLOW_OR_UNDEFINED'
    return centered, norm, None


def _mode_matches(value, mode):
    if mode not in ('synthetic', 'real') or value.get('synthetic') is not (mode == 'synthetic') or value.get('mode', 'synthetic') != mode:
        raise ReferenceContractError('PEARSON_DATA_MODE_MISMATCH')


def build_reference(prepared, cache_dir, mode='synthetic'):
    _mode_matches({'synthetic': prepared.synthetic, 'mode': prepared.mode}, mode)
    references = []
    for values in prepared.label_means:
        centered, norm, reason = _centered(values)
        if centered is None or norm is None or not math.isfinite(norm) or norm <= 0:
            raise ReferenceContractError('INVALID_TRAINING_REFERENCE_NORM', reason or '')
        references.append({'centered': centered, 'norm': norm})
    result = {'method_id': METHOD_ID, 'class_ids': list(prepared.supported_labels), 'gene_ids': list(prepared.gene_ids),
              'references': references, 'synthetic': prepared.synthetic, 'mode': mode, 'source_scope': prepared.source_scope,
              'binding_provenance': prepared.binding_provenance, 'reference_type': 'equal_contributing_donor_label_centroids'}
    path = Path(cache_dir)
    path.mkdir(parents=True, exist_ok=True)
    (path / 'PEARSON_REFERENCE.json').write_text(json.dumps(result, indent=2, allow_nan=False), encoding='utf-8')
    return result


def load_reference(path, mode='synthetic'):
    reference = json.loads(Path(path).read_text(encoding='utf-8'))
    _mode_matches(reference, mode)
    return reference


def score_query(reference, values_on_G, mode='synthetic'):
    _mode_matches(reference, mode)
    if tuple(reference['class_ids']) != CLASS_IDS or len(values_on_G) != len(reference['gene_ids']):
        raise ReferenceContractError('PEARSON_AXIS_MISMATCH')
    centered, norm, reason = _centered(values_on_G)
    if centered is None:
        raise QueryFailure('QUERY_' + reason)
    if norm is None or not math.isfinite(norm) or norm <= 0:
        raise QueryFailure('QUERY_ZERO_OR_INVALID_NORM')
    scores = []
    for ref in reference['references']:
        denominator_squared = norm * ref['norm']
        if not math.isfinite(denominator_squared) or denominator_squared <= 0:
            raise QueryFailure('QUERY_UNDEFINED_SCORE')
        try:
            numerator = math.fsum(a * b for a, b in zip(centered, ref['centered']))
        except (OverflowError, ValueError) as error:
            raise QueryFailure('QUERY_UNDEFINED_SCORE', str(error)) from error
        value = numerator / math.sqrt(denominator_squared)
        if not math.isfinite(value):
            raise QueryFailure('QUERY_UNDEFINED_SCORE')
        scores.append(value)
    return {'method_id': METHOD_ID, 'class_ids': list(CLASS_IDS), 'score9': scores,
            'score_meaning': 'Pearson correlations; not probabilities', 'source_only_reference': True,
            'synthetic': reference['synthetic'], 'mode': mode, 'source_scope': reference.get('source_scope')}
