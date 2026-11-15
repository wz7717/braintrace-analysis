"""Summarize paired UKBEC RANK_FIXED out-of-fold predictions from S2."""

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path
from statistics import mean


METHODS = {'logistic_L2_C0.1': 'braintrace', 'RF500': 'rf500'}
HIP = '4219'


def paired_predictions(predictions_path):
    pairs = {}
    with predictions_path.open(encoding='utf-8-sig', newline='') as source:
        for row in csv.DictReader(source):
            if row['cohort'] != 'human' or row['route'] != 'RANK_FIXED':
                continue
            if row['method'] not in METHODS:
                continue
            identity = tuple(row[key] for key in ('sample_id', 'donor', 'truth', 'fold'))
            method = METHODS[row['method']]
            pair = pairs.setdefault(identity, {})
            assert method not in pair, (identity, method)
            pair[method] = row

    rows = []
    for (sample_id, donor, truth, fold), pair in pairs.items():
        result = {'sample_id': sample_id, 'donor': donor, 'truth': truth, 'fold': fold}
        for method in METHODS.values():
            prediction = pair[method]
            result[f'{method}_pred_top1'] = prediction['pred_top1']
            result[f'{method}_Top1'] = int(prediction['h1'])
            result[f'{method}_Top3'] = int(prediction['h3'])
            result[f'{method}_MRR'] = float(prediction['mrr'])
            result[f'{method}_nonHIP_pred_HIP'] = int(truth != HIP and prediction['pred_top1'] == HIP)
            result[f'{method}_failure'] = prediction['reason']
        rows.append(result)

    assert len(rows) == len({row['sample_id'] for row in rows}) == 1083
    assert len({row['donor'] for row in rows}) == 132
    return rows


def summarize(stratum, rows, endpoints):
    donor_groups = defaultdict(list)
    for row in rows:
        donor_groups[row['donor']].append(row)

    for endpoint in endpoints:
        a, b = f'braintrace_{endpoint}', f'rf500_{endpoint}'
        braintrace = mean(mean(row[a] for row in group) for group in donor_groups.values())
        rf500 = mean(mean(row[b] for row in group) for group in donor_groups.values())
        counts = Counter((row[a], row[b]) for row in rows)
        binary = endpoint != 'MRR'
        yield {
            'stratum': stratum,
            'endpoint': endpoint,
            'specimens': len(rows),
            'reported_donors': len(donor_groups),
            'braintrace_failures': sum(bool(row['braintrace_failure']) for row in rows),
            'rf500_failures': sum(bool(row['rf500_failure']) for row in rows),
            'braintrace_donor_equal': braintrace,
            'rf500_donor_equal': rf500,
            'rf500_minus_braintrace_donor_equal': rf500 - braintrace,
            'braintrace_sample_equal': mean(row[a] for row in rows),
            'rf500_sample_equal': mean(row[b] for row in rows),
            'braintrace_event_count': counts[1, 1] + counts[1, 0] if binary else None,
            'rf500_event_count': counts[1, 1] + counts[0, 1] if binary else None,
            'both_event': counts[1, 1] if binary else None,
            'braintrace_only_event': counts[1, 0] if binary else None,
            'rf500_only_event': counts[0, 1] if binary else None,
            'neither_event': counts[0, 0] if binary else None,
        }


def write_csv(output, name, rows):
    with (output / name).open('w', encoding='utf-8', newline='') as destination:
        writer = csv.DictWriter(destination, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--predictions", type=Path, required=True, help="Frozen S2 ALL_PREDICTIONS.csv")
    parser.add_argument("--output", type=Path, required=True, help="New output directory")
    args = parser.parse_args()
    rows = paired_predictions(args.predictions)
    args.output.mkdir(parents=True)
    truth_classes = sorted({row['truth'] for row in rows})
    nonhip = [row for row in rows if row['truth'] != HIP]
    strata = {
        'all_nine_labels': rows,
        'non_HIP': nonhip,
        **{label: [row for row in rows if row['truth'] == label] for label in truth_classes},
    }
    metrics = [
        summary for stratum, subset in strata.items()
        for summary in summarize(stratum, subset, ('Top1', 'Top3', 'MRR'))
    ]
    metrics.extend(summarize('non_HIP', nonhip, ('nonHIP_pred_HIP',)))
    write_csv(args.output, 'UKBEC_RANK_RF500_PAIRED_ENDPOINTS_20261114.csv', metrics)

    counts = Counter(
        (method, row['truth'], row[f'{method}_pred_top1'])
        for row in rows for method in METHODS.values()
    )
    outcomes = [
        {'method': method, 'truth': truth, 'pred_top1': prediction, 'n': count}
        for (method, truth, prediction), count in sorted(counts.items())
    ]
    write_csv(args.output, 'UKBEC_RANK_RF500_TOP1_FLOWS_20261114.csv', outcomes)
    write_csv(args.output, 'UKBEC_RANK_RF500_PAIRED_SPECIMENS_20261114.csv', rows)


if __name__ == '__main__':
    main()
