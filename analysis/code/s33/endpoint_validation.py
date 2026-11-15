"""Frozen single-label endpoints: unknown truth is an error, invalid rows are misses."""
import numpy as np


def validate_class_truth(classes, truth):
    classes = np.asarray(classes, dtype=str)
    truth = np.asarray(truth, dtype=str)
    if classes.ndim != 1 or not len(classes) or len(set(classes)) != len(classes):
        raise ValueError('Class axis must be nonempty and unique')
    if truth.ndim != 1:
        raise ValueError('Truth must have one label per sample')
    unknown = sorted(set(truth) - set(classes))
    if unknown:
        raise ValueError('Unsupported truth label(s) in frozen endpoint: ' + ', '.join(unknown))
    return classes, truth


def single_label_ranks(order, classes, truth, valid=None):
    classes, truth = validate_class_truth(classes, truth)
    order = np.asarray(order)
    if order.shape != (len(truth), len(classes)):
        raise ValueError('Ranking, class and truth axes disagree')
    valid = np.ones(len(truth), dtype=bool) if valid is None else np.asarray(valid, dtype=bool)
    if valid.shape != (len(truth),):
        raise ValueError('Validity mask must match the sample axis')
    ranks = np.argmax(classes[order] == truth[:, None], axis=1) + 1
    return np.where(valid, ranks, 0)
