"""Fit, select and export the source-only Network variants."""
from itertools import combinations
from concurrent.futures import ThreadPoolExecutor
import json
import warnings

import joblib
import numpy as np
import pandas as pd
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing._data import _is_constant_feature

from source_data import Fisher

CS = (0.001, 0.01, 0.1)
PARAMS = dict(solver="lbfgs", penalty="l2", multi_class="multinomial", class_weight="balanced",
              fit_intercept=True, max_iter=5000, tol=1e-6, random_state=20261015)


class NetworkSplit:
    def __init__(self, data, held_out, variant):
        self.data = data
        self.train, self.query, self.eligible = data.split(held_out)
        fisher = Fisher(data.x, data.networks, self.train)
        self.panel = fisher.panel(self.eligible, fisher.labels, 5000)
        self.genes = data.symbols[self.panel]
        self.train_z = data.x[np.ix_(self.panel, self.train)].T.copy()
        self.query_z = data.x[np.ix_(self.panel, self.query)].T.copy()
        self.scaler = StandardScaler().fit(self.train_z)
        sd = np.sqrt(self.scaler.var_)
        special = _is_constant_feature(self.scaler.var_, self.scaler.mean_, len(self.train))
        denominator = self.scaler.scale_.copy()
        tau = .25 * np.median(sd[sd > 0])
        floor = variant in ("floor_clip","floor_only")
        clip = variant in ("floor_clip","clip_only")
        if floor:
            eligible = (sd > 0) & ~special
            denominator[eligible] = np.maximum(denominator[eligible], tau)
        self.state = dict(train_mean=self.scaler.mean_, raw_var=self.scaler.var_, raw_sd=sd,
                          denominator=denominator, special_mask=special, constant_mask=sd==0,
                          positive_near_constant_mask=special & (sd>0), original_scale=self.scaler.scale_,
                          tau=np.array(tau), floor_enabled=np.array(floor), clip_enabled=np.array(clip))
        for values in (self.train_z,self.query_z):
            values -= self.scaler.mean_
            values /= denominator
            if clip:
                np.clip(values, -5, 5, out=values)

    def fit(self, c):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=FutureWarning)
            warnings.filterwarnings("error", category=ConvergenceWarning)
            model = LogisticRegression(C=c, **PARAMS).fit(self.train_z, self.data.networks[self.train])
        return model

    def score(self, model):
        probabilities = model.predict_proba(self.query_z)
        ranks = model.classes_[np.argsort(-probabilities, axis=1, kind="stable")]
        return probabilities, ranks

    def export(self, folder, model):
        folder.mkdir()
        np.savez_compressed(folder / "state.npz", **self.state)
        joblib.dump(model, folder / "model.joblib")
        (folder / "genes.json").write_text(json.dumps(self.genes.tolist(),indent=2),encoding="utf8")
        (folder / "classes.json").write_text(json.dumps(model.classes_.tolist(),indent=2),encoding="utf8")


def score_grid(data, held, variant):
    split = NetworkSplit(data, held, variant)
    records = []
    for c in CS:
        try:
            _, ranks = split.score(split.fit(c))
        except ConvergenceWarning:
            ranks = None
        for donor in held:
            other = next((d for d in held if d != donor),"")
            member = data.donors[split.query] == donor
            if ranks is None:
                top1 = top3 = np.nan
            else:
                matches = ranks[member] == data.networks[split.query[member],None]
                top1 = matches[:,0].mean()
                top3 = matches[:,:3].any(axis=1).mean()
            records.append(dict(outer_holdout=other, validation_donor=donor, C=c,
                                n=int(member.sum()), top1=top1, top3=top3))
    return records


def all_scores(data, variant, workers):
    donors = sorted(set(data.donors))
    held = list(combinations(donors,2)) + [(donor,) for donor in donors]
    with ThreadPoolExecutor(max_workers=workers) as pool:
        parts = pool.map(lambda pair: score_grid(data,pair,variant),held)
        return pd.DataFrame(row for group in parts for row in group)


def choose_c(scores, holdout):
    # Requiring valid scores in every validation donor fold reproduces the source selection rule.
    selected = scores[scores.outer_holdout == holdout].groupby("C")
    means = selected[["top1","top3"]].mean()
    valid = selected.top1.count() == selected.size()
    means = means.loc[valid]
    if means.empty:
        raise ValueError(f"All C candidates failed for holdout {holdout!r}")
    best = means[means.top1 >= means.top1.max()-1e-12]
    return float(best[best.top3 >= best.top3.max()-1e-12].index.min())

