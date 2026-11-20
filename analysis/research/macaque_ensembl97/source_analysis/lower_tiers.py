"""Build source-fold fine-region references and score the true Network gate."""
from itertools import combinations
import json
from pathlib import Path

import numpy as np

from source_data import Fisher


def write_json(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2), encoding="utf8")


class FineReferences:
    def __init__(self, data, train, eligible, output):
        self.data = data
        self.eligible = eligible
        self.fisher = Fisher(data.x, data.regions, train)
        self.parents = {r: data.networks[ix[0]] for r, ix in self.fisher.members.items()}
        self.output = Path(output)
        self.output.mkdir()
        self.lookup = {}
        self.index = []

    def reference(self, networks):
        from old_lower import build_resolution_groups, build_region_reference

        key = tuple(sorted(networks))
        if key in self.lookup:
            return self.lookup[key]
        candidates = sorted(r for r, network in self.parents.items() if network in key)
        members = {r: self.fisher.members[r] for r in candidates}
        panel = self.fisher.panel(self.eligible, candidates, 200)
        genes = self.data.symbols[panel]
        expression = self.data.x[panel]
        rows = np.arange(len(panel))
        ref = build_region_reference(expression, candidates, members)
        annotations, _ = build_resolution_groups(
            expression, candidates, members, self.parents, rows,
            min_resolution_samples=8, min_merge_samples=3, min_pair_errors=2,
            min_confusion_rate=.15, similarity_threshold=.95,
            merge_similarity_threshold=.90, max_group_size=8,
        )
        folder = self.output / f"beam_{len(self.lookup)+1:03d}"
        folder.mkdir()
        np.savez_compressed(folder / "INPUT.npz", genes=genes, gene_rows=np.searchsorted(self.eligible, panel),
                            regions=np.array(candidates), reference=ref)
        write_json(folder / "RULES_AND_MEMBERS.json", {"beam":key, "candidates":candidates,
                   "annotations":annotations, "training_sample_indices":{r:ix.tolist() for r,ix in members.items()},
                   "n_groups":len({a["resolution_group"] for a in annotations.values()})})
        # Evaluation looks up selected rows on the complete, count-derived expression axis.
        model = {"candidates":candidates, "order":panel, "ref":ref, "ann":annotations, "folder":str(folder.resolve())}
        self.lookup[key] = model
        self.index.append({"id":folder.name, "networks":key})
        return model

    def evaluate(self, idx, gate):
        from old_lower import lower_predict

        ref = self.reference(gate)
        return lower_predict(ref, self.data.x[:, idx], self.data.symbols, self.data.regions[idx])

    def complete_final(self, classes):
        for key in combinations(sorted(classes), 3):
            self.reference(key)

    def save_index(self):
        write_json(self.output / "BEAMS.json", self.index)

