import json, math, sys
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code" / "b17"))

from shared_core import CLASS_IDS
from pearson_method import build_reference, score_query

labels = CLASS_IDS
means = tuple(tuple(float(i == j) for j in range(len(labels))) for i in range(len(labels)))
source = SimpleNamespace(synthetic=True, mode="synthetic", label_means=means,
                         supported_labels=labels, gene_ids=tuple(f"synthetic_{n}" for n in range(len(labels))),
                         source_scope="SYNTHETIC_ONLY", binding_provenance={})
with TemporaryDirectory(prefix="b17_synthetic_") as directory:
    reference = build_reference(source, Path(directory), mode="synthetic")
    answer = score_query(reference, means[0], mode="synthetic")
scores = answer["score9"]
assert len(scores) == 9 and answer["synthetic"]
assert math.isclose(scores[0], 1.0, abs_tol=1e-12)
assert all(math.isclose(x, -0.125, abs_tol=1e-12) for x in scores[1:])
print(json.dumps({"scores": scores, "scope": "synthetic-only"}))
