"""Replay this frozen S33 cohort using the four-column ACT roster."""

import argparse
import csv
import hashlib
import io
import json
import os
import platform
import time
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime
from fractions import Fraction
from pathlib import Path
from zoneinfo import ZoneInfo

THREAD_KEYS = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS")
for setting in THREAD_KEYS:
    os.environ[setting] = "1"

import numpy as np
import requests
import scipy
from scipy.stats import rankdata


def now():
    return datetime.now(ZoneInfo("Asia/Shanghai")).isoformat(timespec="milliseconds")


def digest(body):
    return hashlib.sha256(body).hexdigest()


def locked_csv(directory, entry):
    body = (directory / entry["file"]).read_bytes()
    if digest(body) != entry["sha256"]:
        raise ValueError(f"Frozen CSV identity differs: {entry['file']}")
    return list(csv.DictReader(io.StringIO(body.decode("utf-8-sig"), newline="")))


def write_json(path, data):
    with path.open("x", encoding="utf-8", newline="\n") as stream:
        json.dump(data, stream, indent=2, ensure_ascii=False, allow_nan=False)
        stream.write("\n")


def write_csv(path, rows):
    with path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


# Parsing and endpoint definitions below derive from the frozen scorer.
class SpecimenFailure(ValueError):
    """A per-specimen failure defined in the frozen protocol."""


def parse_tpm(body, expected_header, positions):
    with io.TextIOWrapper(io.BytesIO(body), encoding="utf-8-sig", newline="") as stream:
        rows = csv.reader(stream, delimiter="\t")
        if next(rows, None) != expected_header:
            raise SpecimenFailure("unexpected_RSEM_header")
        gene_col = expected_header.index("gene_id")
        tpm_col = expected_header.index("TPM")
        values = np.empty(len(positions), dtype=np.float64)
        seen = set()
        for row in rows:
            if len(row) <= gene_col or row[gene_col] not in positions:
                continue
            gene = row[gene_col]
            if gene in seen:
                raise SpecimenFailure(f"duplicate_required_gene:{gene}")
            try:
                values[positions[gene]] = float(row[tpm_col])
            except (IndexError, ValueError) as exc:
                raise SpecimenFailure(f"nonnumeric_or_missing_TPM:{gene}") from exc
            seen.add(gene)
    if len(seen) != len(positions):
        raise SpecimenFailure(f"missing_required_genes:{len(positions) - len(seen)}")
    if not np.isfinite(values).all() or np.any(values < 0):
        raise SpecimenFailure("nonfinite_or_negative_TPM")
    if np.ptp(values) == 0:
        raise SpecimenFailure("constant_required_TPM")
    return values

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

def summarize(predictions, classes):
    grouped = defaultdict(list)
    for row in predictions:
        grouped[row["donor_group"]].append(row)
    donor_rows = []
    donor_fractions = []
    for donor, rows in grouped.items():
        hits = sum(row["top1"] for row in rows)
        rate = Fraction(hits, len(rows))
        donor_fractions.append(rate)
        donor_rows.append({
            "donor_group": donor, "specimens": len(rows), "hits": hits, "top1": float(rate),
            "technical_failures": sum(row["status"] == "TECHNICAL_FAILURE" for row in rows),
        })
    donor_equal = sum(donor_fractions) / len(donor_fractions)
    region_rows = []
    confusion = []
    for region in ("HIP", "TCx"):
        subset = [row for row in predictions if row["region"] == region]
        truth = subset[0]["truth_class_id"]
        counts = Counter(row["predicted_class_id"] or "TECHNICAL_FAILURE" for row in subset)
        destinations = {
            label: counts[label] for label in classes if label != truth and counts[label]
        }
        region_rows.append({
            "region": region, "truth_class_id": truth, "truth_name": subset[0]["truth_name"],
            "specimens": len(subset), "donors": len(subset), "hits": counts[truth],
            "misclassified": sum(destinations.values()), "technical_failures": counts["TECHNICAL_FAILURE"],
            "top1": counts[truth] / len(subset), "misclassification_destinations": destinations,
        })
        confusion.append({
            "region": region, "truth_class_id": truth,
            **{label: counts[label] for label in (*classes, "TECHNICAL_FAILURE")},
        })
    return donor_equal, donor_rows, region_rows, confusion


def fetch(url, expected_sha256):
    receipt = {"url": url, "expected_sha256": expected_sha256, "attempts": []}
    with requests.Session() as session:
        for number in range(1, 4):
            attempt = {"number": number, "started_Asia_Shanghai": now()}
            receipt["attempts"].append(attempt)
            try:
                with session.get(url, timeout=(15, 45), stream=True) as response:
                    status = response.status_code
                    attempt.update(http_status=status, response_url=response.url)
                    if status == 200:
                        body = response.content
                        attempt.update(body_bytes=len(body), sha256=digest(body))
                        if attempt["sha256"] == expected_sha256:
                            return body, dict(receipt, status="MATCH", failure_reason="")
                        reason = "unexpected_public_object_sha256"
                        retry = False
                    else:
                        reason = f"HTTP_{status}"
                        retry = status in (408, 429) or 500 <= status < 600
            except requests.RequestException as exc:
                reason = f"{type(exc).__name__}: {exc}"
                retry = True
            attempt["failure_reason"] = reason
            if not retry:
                break
            if number < 3:
                time.sleep(2 * number)
    return None, dict(receipt, status="TECHNICAL_FAILURE", failure_reason=reason)


def fetch_tpm(row, header, positions):
    body, receipt = fetch(row["tpm_url"], row["sha256"])
    if body is None:
        return None, receipt
    try:
        return parse_tpm(body, header, positions), receipt
    except (SpecimenFailure, UnicodeError, csv.Error) as exc:
        return None, dict(receipt, status="TECHNICAL_FAILURE", failure_reason=str(exc))


def check_roster(roster, config):
    expected_columns = ["tpm_url", "sha256", "brain_region", "donor_group"]
    if len(roster) != config["cohort"]["specimens"] or list(roster[0]) != expected_columns:
        raise ValueError("Unexpected four-column frozen roster shape")
    groups = {row["donor_group"] for row in roster}
    specimens = {(row["donor_group"], row["brain_region"]) for row in roster}
    urls = {row["tpm_url"] for row in roster}
    if len(groups) != config["cohort"]["donors"] or len(specimens) != len(roster) or len(urls) != len(roster):
        raise ValueError("Unexpected frozen donor grouping or duplicate specimen")
    if any(not url.startswith("https://aging.brain-map.org/") for url in urls):
        raise ValueError("TPM links must resolve from the frozen Allen source")


def score_rows(roster, tpm, valid, background, model, feature_indices, names, truth, receipts):
    classes = model["class_ids"].astype(str)
    logits = np.full((len(roster), len(classes)), np.nan)
    scores = np.full_like(logits, np.nan)
    truth_ranks = np.zeros(len(roster), dtype=int)
    sample_to_order = {}
    if valid.any():
        percentile_ranks = (rankdata(tpm[valid], axis=1, method="average") - 1) / (len(background) - 1)
        ranks = percentile_ranks[:, feature_indices]
        standardized = (ranks - model["mean"]) / model["scale"]
        logits[valid] = standardized @ model["coef"].T + model["intercept"]
        shifted = logits[valid] - logits[valid].max(axis=1, keepdims=True)
        exp_logits = np.exp(shifted)
        scores[valid] = exp_logits / exp_logits.sum(axis=1, keepdims=True)
        ordered = np.argsort(-logits[valid], axis=1, kind="stable")
        truth_ranks[valid] = single_label_ranks(ordered, classes, truth[valid])
        sample_to_order = {i: list(classes[order]) for i, order in zip(np.flatnonzero(valid), ordered)}

    predictions = []
    for i, row in enumerate(roster):
        order = sample_to_order.get(i, [])
        top_class = order[0] if order else ""
        result = {
            "specimen_key": f"S{i + 1:02d}",
            "donor_group": row["donor_group"],
            "region": row["brain_region"],
            "truth_class_id": truth[i],
            "truth_name": names[truth[i]],
            "status": "SCORED" if valid[i] else "TECHNICAL_FAILURE",
            "failure_reason": receipts[i]["failure_reason"],
            "predicted_class_id": top_class,
            "predicted_name": names.get(top_class, ""),
            "truth_rank": int(truth_ranks[i]) if valid[i] else "",
            "top1": int(truth_ranks[i] == 1),
            "ranked_class_ids": "|".join(order),
        }
        for j, label in enumerate(classes):
            result[f"logit_{label}"] = float(logits[i, j]) if valid[i] else ""
            result[f"score_{label}"] = float(scores[i, j]) if valid[i] else ""
        predictions.append(result)
    return predictions


def replay(config_path, output):
    config_bytes = config_path.read_bytes()
    config = json.loads(config_bytes)
    runtime = {
        "python": platform.python_version(), "numpy": np.__version__,
        "scipy": scipy.__version__, "requests": requests.__version__,
    }
    if runtime != config["runtime"]:
        raise ValueError(f"Unexpected runtime versions: {runtime}")
    roster = locked_csv(config_path.parent, config["roster"])
    background = locked_csv(config_path.parent, config["background"])
    check_roster(roster, config)
    model_to_position = {row["model_gene_id"]: i for i, row in enumerate(background)}
    positions = {row["ACT_gene_id"]: i for i, row in enumerate(background)}
    adapter = config["rank_adapter"]
    if len(positions) != adapter["background_n"] or len(model_to_position) != len(positions):
        raise ValueError("Frozen ACT rank background differs")

    mapping = config["truth_mapping"]
    truth = np.array([mapping[row["brain_region"]] for row in roster], dtype=str)
    names = {row["id"]: row["name"] for row in config["class_order"]}

    started = now()
    output.mkdir(parents=True, exist_ok=False)
    public_model = config["model"]
    body, model_receipt = fetch(public_model["url"], public_model["sha256"])
    if body is None:
        write_json(output / "SOURCE_RECEIPTS.json", {"model": model_receipt, "tpm": []})
        raise RuntimeError("The public fixed model could not be retrieved at its frozen hash")

    model_keys = ("selected_gene_ids", "mean", "scale", "coef", "intercept", "class_ids")
    with np.load(io.BytesIO(body), allow_pickle=False) as source:
        model = {key: source[key] for key in model_keys}
    classes = model["class_ids"].astype(str)
    if list(classes) != list(names):
        raise ValueError("Stored nine-class order differs from the frozen model order")
    validate_class_truth(classes, truth)
    feature_indices = np.array([model_to_position[gene] for gene in model["selected_gene_ids"]])
    if len(set(feature_indices)) != adapter["features_n"]:
        raise ValueError("Frozen RANK feature mapping is incomplete")

    tpm = np.full((len(roster), len(background)), np.nan)
    valid = np.zeros(len(roster), dtype=bool)
    receipts = []
    print("Fetching", len(roster), "Allen TPM objects", flush=True)
    with ThreadPoolExecutor(max_workers=6) as executor:
        jobs = [executor.submit(fetch_tpm, row, adapter["file_header"], positions) for row in roster]
        for i, job in enumerate(jobs):
            vector, receipt = job.result()
            receipts.append(dict(receipt, specimen_key=f"S{i + 1:02d}"))
            if vector is not None:
                tpm[i] = vector
                valid[i] = True
            print(receipts[-1]["specimen_key"], receipt["status"], flush=True)
    write_json(output / "SOURCE_RECEIPTS.json", {"model": model_receipt, "tpm": receipts})

    inference_started = now()
    predictions = score_rows(roster, tpm, valid, background, model, feature_indices, names, truth, receipts)
    inference_finished = now()
    donor_equal, donors, regions, confusion = summarize(predictions, list(classes))
    write_csv(output / "PREDICTIONS.csv", predictions)
    write_csv(output / "BY_DONOR.csv", donors)
    write_csv(output / "BY_REGION.csv", [
        {**row, "misclassification_destinations": json.dumps(row["misclassification_destinations"])}
        for row in regions
    ])
    write_csv(output / "CONFUSION.csv", confusion)

    hits = sum(row["top1"] for row in predictions)
    summary = {
        "start_Asia_Shanghai": started,
        "inference_started_Asia_Shanghai": inference_started,
        "inference_finished_Asia_Shanghai": inference_finished,
        "config_sha256": digest(config_bytes),
        "runner_sha256": digest(Path(__file__).read_bytes()),
        "roster_sha256": config["roster"]["sha256"],
        "background_sha256": config["background"]["sha256"],
        "model_sha256": public_model["sha256"],
        "runtime": runtime,
        "platform": platform.platform(),
        "numeric_threads": {key: os.environ[key] for key in THREAD_KEYS},
        "inference_passes": int(valid.any()),
        "classifier_fits": 0,
        "tuning_runs": 0,
        "cohort": config["cohort"],
        "scored_specimens": int(valid.sum()),
        "technical_failures": int((~valid).sum()),
        "model_misclassifications": int(valid.sum()) - hits,
        "primary_donor_equal_top1": float(donor_equal),
        "primary_exact_fraction": str(donor_equal),
        "pooled_top1_hits": hits,
        "pooled_top1": hits / len(roster),
        "by_region": regions,
    }
    write_json(output / "RESULTS.json", summary)
    print("Scored", int(valid.sum()), "of", len(roster), "specimens", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True, help="New output directory")
    args = parser.parse_args()
    replay(args.config.resolve(), args.output.resolve())
