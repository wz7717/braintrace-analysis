"""Arguments and CSV reading for the standalone figure entries."""

import argparse
import csv
from pathlib import Path


def renderer_args(number, description):
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        "--inputs",
        type=Path,
        default=Path(__file__).resolve().parent / "inputs" / f"figure_{number}",
        help="Directory containing the saved plot tables",
    )
    parser.add_argument("--output", type=Path, required=True, help="Destination PNG path")
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    return args


def read_rows(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))
