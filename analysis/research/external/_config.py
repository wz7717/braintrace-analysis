import argparse
import json
from pathlib import Path


def arguments(description):
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    config = args.config.resolve()
    values = json.loads(config.read_text(encoding="utf-8-sig"))
    args.output.mkdir(parents=True, exist_ok=False)
    return values, config.parent, args.output.resolve()


def location(config, base, name):
    return (base / config[name]).resolve()
