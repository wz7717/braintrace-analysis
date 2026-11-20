"""BrainTrace numerical helpers used by the Ensembl 97 reproduction route."""
from __future__ import annotations
import hashlib
from pathlib import Path

def sha256(path: Path) -> str:
    result = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


def require_hash(path: Path, expected: str) -> None:
    if sha256(path) != expected:
        raise ValueError(f"SHA256 mismatch: {path}")
