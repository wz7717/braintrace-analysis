"""Stage the locked AHBA ID maps from verified public Git bytes."""
from argparse import ArgumentParser
from hashlib import sha256
from pathlib import Path

FILES = {
    "author_map": (
        "bo2023_bulk_atlas_buildkit/04_expressed_genes_neocortex_plus_subcortical.cleaned_symbols.csv",
        1328791, "6d1a3aa59a67c12453fcf37d82bbd64295f17507badce294807f890ccdf8fa9c",
        1352397, "24e0545a478ed2643322f994627a0d1c8bfac3061f6a755eae49080b2b92a78a",
    ),
    "orthology": (
        "data/orthology/ensembl_mfascicularis_hsapiens_homology.tsv",
        2014062, "cf4ffcd8ad88fe289d7be78c9baaff06b6d6e73a7a07496af11af7464d37d109",
        2054920, "f88271ed0136a949f096219b7094573b586ac4b0f7dca76f1cfcaf3a45288693",
    ),
}
COMMIT = "9e0044fe52c7ba455605108e0c06616d0a223536"


def checked_stage(kind, source, output):
    _, raw_n, raw_sha, locked_n, locked_sha = FILES[kind]
    raw = source.read_bytes()
    if (len(raw), sha256(raw).hexdigest()) != (raw_n, raw_sha):
        raise ValueError("The supplied file differs from the pinned LF Git artifact")
    staged = raw.replace(b"\n", b"\r\n")
    if (len(staged), sha256(staged).hexdigest()) != (locked_n, locked_sha):
        raise ValueError("Staged CRLF bytes differ from the locked analysis input")
    with output.open("xb") as stream:
        stream.write(staged)


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("kind", choices=FILES)
    parser.add_argument("source", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    checked_stage(args.kind, args.source, args.output)


if __name__ == "__main__":
    main()
