"""Generate the cohort-independent symbol table from the complete Ensembl 97 gene annotation."""
import argparse
import csv
import gzip
from pathlib import Path
import re


def build(gtf, output):
    genes = {}
    with gzip.open(gtf, "rt", encoding="utf8") as stream:
        for line in stream:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if fields[2] == "gene":
                attributes = dict(re.findall(r'(\w+) "([^"]*)";', fields[8]))
                gid = attributes["gene_id"]
                genes[gid] = attributes.get("gene_name") or gid
    with output.open("w", newline="", encoding="utf8") as stream:
        table = csv.writer(stream, lineterminator="\n")
        table.writerow(["Gene.stable.ID","Gene.name"])
        table.writerows(sorted(genes.items()))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("gtf", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    build(args.gtf, args.output)
