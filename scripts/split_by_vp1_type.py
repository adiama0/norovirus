import argparse
import csv
import os
import re
from collections import defaultdict
from Bio import SeqIO


def safe_name(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", s)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("metadata_tsv")
    parser.add_argument("input_fasta")
    parser.add_argument("outdir")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    with open(args.metadata_tsv, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        rows = list(reader)

    if "name" not in rows[0] or "VP1_type" not in rows[0]:
        raise ValueError("metadata TSV must contain 'name' and 'VP1_type' columns")

    by_type = defaultdict(list)
    for row in rows:
        vp1_type = row["VP1_type"].strip()
        if not vp1_type:
            continue
        by_type[vp1_type].append(row)

    fasta_records = SeqIO.to_dict(SeqIO.parse(args.input_fasta, "fasta"))

    for vp1_type, type_rows in by_type.items():
        dirname = safe_name(vp1_type)
        type_dir = os.path.join(args.outdir, dirname)
        os.makedirs(type_dir, exist_ok=True)

        with open(os.path.join(type_dir, "type_name.txt"), "w") as f:
            f.write(vp1_type + "\n")

        with open(os.path.join(type_dir, "metadata.tsv"), "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=type_rows[0].keys(), delimiter="\t")
            writer.writeheader()
            writer.writerows(type_rows)

        selected = []
        for row in type_rows:
            seq_id = row["name"]
            if seq_id in fasta_records:
                selected.append(fasta_records[seq_id])

        SeqIO.write(selected, os.path.join(type_dir, "sequences.fasta"), "fasta")


if __name__ == "__main__":
    main()