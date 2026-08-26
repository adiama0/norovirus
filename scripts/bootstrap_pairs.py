import argparse
import gzip
import numpy as np
from Bio import SeqIO


parser = argparse.ArgumentParser()
parser.add_argument("--r1", required=True)
parser.add_argument("--r2", required=True)
parser.add_argument("--num-pairs", type=int, required=True)
parser.add_argument("--seed", type=int, required=True)
parser.add_argument("--out-r1", required=True)
parser.add_argument("--out-r2", required=True)

args = parser.parse_args()

def get_format(path):
    path = path.lower()

    if path.endswith(".gz"):
        path = path[:-3]

    if path.endswith(".fastq") or path.endswith(".fq"):
        return "fastq"

    if path.endswith(".fasta") or path.endswith(".fa"):
        return "fasta"

    raise ValueError(
        f"Could not determine sequence format for {path}"
    )


def open_input(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")

    return open(path, "r")

r1_format = get_format(args.r1)
r2_format = get_format(args.r2)

if r1_format != r2_format:
    raise ValueError(
        "R1 and R2 must have the same sequence format."
    )

# Generate bootstrap indices
rng = np.random.default_rng(args.seed)

indices = rng.integers(
    low=0,
    high=args.num_pairs,
    size=args.num_pairs,
    dtype=np.int64
)

indices.sort()

# Read R1 and R2 together and write selected pairs
with open_input(args.r1) as r1_handle, \
     open_input(args.r2) as r2_handle, \
     open(args.out_r1, "w") as out_r1_handle, \
     open(args.out_r2, "w") as out_r2_handle:

    r1_records = SeqIO.parse(
        r1_handle,
        r1_format
    )

    r2_records = SeqIO.parse(
        r2_handle,
        r2_format
    )

    bootstrap_position = 0

    for pair_index, (r1_record, r2_record) in enumerate(
        zip(r1_records, r2_records)
    ):
        while (
            bootstrap_position < args.num_pairs
            and indices[bootstrap_position] == pair_index
        ):
            out_r1_handle.write(
                r1_record.format(r1_format)
            )

            out_r2_handle.write(
                r2_record.format(r2_format)
            )

            bootstrap_position += 1

        if bootstrap_position == args.num_pairs:
            break


if bootstrap_position != args.num_pairs:
    raise ValueError(
        f"Expected {args.num_pairs} read pairs, "
        f"but only sampled {bootstrap_position}."
    )