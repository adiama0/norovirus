import argparse
from collections import Counter
from Bio import SeqIO


def consensus_base(column):
    counts = Counter(base.upper() for base in column)

    # Ignore gaps and missing/ambiguous unknowns for voting
    for bad in ["-", "N", "?"]:
        counts.pop(bad, None)

    if not counts:
        return "N"

    most_common = counts.most_common()
    top_count = most_common[0][1]
    top_bases = [base for base, count in most_common if count == top_count]

    # If there is a tie, return N conservatively
    if len(top_bases) > 1:
        return "N"

    return top_bases[0]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("alignment_fasta")
    parser.add_argument("output_fasta")
    parser.add_argument("seq_name")
    args = parser.parse_args()

    records = list(SeqIO.parse(args.alignment_fasta, "fasta"))
    if not records:
        raise ValueError("No sequences found in alignment")

    length = len(records[0].seq)
    for r in records:
        if len(r.seq) != length:
            raise ValueError("Alignment sequences are not all the same length")

    consensus = []
    for i in range(length):
        column = [str(r.seq[i]) for r in records]
        consensus.append(consensus_base(column))

    seq = "".join(consensus)

    with open(args.output_fasta, "w") as out:
        out.write(f">{args.seq_name}\n")
        out.write(seq + "\n")


if __name__ == "__main__":
    main()