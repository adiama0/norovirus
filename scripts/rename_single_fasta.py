import argparse
from Bio import SeqIO


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_fasta")
    parser.add_argument("output_fasta")
    parser.add_argument("new_id")
    args = parser.parse_args()

    record = SeqIO.read(args.input_fasta, "fasta")
    record.id = args.new_id
    record.name = args.new_id
    record.description = ""
    SeqIO.write(record, args.output_fasta, "fasta")


if __name__ == "__main__":
    main()