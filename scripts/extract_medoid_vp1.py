#!/usr/bin/env python3

import argparse

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


# GenBank records that contain the nucleotide sequence but do not
# contain a usable VP1 annotation can be handled here.
#
# Coordinates are GenBank-style:
#   1-based
#   inclusive start and end
MANUAL_VP1_COORDINATES = {
    "LC122908": (5067, 6710),
}


SEARCH_TERMS = (
    "vp1",
    "orf2",
    "major capsid",
    "capsid protein",
)


def feature_text(feature):
    values = []

    for key in (
        "gene",
        "product",
        "note",
        "label",
        "standard_name",
    ):
        values.extend(
            str(value)
            for value in feature.qualifiers.get(key, [])
        )

    return " ".join(values).lower()


def find_vp1_feature(record):
    candidates = []

    for feature in record.features:

        if feature.type != "CDS":
            continue

        text = feature_text(feature)

        if any(term in text for term in SEARCH_TERMS):
            candidates.append(feature)

    if not candidates:
        return None

    def score(feature):
        text = feature_text(feature)

        value = 0

        if "vp1" in text:
            value += 100

        if "orf2" in text:
            value += 50

        if "major capsid" in text:
            value += 25

        if "capsid protein" in text:
            value += 10

        return value

    candidates.sort(
        key=score,
        reverse=True
    )

    return candidates[0]


def main():
    parser = argparse.ArgumentParser(
        description="Extract VP1 from a medoid GenBank record."
    )

    parser.add_argument("--genbank", required=True)
    parser.add_argument("--accession", required=True)
    parser.add_argument("--group", required=True)
    parser.add_argument("--vp1-type", required=True)
    parser.add_argument("--output", required=True)

    args = parser.parse_args()

    record = SeqIO.read(
        args.genbank,
        "genbank"
    )

    vp1_feature = find_vp1_feature(record)

    if vp1_feature is not None:

        vp1_sequence = vp1_feature.extract(
            record.seq
        )

        print(
            f"VP1 annotation found: "
            f"{vp1_feature.location}"
        )

    elif args.accession in MANUAL_VP1_COORDINATES:

        start, end = MANUAL_VP1_COORDINATES[
            args.accession
        ]

        # Python is 0-based and end-exclusive.
        # GenBank 5067..6710 becomes [5066:6710].
        vp1_sequence = record.seq[
            start - 1:end
        ]

        print(
            f"No usable VP1 annotation found for "
            f"{args.accession}."
        )

        print(
            f"Using manual VP1 coordinates "
            f"{start}..{end}."
        )

    else:

        raise ValueError(
            f"Could not identify VP1 in {args.genbank}. "
            f"If this GenBank record is unannotated, add verified "
            f"coordinates for {args.accession} to "
            f"MANUAL_VP1_COORDINATES."
        )

    if not 1500 <= len(vp1_sequence) <= 1800:
        raise ValueError(
            f"VP1 length for {args.accession} is "
            f"{len(vp1_sequence)} nt, outside the expected "
            f"1500-1800 nt range."
        )

    output_record = SeqRecord(
        vp1_sequence,
        id=(
            f"GROUP_{args.group}|"
            f"{args.vp1_type}|"
            f"{args.accession}"
        ),
        description=""
    )

    SeqIO.write(
        [output_record],
        args.output,
        "fasta"
    )

    print(f"VP1 length: {len(vp1_sequence)} nt")
    print(f"Wrote: {args.output}")


if __name__ == "__main__":
    main()
