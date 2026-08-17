#!/usr/bin/env python3

import argparse

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

SEARCH_TERMS = (
    "vp1",
    "orf2",
    "major capsid",
    "capsid protein",
)


STOP_CODONS = {
    "TAA",
    "TAG",
    "TGA",
}


def feature_text(feature):
    """
    Combine useful GenBank feature qualifiers into one lowercase string.
    """

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
    """
    Search annotated CDS features for VP1 / ORF2 / capsid terms.
    """

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


def find_complete_orfs(sequence, min_nt=1500, max_nt=1800):
    """
    Find maximal complete forward-strand ORFs in the expected VP1 size range.
    """

    sequence = str(sequence).upper()

    candidates = []

    for frame in range(3):

        start = None

        for i in range(frame, len(sequence) - 2, 3):

            codon = sequence[i:i + 3]

            if start is None:

                if codon == "ATG":
                    start = i

            else:

                if codon in STOP_CODONS:

                    end = i + 3
                    length = end - start

                    if min_nt <= length <= max_nt:

                        candidates.append(
                            {
                                "start_python": start,
                                "end_python": end,
                                "start_genbank": start + 1,
                                "end_genbank": end,
                                "length": length,
                                "frame": frame,
                            }
                        )

                    start = None

    return candidates


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

    # Option 1: Use an actual GenBank VP1/ORF2 CDS annotation

    if vp1_feature is not None:

        vp1_sequence = vp1_feature.extract(
            record.seq
        )

        print(
            f"VP1 annotation found: "
            f"{vp1_feature.location}"
        )

    # Option 2: If there is no annotation or manual entry, search for a complete forward ORF with the expected VP1 length.
    else:

        candidates = find_complete_orfs(
            record.seq,
            min_nt=1500,
            max_nt=1800,
        )

        if len(candidates) == 1:

            candidate = candidates[0]

            vp1_sequence = record.seq[
                candidate["start_python"]:
                candidate["end_python"]
            ]

            print(
                f"No usable VP1 CDS annotation found for "
                f"{args.accession}."
            )

            print(
                "Automatically detected one complete "
                "VP1-sized ORF."
            )

            print(
                f"Using inferred coordinates: "
                f"{candidate['start_genbank']}.."
                f"{candidate['end_genbank']}"
            )

            print(
                f"Reading frame: {candidate['frame']}"
            )

        elif len(candidates) == 0:

            raise ValueError(
                f"Could not identify VP1 in {args.genbank}. "
                f"No annotated VP1 CDS, no manually verified "
                f"coordinates, and no unique complete 1500-1800 nt "
                f"forward ORF was found."
            )

        else:

            candidate_text = ", ".join(
                f"{candidate['start_genbank']}.."
                f"{candidate['end_genbank']} "
                f"({candidate['length']} nt)"
                for candidate in candidates
            )

            raise ValueError(
                f"Could not identify VP1 unambiguously in "
                f"{args.genbank}. Multiple VP1-sized complete ORFs "
                f"were found: {candidate_text}. "
                f"Inspect the record and add the verified VP1 "
                f"coordinates to MANUAL_VP1_COORDINATES."
            )


    # Checks:
    #Lengtj 
    if not 1500 <= len(vp1_sequence) <= 1800:

        raise ValueError(
            f"VP1 length for {args.accession} is "
            f"{len(vp1_sequence)} nt, outside the expected "
            f"1500-1800 nt range."
        )

    # Error in reading codons
    if len(vp1_sequence) % 3 != 0:
        raise ValueError(
            f"VP1 length for {args.accession} is "
            f"{len(vp1_sequence)} nt and is not divisible by 3."
        )

    protein = vp1_sequence.translate()
    internal_stops = protein[:-1].count("*")

    print(f"VP1 length: {len(vp1_sequence)} nt")

    print(f"Translation length: {len(protein)} aa "
        f"(including terminal stop if present)")

    # Write the FASTA 
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

    print(
        f"Wrote: {args.output}"
    )

if __name__ == "__main__":
    main()