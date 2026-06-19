from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
import argparse
import shutil
import sys

SEARCH_TERMS = ("vp1", "orf2", "capsid", "major capsid")

def qualifier_strings(feature):
    values = []
    for key in ("gene", "product", "note", "label", "standard_name"):
        for value in feature.qualifiers.get(key, []):
            values.append(str(value).lower())
    values.append(feature.type.lower())
    return values

def feature_matches_target(feature, target):
    target = target.lower()
    values = qualifier_strings(feature)
    haystack = " ".join(values)

    if target == "vp1":
        return any(term in haystack for term in SEARCH_TERMS)

    return "gene" in feature.qualifiers and feature.qualifiers["gene"][0].lower() == target

def feature_score(feature):
    haystack = " ".join(qualifier_strings(feature))
    score = 0
    if feature.type == "CDS":
        score += 5
    if "vp1" in haystack:
        score += 20
    if "orf2" in haystack:
        score += 10
    if "major capsid" in haystack:
        score += 8
    elif "capsid" in haystack:
        score += 5
    return score

def choose_feature(record, target):
    candidates = [f for f in record.features if feature_matches_target(f, target)]
    if not candidates:
        return None
    candidates.sort(key=feature_score, reverse=True)
    return candidates[0]

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference", required=True, help="Input GenBank file")
    parser.add_argument("--gene", required=True, help="Gene to extract, e.g. VP1")
    parser.add_argument("--output", required=True, help="Output GenBank file")
    args = parser.parse_args()

    if args.gene.lower() == "genome":
        shutil.copyfile(args.reference, args.output)
        return 0

    found = False

    for seq_record in SeqIO.parse(args.reference, "genbank"):
        feature = choose_feature(seq_record, args.gene)
        if feature is None:
            continue

        gene_seq = feature.location.extract(seq_record.seq)

        gene_record = SeqRecord(
            gene_seq,
            id=seq_record.id.split(".")[0],
            name=args.gene,
            description=""
        )

        source_feature = SeqFeature(
            FeatureLocation(0, len(gene_seq)),
            type="source",
            qualifiers={"mol_type": "genomic RNA", "gene": args.gene}
        )
        gene_record.features.append(source_feature)

        cds_feature = SeqFeature(
            FeatureLocation(0, len(gene_seq)),
            type="CDS",
            qualifiers={"gene": args.gene}
        )
        gene_record.features.append(cds_feature)

        gene_record.annotations["molecule_type"] = "RNA"

        SeqIO.write(gene_record, args.output, "genbank")
        found = True
        break

    if not found:
        print(f"ERROR: could not find {args.gene} in {args.reference}", file=sys.stderr)
        return 1

    return 0

if __name__ == "__main__":
    raise SystemExit(main())
