#!/bin/bash

set -euo pipefail

THREADS="${THREADS:-8}"

# Uncomment for full run
rm -rf tmp
rm -rf data/references

mkdir -p tmp
mkdir -p data/references

# Download all norovirus sequences from Nextstrain

echo "Downloading Nextstrain sequences..."

wget --show-progress \
  https://data.nextstrain.org/files/workflows/norovirus/sequences.fasta.zst \
  -O tmp/nextstrain_sequences.fasta.zst

zstd -d -f \
  tmp/nextstrain_sequences.fasta.zst \
  -o tmp/nextstrain_sequences.fasta


# Extract VP1 metadata and VP1 phylogenetic tree

echo "Extracting VP1 metadata and tree..."

barcodeforge extract-auspice-data \
  data/norovirus_all_VP1.json \
  --attributes "VP1_type" \
  --output_metadata_path tmp/metadata.tsv

mv tree.nwk \
  tmp/VP1_tree.nwk


# Keep only samples with clean GII VP1 labels

echo "Cleaning metadata..."

awk -F'\t' '
  NR == 1 {
    print
    next
  }

  $2 ~ /^GII\.[0-9]+$/ {
    print
  }
' tmp/metadata.tsv \
  > tmp/metadata_filtered.tsv

echo "Finding individual GII VP1 types..."
GII_TYPES=()

while IFS= read -r vp1_type
do
  GII_TYPES+=("${vp1_type}")
done < <(
  awk -F'\t' '
    NR > 1 {
      print $2
    }
  ' tmp/metadata_filtered.tsv \
    | sort -t. -k2,2n -u
)


if [[ "${#GII_TYPES[@]}" -eq 0 ]]
then
  echo "ERROR: No GII VP1 types were found." >&2
  exit 1
fi

# Create one FASTA file for every type
echo "Creating FASTA files..."

TYPE_LABELS=()

for vp1_type in "${GII_TYPES[@]}"
do

  type="${vp1_type//./_}"

  TYPE_LABELS+=("${type}")


  awk -F'\t' \
    -v type="${vp1_type}" '
      NR == 1 || $2 == type {
        print
      }
    ' tmp/metadata_filtered.tsv \
    > "tmp/metadata_group_${type}.tsv"


  python scripts/downsample_seqs.py \
    tmp/nextstrain_sequences.fasta \
    "tmp/metadata_group_${type}.tsv" \
    "tmp/sequences_${type}.fasta"

done


# Calculate one VP1-tree medoid per group
echo "Finding VP1-tree medoids..."

python scripts/find_tree_medoids.py \
  --tree tmp/VP1_tree.nwk \
  --metadata tmp/metadata_filtered.tsv \
  --sequences tmp/nextstrain_sequences.fasta \
  --min-genome-length 7000 \
  --group-dir tmp \
  --groups "${TYPE_LABELS[@]}" \
  --output data/references/medoids.tsv

# Download each medoid and extract VP1
echo "Retrieving and extracting medoid VP1 sequences..."

> data/references/reference_multi.fasta

while IFS=$'\t' read -r \
  group \
  medoid_accession \
  vp1_type \
  group_tree_tips \
  full_genome_candidates \
  medoid_genome_length \
  mean_distance

do

  # Skip header.
  if [[ "${group}" == "group" ]]
  then
    continue
  fi

  # Skip empty lines.
  if [[ -z "${group}" ]]
  then
    continue
  fi

  # Remove possible carriage returns.
  group="${group%$'\r'}"
  medoid_accession="${medoid_accession%$'\r'}"
  vp1_type="${vp1_type%$'\r'}"

  echo "Processing ${vp1_type}"
  echo "Medoid accession: ${medoid_accession}"
  echo "Genome length: ${medoid_genome_length}"

  curl -L \
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${medoid_accession}&rettype=gbwithparts&retmode=text" \
    -o "data/references/medoid_${group}.gb"

  python scripts/extract_medoid_vp1.py \
    --genbank "data/references/medoid_${group}.gb" \
    --accession "${medoid_accession}" \
    --group "${group}" \
    --vp1-type "${vp1_type}" \
    --output "data/references/medoid_${group}_VP1.fasta"

  cat "data/references/medoid_${group}_VP1.fasta" \
    >> data/references/reference_multi.fasta


  sleep 1

done < data/references/medoids.tsv

echo
echo "Selected medoids:"
column -t -s $'\t' \
  data/references/medoids.tsv
echo
echo "Final individual-GII VP1 multifasta headers:"
grep "^>" \
  data/references/reference_multi.fasta
echo
echo "Reference statistics:"
seqkit stats \
  data/references/reference_multi.fasta
echo
echo "Number of individual GII references:"
grep -c "^>" \
  data/references/reference_multi.fasta
echo
echo "Medoid workflow completed successfully."