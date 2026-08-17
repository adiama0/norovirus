#!/bin/bash

set -euo pipefail

THREADS="${THREADS:-8}"

# Group A: GII.4, GII.20
# Group B: GII.3, GII.7, GII.8, GII.14
# Group C: GII.6, GII.11, GII.18
# Group D: GII.1, GII.2, GII.5, GII.10, GII.12, GII.16, GII.22, GII.23, GII.24, GII.25, GII.26, GII.27
# Group E: GII.13, GII.17, GII.21

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


# Split metadata into 5 groups

echo "Creating group-specific metadata files..."


# Group A:

awk -F'\t' '
  NR == 1 ||
  $2 ~ /^GII\.(4|20)$/
' tmp/metadata_filtered.tsv \
  > tmp/metadata_group_A.tsv


# Group B:

awk -F'\t' '
  NR == 1 ||
  $2 ~ /^GII\.(3|7|8|14)$/
' tmp/metadata_filtered.tsv \
  > tmp/metadata_group_B.tsv


# Group C:
awk -F'\t' '
  NR == 1 ||
  $2 ~ /^GII\.(6|11|18)$/
' tmp/metadata_filtered.tsv \
  > tmp/metadata_group_C.tsv


# Group D:
awk -F'\t' '
  NR == 1 ||
  $2 ~ /^GII\.(1|2|5|10|12|16|22|23|24|25|26|27)$/
' tmp/metadata_filtered.tsv \
  > tmp/metadata_group_D.tsv


# Group E:

awk -F'\t' '
  NR == 1 ||
  $2 ~ /^GII\.(13|17|21)$/
' tmp/metadata_filtered.tsv \
  > tmp/metadata_group_E.tsv


# Create one FASTA file for every group

echo "Creating group-specific FASTA files..."

for group in A B C D E
do

  python scripts/downsample_seqs.py \
    tmp/nextstrain_sequences.fasta \
    "tmp/metadata_group_${group}.tsv" \
    "tmp/sequences_group_${group}.fasta"

done


# Calculate one VP1-tree medoid per group

echo "Finding VP1-tree medoids..."

python scripts/find_tree_medoids.py \
  --tree tmp/VP1_tree.nwk \
  --metadata tmp/metadata_filtered.tsv \
  --group-dir tmp \
  --groups A B C D E \
  --output data/references/medoids.tsv


# Download each medoid and extract VP1

echo "Retrieving and extracting medoid VP1 sequences..."

> data/references/reference_multi.fasta


while IFS=$'\t' read -r \
  group \
  medoid_accession \
  vp1_type \
  group_tree_tips \
  mean_distance

do

  # Skip TSV header.
  if [[ "${group}" == "group" ]]
  then
    continue
  fi

  # Skip an empty line.
  if [[ -z "${group}" ]]
  then
    continue
  fi

  # Remove possible Windows carriage returns.
  group="${group%$'\r'}"
  medoid_accession="${medoid_accession%$'\r'}"
  vp1_type="${vp1_type%$'\r'}"

  echo
  echo "Processing Group ${group}"
  echo "Medoid accession: ${medoid_accession}"
  echo "VP1 type: ${vp1_type}"


  # Download the GenBank record.
  curl -L \
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${medoid_accession}&rettype=gbwithparts&retmode=text" \
    -o "data/references/medoid_group_${group}.gb"


  # Extract VP1 and assign a clear FASTA header.
  python scripts/extract_medoid_vp1.py \
    --genbank "data/references/medoid_group_${group}.gb" \
    --accession "${medoid_accession}" \
    --group "${group}" \
    --vp1-type "${vp1_type}" \
    --output "data/references/medoid_group_${group}_VP1.fasta"


  # Append this VP1 sequence to the final multifasta.
  cat "data/references/medoid_group_${group}_VP1.fasta" \
    >> data/references/reference_multi.fasta


  # Avoid sending NCBI requests too rapidly.
  sleep 1

done < data/references/medoids.tsv

echo
echo "Medoid workflow completed successfully."
