#!/bin/bash

set -euo pipefail

THREADS="${THREADS:-8}"

# Uncomment to rerun full workflow
rm -rf tmp_consensus 
rm -f data/sequences.fasta data/sequences.fasta.zst

mkdir -p tmp_consensus
mkdir -p tmp_consensus/by_type

# Pull sequences from Nextstrain
wget https://data.nextstrain.org/files/workflows/norovirus/sequences.fasta.zst \
  -O data/sequences.fasta.zst && zstd -d data/sequences.fasta.zst

# Extract metadata from the Auspice JSON
echo "Extracting metadata..."
barcodeforge extract-auspice-data \
  data/norovirus_all_VP1.json \
  --attributes "VP1_type" \
  --attributes "region" \
  --output_metadata_path tmp_consensus/metadata.tsv

# Keep only North America samples with a valid VP1 type
echo "Cleaning metadata..."
awk -F'\t' 'NR==1 || ($2 != "")' \
  tmp_consensus/metadata.tsv > tmp_consensus/metadata_filtered.tsv

# Downsample the raw sequence file to only the samples in the metadata
echo "Creating FASTA file..."
python scripts/downsample_seqs.py \
  data/sequences.fasta \
  tmp_consensus/metadata_filtered.tsv \
  tmp_consensus/downsample_seqs.fasta

# Create the norovirus VP1 reference GenBank file
echo "Creating VP1 reference..."
python scripts/reference_parsing.py \
  --reference data/reference.gb \
  --gene VP1 \
  --output tmp_consensus/norovirus_reference_VP1.gb

# Split sequences and metadata by VP1 type
echo "Splitting sequences by VP1 type..."
python scripts/split_by_vp1_type.py \
  tmp_consensus/metadata_filtered.tsv \
  tmp_consensus/downsample_seqs.fasta \
  tmp_consensus/by_type

# Build one consensus sequence per VP1 type
echo "Building consensus sequence per VP1 type..."
for TYPE_DIR in tmp_consensus/by_type/*; do
  [ -d "${TYPE_DIR}" ] || continue

  VP1_TYPE=$(cat "${TYPE_DIR}/type_name.txt")
  echo "Processing ${VP1_TYPE}"

  NSEQ=$(grep -c '^>' "${TYPE_DIR}/sequences.fasta" || true)

  if [ "${NSEQ}" -eq 0 ]; then
    echo "Skipping ${VP1_TYPE}: no sequences"
    continue
  fi

  # Align sequences within this VP1 type to VP1 coordinate space
  augur align \
    --sequences "${TYPE_DIR}/sequences.fasta" \
    --reference-sequence tmp_consensus/norovirus_reference_VP1.gb \
    --output "${TYPE_DIR}/alignment.fasta" \
    --fill-gaps \
    --remove-reference \
    --nthreads "${THREADS}"

  # Build a consensus from the alignment
  python scripts/make_consensus.py \
    "${TYPE_DIR}/alignment.fasta" \
    "${TYPE_DIR}/type_consensus.fasta" \
    "${VP1_TYPE}"

done

# Collect all type-level consensus sequences into one FASTA
echo "Collecting type-level consensus sequences..."
find tmp_consensus/by_type -name "type_consensus.fasta" -exec cat {} + > tmp_consensus/type_consensus.fasta