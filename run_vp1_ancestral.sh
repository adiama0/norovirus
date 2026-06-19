#!/bin/bash

set -euo pipefail

THREADS="${THREADS:-8}"

# Uncomment to rerun full workflow
rm -rf tmp_ancestral barcodeforge_workdir
rm -f data/sequences.fasta data/sequences.fasta.zst

mkdir -p tmp_ancestral
mkdir -p tmp_ancestral/by_type

# Pull sequences from Nextstrain
wget https://data.nextstrain.org/files/workflows/norovirus/sequences.fasta.zst \
  -O data/sequences.fasta.zst && zstd -d data/sequences.fasta.zst

# Extract metadata from the Auspice JSON
echo "Extracting metadata..."
barcodeforge extract-auspice-data \
  data/norovirus_all_VP1.json \
  --attributes "VP1_type" \
  --attributes "region" \
  --output_metadata_path tmp_ancestral/metadata.tsv

# Keep only North America samples with a valid VP1 type
echo "Cleaning metadata..."
awk -F'\t' 'NR==1 || ($2 != "")' \
  tmp_ancestral/metadata.tsv > tmp_ancestral/metadata_filtered.tsv

# Downsample the raw sequence file to only the samples in the metadata
echo "Creating FASTA file..."
python scripts/downsample_seqs.py \
  data/sequences.fasta \
  tmp_ancestral/metadata_filtered.tsv \
  tmp_ancestral/downsample_seqs.fasta

# Create the norovirus VP1 reference GenBank file
echo "Creating VP1 reference..."
python scripts/reference_parsing.py \
  --reference data/reference.gb \
  --gene VP1 \
  --output tmp_ancestral/norovirus_reference_VP1.gb

# Split sequences and metadata by VP1 type
echo "Splitting sequences by VP1 type..."
python scripts/split_by_vp1_type.py \
  tmp_ancestral/metadata_filtered.tsv \
  tmp_ancestral/downsample_seqs.fasta \
  tmp_ancestral/by_type

# Infer one ancestral sequence per VP1 type
echo "Inferring ancestral sequence per VP1 type..."
for TYPE_DIR in tmp_ancestral/by_type/*; do
  [ -d "${TYPE_DIR}" ] || continue

  VP1_TYPE=$(cat "${TYPE_DIR}/type_name.txt")
  echo "Processing ${VP1_TYPE}"

  NSEQ=$(grep -c '^>' "${TYPE_DIR}/sequences.fasta" || true)

  if [ "${NSEQ}" -eq 0 ]; then
    echo "Skipping ${VP1_TYPE}: no sequences"
    continue
  fi

  if [ "${NSEQ}" -eq 1 ]; then
    echo "Single-sequence type for ${VP1_TYPE}; extracting VP1-aligned region only"

    augur align \
    --sequences "${TYPE_DIR}/sequences.fasta" \
    --reference-sequence tmp_ancestral/norovirus_reference_VP1.gb \
    --output "${TYPE_DIR}/single_alignment.fasta" \
    --fill-gaps \
    --remove-reference \
    --nthreads "${THREADS}"

    python scripts/rename_single_fasta.py \
    "${TYPE_DIR}/single_alignment.fasta" \
    "${TYPE_DIR}/type_ancestral.fasta" \
    "${VP1_TYPE}"
    continue
  fi

  if [ "${NSEQ}" -lt 3 ]; then
    echo "Skipping ${VP1_TYPE}: not enough sequences to build a tree (${NSEQ})"
    continue
  fi

  # Align within this VP1 type
  augur align \
    --sequences "${TYPE_DIR}/sequences.fasta" \
    --reference-sequence tmp_ancestral/norovirus_reference_VP1.gb \
    --output "${TYPE_DIR}/alignment.fasta" \
    --fill-gaps \
    --remove-reference \
    --nthreads "${THREADS}"

  # Build tree within this VP1 type
  augur tree \
    --alignment "${TYPE_DIR}/alignment.fasta" \
    --output "${TYPE_DIR}/tree_raw.nwk" \
    --nthreads "${THREADS}"

  # Refine within type
  augur refine \
    --tree "${TYPE_DIR}/tree_raw.nwk" \
    --alignment "${TYPE_DIR}/alignment.fasta" \
    --output-tree "${TYPE_DIR}/tree_refined.nwk" \
    --output-node-data "${TYPE_DIR}/tree_refined.node_data.json" \
    --root mid_point \
    --branch-length-inference joint \
    --keep-polytomies

  # Infer ancestral sequences within type
  augur ancestral \
    --tree "${TYPE_DIR}/tree_refined.nwk" \
    --alignment "${TYPE_DIR}/alignment.fasta" \
    --output-node-data "${TYPE_DIR}/ancestral.json" \
    --output-sequences "${TYPE_DIR}/ancestral_sequences.fasta" \
    --inference marginal

  # Extract the root ancestral sequence for this type
  python scripts/extract_root_ancestral.py \
    "${TYPE_DIR}/tree_refined.nwk" \
    "${TYPE_DIR}/ancestral_sequences.fasta" \
    "${TYPE_DIR}/type_ancestral_raw.fasta"

  # Rename to the VP1 type
  python scripts/rename_single_fasta.py \
    "${TYPE_DIR}/type_ancestral_raw.fasta" \
    "${TYPE_DIR}/type_ancestral.fasta" \
    "${VP1_TYPE}"
done

# Collect all type-level consensus sequences into one FASTA
echo "Collecting type-level consensus sequences..."
find tmp_ancestral/by_type -name "type_ancestral.fasta" -exec cat {} + > tmp_ancestral/type_ancestral.fasta