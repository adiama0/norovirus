# !/bin/bash

set -euo pipefail

THREADS="${THREADS:-8}"

# uncomment to rerun full workflow
rm -rf tmp barcodeforge_workdir
rm -f data/sequences.fasta data/sequences.fasta.zst

mkdir -p tmp

# Pull sequences from Nextstrain
wget https://data.nextstrain.org/files/workflows/norovirus/sequences.fasta.zst \
  -O data/sequences.fasta.zst && zstd -d data/sequences.fasta.zst

# Extract the tree and metadata from the auspice JSON
echo "Extracting metadata..."
barcodeforge extract-auspice-data \
  data/norovirus_all_VP1.json \
  --attributes "VP1_type" \
  --output_metadata_path tmp/metadata.tsv

# Remove unlabeled samples to avoid noise
echo "Cleaning metadata.."
awk -F'\t' 'NR==1 || $2 != ""' tmp/metadata.tsv > tmp/metadata_filtered.tsv

# Downsample the sequences from Nextstrain to the sequences present in the metadata file
echo "Creating fasta file..."
python downsample_seqs.py data/sequences.fasta tmp/metadata_filtered.tsv tmp/downsample_seqs.fasta

# Create the VP1 reference GenBank file using Nextstrain script
echo "Creating VP1 reference..."
python reference_parsing.py \
  --reference data/reference.gb \
  --gene VP1 \
  --output tmp/norovirus_reference_all_VP1.gb

# Align sequences with augur align, matching the Nextstrain strategy
echo "Aligning sequences..."
augur align \
  --sequences tmp/downsample_seqs.fasta \
  --reference-sequence tmp/norovirus_reference_all_VP1.gb \
  --output tmp/alignment.fasta \
  --fill-gaps \
  --remove-reference \
  --nthreads "${THREADS}"

# Build a tree from the same alignment
echo "Building tree..."
augur tree \
  --alignment tmp/alignment.fasta \
  --output tmp/tree_raw.nwk \
  --nthreads "${THREADS}"

echo "Rooting/refining tree..."
augur refine \
  --tree tmp/tree_raw.nwk \
  --alignment tmp/alignment.fasta \
  --output-tree tmp/tree_refined.nwk \
  --output-node-data tmp/tree_refined.node_data.json \
  --root mid_point

echo "Inferring ancestral sequences..."
augur ancestral \
  --tree tmp/tree_refined.nwk \
  --alignment tmp/alignment.fasta \
  --output-node-data tmp/ancestral.json \
  --output-sequences tmp/ancestral_sequences.fasta \
  --inference joint

python extract_root_ancestral.py \
  tmp/tree_refined.nwk \
  tmp/ancestral_sequences.fasta \
  tmp/norovirus_all_VP1_ancestral.fasta

# Create lineage map for BarcodeForge
awk -F'\t' 'NR>1{print $2"\t"$1}' tmp/metadata_filtered.tsv > tmp/lineage_map.tsv

# python extract_fasta.py tmp/norovirus_reference_all_VP1.gb tmp/norovirus_reference_all_VP1.fasta

# Run BarcodeForge 
echo "Running Barcodeforge..."
barcodeforge barcode \
  tmp/norovirus_all_VP1_ancestral.fasta \
  tmp/alignment.fasta \
  tmp/tree_refined.nwk \
  tmp/lineage_map.tsv \
  --tree-format newick \
  --prefix "NV"