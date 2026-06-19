#!/bin/bash

# CHANGE DIRECTORY TO tmp_consensus/ where tmp_ancestral/ if using consensus sequences, vice versa

set -euo pipefail

THREADS="${THREADS:-8}"

rm -rf barcodeforge_workdir

echo "Creating sapovirus outgroup sequence..."
python scripts/reference_parsing.py \
  --reference data/sapovirus_outgroup.gb \
  --gene VP1 \
  --output data/sapovirus_outgroup_VP1.gb

python scripts/extract_fasta.py \
  data/sapovirus_outgroup_VP1.gb \
  data/sapovirus_outgroup_VP1.fasta

echo "Aligning type-level ancestral sequences..."
augur align \
  --sequences tmp_consensus/type_consensus.fasta \
  --reference-sequence tmp_consensus/norovirus_reference_VP1.gb \
  --output tmp_consensus/alignment.fasta \
  --fill-gaps \
  --remove-reference \
  --nthreads "${THREADS}"

echo "Adding sapovirus outgroup to type-level alignment..."
augur align \
  --existing-alignment tmp_consensus/alignment.fasta \
  --sequences data/sapovirus_outgroup_VP1.fasta \
  --output tmp_consensus/alignment_with_outgroup.fasta \
  --fill-gaps \
  --nthreads "${THREADS}"

echo "Building global tree from type-level ancestors..."
augur tree \
  --alignment tmp_consensus/alignment_with_outgroup.fasta \
  --output tmp_consensus/tree_raw_with_outgroup.nwk \
  --nthreads "${THREADS}" \
  --tree-builder-args="-seed 42"

echo "Rooting/refining global tree..."
augur refine \
  --tree tmp_consensus/tree_raw_with_outgroup.nwk \
  --alignment tmp_consensus/alignment_with_outgroup.fasta \
  --output-tree tmp_consensus/tree_refined_rooted.nwk \
  --output-node-data tmp_consensus/tree_refined_rooted.node_data.json \
  --root NC_075724 \
  --remove-outgroup \
  --branch-length-inference joint \
  --keep-polytomies

echo "Inferring global ancestral sequence..."
augur ancestral \
  --tree tmp_consensus/tree_refined_rooted.nwk \
  --alignment tmp_consensus/alignment.fasta \
  --output-node-data tmp_consensus/ancestral.json \
  --output-sequences tmp_consensus/ancestral_sequences.fasta \
  --inference marginal

python scripts/extract_root_ancestral.py \
  tmp_consensus/tree_refined_rooted.nwk \
  tmp_consensus/ancestral_sequences.fasta \
  tmp_consensus/norovirus_all_VP1_ancestral.fasta

# Create lineage map for BarcodeForge
awk -F'\t' 'NR>1 { print $2 "\t" $1 }' tmp_consensus/metadata_filtered.tsv > tmp_consensus/lineage_map.tsv

# echo "Running BarcodeForge..."
# barcodeforge barcode \
#   tmp/norovirus_all_VP1_ancestral.fasta \
#   tmp/alignment.fasta \
#   tmp/tree_refined_rooted.nwk \
#   tmp/lineage_map.tsv \
#   --tree-format newick \
#   --prefix "NV"