#!/bin/bash

set -euo pipefail

# Group A: GII.4, GII.20
# Group B: GII.3, GII.7, GII.8, GII.14
# Group C: GII.6, GII.11, GII.18
# Group D: GII.1, GII.2, GII.5, GII.10, GII.12, GII.16, GII.22, GII.23, GII.24, GII.25, GII.26, GII.27
# Group E: GII.13, GII.17, GII.21

rm -rf data/competitive_mapping_results

mkdir -p data/competitive_mapping_results/bams
mkdir -p data/competitive_mapping_results/group_bams


# Define and index the 5-reference VP1 multifasta
REFERENCE="data/references/reference_multi.fasta"
echo "Indexing multireference FASTA..."

bwa-mem2 index "${REFERENCE}"
samtools faidx "${REFERENCE}"

# Retreive reference names
REF_A=$(grep -m1 '^>GROUP_A' "${REFERENCE}" | awk '{sub(/^>/,""); print $1}')
REF_B=$(grep -m1 '^>GROUP_B' "${REFERENCE}" | awk '{sub(/^>/,""); print $1}')
REF_C=$(grep -m1 '^>GROUP_C' "${REFERENCE}" | awk '{sub(/^>/,""); print $1}')
REF_D=$(grep -m1 '^>GROUP_D' "${REFERENCE}" | awk '{sub(/^>/,""); print $1}')
REF_E=$(grep -m1 '^>GROUP_E' "${REFERENCE}" | awk '{sub(/^>/,""); print $1}')

# Create global output tables
printf "sample\tGroup_A\tGroup_B\tGroup_C\tGroup_D\tGroup_E\tTotal_unique_fragments\n" \
  > data/competitive_mapping_results/competitive_mapping_counts.tsv

printf "sample\tGroup_A\tGroup_B\tGroup_C\tGroup_D\tGroup_E\n" \
  > data/competitive_mapping_results/competitive_mapping_proportions.tsv

printf "sample\tInput_read_pairs\tUnique_VP1_fragments\tNot_uniquely_assigned\tAssigned_fraction\n" \
  > data/competitive_mapping_results/competitive_mapping_qc.tsv

echo "Processing Samples"
# 33 thru 50 newer samples
for number in {33..50}
do

  directory="data/sequences/sample_${number}"
  sample="sample_${number}"

  echo
  echo "Processing ${sample}..."

  R1="${directory}/${sample}_R1.fastq"
  R2="${directory}/${sample}_R2.fastq"

  if [[ ! -s "${R1}" ]]; then
    echo "ERROR: Missing or empty R1 file: ${R1}" >&2
    exit 1
  fi

  if [[ ! -s "${R2}" ]]; then
    echo "ERROR: Missing or empty R2 file: ${R2}" >&2
    exit 1
  fi


  BAM="data/competitive_mapping_results/bams/${sample}.multireference.sorted.bam"

  # Competitively align each paired fragment against all five references at the same time.
  echo "Aligning ${sample}..."

  bwa-mem2 mem \
    -t 8 \
    -R "@RG\tID:${sample}\tSM:${sample}\tPL:ILLUMINA" \
    "${REFERENCE}" \
    "${R1}" \
    "${R2}" \
    | samtools sort \
        -@ 4 \
        -o "${BAM}"

  samtools index "${BAM}"

  samtools view \
    -b \
    -f 2 \
    -F 2308 \
    "${BAM}" \
    ${REF_A} \
    > "data/competitive_mapping_results/group_bams/${sample}.group_A.unique.bam"

  samtools index \
    "data/competitive_mapping_results/group_bams/${sample}.group_A.unique.bam"

  samtools view \
    -b \
    -f 2 \
    -F 2308 \
    "${BAM}" \
    ${REF_B} \
    > "data/competitive_mapping_results/group_bams/${sample}.group_B.unique.bam"

  samtools index \
    "data/competitive_mapping_results/group_bams/${sample}.group_B.unique.bam"

  samtools view \
    -b \
    -f 2 \
    -F 2308 \
    "${BAM}" \
    ${REF_C} \
    > "data/competitive_mapping_results/group_bams/${sample}.group_C.unique.bam"

  samtools index \
    "data/competitive_mapping_results/group_bams/${sample}.group_C.unique.bam"

  samtools view \
    -b \
    -f 2 \
    -F 2308 \
    "${BAM}" \
    ${REF_D} \
    > "data/competitive_mapping_results/group_bams/${sample}.group_D.unique.bam"

  samtools index \
    "data/competitive_mapping_results/group_bams/${sample}.group_D.unique.bam"

  samtools view \
    -b \
    -f 2 \
    -F 2308 \
    "${BAM}" \
    ${REF_E} \
    > "data/competitive_mapping_results/group_bams/${sample}.group_E.unique.bam"

  samtools index \
    "data/competitive_mapping_results/group_bams/${sample}.group_E.unique.bam"

  # Count uniquely assigned FRAGMENTS.
  count_A=$(
    samtools view \
      -c \
      -q 20 \
      -f 66 \
      -F 2308 \
      "${BAM}" \
      ${REF_A}
  )

  count_B=$(
    samtools view \
      -c \
      -q 20 \
      -f 66 \
      -F 2308 \
      "${BAM}" \
      ${REF_B}
  )

  count_C=$(
    samtools view \
      -c \
      -q 20 \
      -f 66 \
      -F 2308 \
      "${BAM}" \
      ${REF_C}
  )

  count_D=$(
    samtools view \
      -c \
      -q 20 \
      -f 66 \
      -F 2308 \
      "${BAM}" \
      ${REF_D}
  )

  count_E=$(
    samtools view \
      -c \
      -q 20 \
      -f 66 \
      -F 2308 \
      "${BAM}" \
      ${REF_E}
  )

  # Sum all confidently assigned VP1 fragments.
  total_unique=$((count_A + count_B + count_C + count_D + count_E))

  # Count the number of original paired fragments.
  input_pairs=$(
    awk 'END {print NR / 4}' "${R1}"
  )

  not_unique=$((input_pairs - total_unique))

  # Normalize the 5 group counts among confidently assigned
  if [[ "${total_unique}" -gt 0 ]]
  then

    proportion_A=$(
      awk \
        -v count="${count_A}" \
        -v total="${total_unique}" \
        'BEGIN {printf "%.6f", count / total}'
    )

    proportion_B=$(
      awk \
        -v count="${count_B}" \
        -v total="${total_unique}" \
        'BEGIN {printf "%.6f", count / total}'
    )

    proportion_C=$(
      awk \
        -v count="${count_C}" \
        -v total="${total_unique}" \
        'BEGIN {printf "%.6f", count / total}'
    )

    proportion_D=$(
      awk \
        -v count="${count_D}" \
        -v total="${total_unique}" \
        'BEGIN {printf "%.6f", count / total}'
    )

    proportion_E=$(
      awk \
        -v count="${count_E}" \
        -v total="${total_unique}" \
        'BEGIN {printf "%.6f", count / total}'
    )

    assigned_fraction=$(
      awk \
        -v assigned="${total_unique}" \
        -v total="${input_pairs}" \
        'BEGIN {printf "%.6f", assigned / total}'
    )

  else
    proportion_A="0.000000"
    proportion_B="0.000000"
    proportion_C="0.000000"
    proportion_D="0.000000"
    proportion_E="0.000000"
    assigned_fraction="0.000000"

  fi


  # Add raw fragment counts to the global count table.
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "${sample}" \
    "${count_A}" \
    "${count_B}" \
    "${count_C}" \
    "${count_D}" \
    "${count_E}" \
    "${total_unique}" \
    >> data/competitive_mapping_results/competitive_mapping_counts.tsv

  # Add normalized proportions to the global table
  printf "%s\t%s\t%s\t%s\t%s\t%s\n" \
    "${sample}" \
    "${proportion_A}" \
    "${proportion_B}" \
    "${proportion_C}" \
    "${proportion_D}" \
    "${proportion_E}" \
    >> data/competitive_mapping_results/competitive_mapping_proportions.tsv


  # Add overall mapping QC
  printf "%s\t%s\t%s\t%s\t%s\n" \
    "${sample}" \
    "${input_pairs}" \
    "${total_unique}" \
    "${not_unique}" \
    "${assigned_fraction}" \
    >> data/competitive_mapping_results/competitive_mapping_qc.tsv

  echo "Group A: ${proportion_A}"
  echo "Group B: ${proportion_B}"
  echo "Group C: ${proportion_C}"
  echo "Group D: ${proportion_D}"
  echo "Group E: ${proportion_E}"

done

echo "Competitive mapping completed."