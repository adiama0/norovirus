#!/bin/bash

set -euo pipefail

rm -rf data/competitive_mapping_results

mkdir -p data/competitive_mapping_results/bams
mkdir -p data/competitive_mapping_results/group_bams


# Define and index the 5-reference VP1 multifasta
REFERENCE="data/references/reference_multi.fasta"
echo "Indexing multireference FASTA..."

bwa-mem2 index "${REFERENCE}"
samtools faidx "${REFERENCE}"

# Retreive reference names
REF_NAMES=()
VP1_TYPES=()

while IFS= read -r ref_name
do
  vp1_type=$(printf "%s\n" "${ref_name}" | cut -d'|' -f2)
  REF_NAMES+=("${ref_name}")
  VP1_TYPES+=("${vp1_type}")
done < <(
  grep '^>' "${REFERENCE}" \
    | sed 's/^>//' \
    | awk '{print $1}'
)

# Create global output tables
COUNTS_FILE="data/competitive_mapping_results/competitive_mapping_counts.tsv"
PROPORTIONS_FILE="data/competitive_mapping_results/competitive_mapping_proportions.tsv"
QC_FILE="data/competitive_mapping_results/competitive_mapping_qc.tsv"

# Counts header.
{
  printf "sample"
  for vp1_type in "${VP1_TYPES[@]}"
  do
    printf "\t%s" "${vp1_type}"
  done
  printf "\tTotal_unique_fragments\n"
} > "${COUNTS_FILE}"

# Proportions header.
{
  printf "sample"
  for vp1_type in "${VP1_TYPES[@]}"
  do
    printf "\t%s" "${vp1_type}"
  done
  printf "\n"
} > "${PROPORTIONS_FILE}"

# QC header.
printf "sample\tInput_read_pairs\tUnique_VP1_fragments\tNot_uniquely_assigned\tAssigned_fraction\n" \
  > "${QC_FILE}"

echo "Processing Samples"
# 33 thru 50 newer samples
for number in {33..50}
do
  directory="data/sequences/sample_${number}"
  sample="sample_${number}"

  echo "Processing ${sample}..."

  R1="${directory}/${sample}_R1.fastq"
  R2="${directory}/${sample}_R2.fastq"
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

  # Counts
  COUNTS=()
  total_unique=0

  for i in "${!REF_NAMES[@]}"
  do
    ref_name="${REF_NAMES[$i]}"
    vp1_type="${VP1_TYPES[$i]}"
    safe_type="${vp1_type//./_}"
    TYPE_BAM="data/competitive_mapping_results/group_bams/${sample}.${safe_type}.unique.bam"

    samtools view \
      -b \
      -f 2 \
      -F 2308 \
      "${BAM}" \
      "${ref_name}" \
      > "${TYPE_BAM}"

    samtools index \
      "${TYPE_BAM}"

    count=$(
      samtools view \
        -c \
        -q 20 \
        -f 66 \
        -F 2308 \
        "${BAM}" \
        "${ref_name}"
    )

    COUNTS+=("${count}")
    total_unique=$((total_unique + count))
  done

  # count original paired fragments 
  input_pairs=$(
    awk 'END {print NR / 4}' "${R1}"
  )

  not_unique=$((input_pairs - total_unique))

  # Normalize genotype counts among VP1 fragments
  PROPORTIONS=()

  if [[ "${total_unique}" -gt 0 ]]
  then
    for count in "${COUNTS[@]}"
    do
      proportion=$(
        awk \
          -v count="${count}" \
          -v total="${total_unique}" \
          'BEGIN {printf "%.6f", count / total}'
      )

      PROPORTIONS+=("${proportion}")
    done

    assigned_fraction=$(
      awk \
        -v assigned="${total_unique}" \
        -v total="${input_pairs}" \
        'BEGIN {printf "%.6f", assigned / total}'
    )

  else
    for vp1_type in "${VP1_TYPES[@]}"
    do
      PROPORTIONS+=("0.000000")
    done

    assigned_fraction="0.000000"
  fi

  # Write raw counts
  {
    printf "%s" "${sample}"
    for count in "${COUNTS[@]}"
    do
      printf "\t%s" "${count}"
    done

    printf "\t%s\n" "${total_unique}"
  } >> "${COUNTS_FILE}"

  # Write normalized proportions
    {
    printf "%s" "${sample}"
    for proportion in "${PROPORTIONS[@]}"
    do
      printf "\t%s" "${proportion}"
    done

    printf "\n"
  } >> "${PROPORTIONS_FILE}"

  # Write QC
  printf "%s\t%s\t%s\t%s\t%s\n" \
    "${sample}" \
    "${input_pairs}" \
    "${total_unique}" \
    "${not_unique}" \
    "${assigned_fraction}" \
    >> "${QC_FILE}"
  
  echo "Competitive mapping proportions:"
  for i in "${!VP1_TYPES[@]}"
  do
    echo "${VP1_TYPES[$i]}: ${PROPORTIONS[$i]}"
  done
done

echo "Competitive mapping completed."