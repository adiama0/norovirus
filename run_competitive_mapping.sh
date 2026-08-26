#!/bin/bash

set -euo pipefail

REFERENCE="data/references/reference_multi.fasta"
SEQUENCE_DIR="data/sequences/clinical"

MIN_MAPQ=20
MIN_BASEQ=20

BOOTSTRAPS=10
BOOTSTRAP_SEED_BASE=20260825

# Genotype filtering thresholds
MIN_PROPORTION=0.001
MIN_DEPTH=5
MIN_BREADTH_5X=0.50
MIN_BOOTSTRAP_SUPPORT=0.95
BOOTSTRAP_DETECTION_THRESHOLD="${MIN_PROPORTION}"

RESULTS_DIR="data/competitive_mapping_results"
BAM_DIR="${RESULTS_DIR}/bams"
GROUP_BAM_DIR="${RESULTS_DIR}/group_bams"
CONFIDENCE_DIR="${RESULTS_DIR}/confidence"
BOOTSTRAP_TMP="${CONFIDENCE_DIR}/tmp"

rm -rf "${BAM_DIR}"
rm -rf "${GROUP_BAM_DIR}"
rm -rf "${CONFIDENCE_DIR}"
rm -rf "${RESULTS_DIR}"

mkdir -p "${BAM_DIR}"
mkdir -p "${GROUP_BAM_DIR}"
mkdir -p "${CONFIDENCE_DIR}"
mkdir -p "${BOOTSTRAP_TMP}"

# Index VP1 multifasta
echo "Indexing multireference FASTA..."

bwa-mem2 index "${REFERENCE}"
samtools faidx "${REFERENCE}"

# Retrieve reference names and build a unique VP1 genotype list.
# Expected FASTA header:
# >GROUP_GII_6|GII.6|LC373450|rep1
REF_NAMES=()
REF_TYPES=()
VP1_TYPES=()

while IFS= read -r ref_name
do
  vp1_type=$(printf "%s\n" "${ref_name}" | cut -d'|' -f2)

  if [[ -z "${vp1_type}" || "${vp1_type}" == "${ref_name}" ]]
  then
    echo "ERROR: Reference header does not match expected pipe-delimited format: ${ref_name}" >&2
    echo "Expected: GROUP_VP1_type|VP1.type|accession|rep#" >&2
    exit 1
  fi

  REF_NAMES+=("${ref_name}")
  REF_TYPES+=("${vp1_type}")

done < <(
  grep '^>' "${REFERENCE}" \
    | sed 's/^>//' \
    | awk '{print $1}'
)

while IFS= read -r vp1_type
do
  VP1_TYPES+=("${vp1_type}")
done < <(
  grep '^>' "${REFERENCE}" \
    | sed 's/^>//' \
    | awk '{print $1}' \
    | cut -d'|' -f2 \
    | awk '!seen[$0]++'
)

echo "Detected ${#REF_NAMES[@]} reference sequences representing ${#VP1_TYPES[@]} VP1 genotypes."

for vp1_type in "${VP1_TYPES[@]}"
do
  reference_count=0

  for i in "${!REF_NAMES[@]}"
  do
    if [[ "${REF_TYPES[$i]}" == "${vp1_type}" ]]
    then
      reference_count=$((reference_count + 1))
    fi
  done

  echo "  ${vp1_type}: ${reference_count} reference(s)"
done

# Output files
COUNTS_FILE="${RESULTS_DIR}/competitive_mapping_counts.tsv"
PROPORTIONS_FILE="${RESULTS_DIR}/competitive_mapping_proportions.tsv"
FILTERED_PROPORTIONS_FILE="${RESULTS_DIR}/competitive_mapping_filtered_proportions.tsv"
QC_FILE="${RESULTS_DIR}/competitive_mapping_qc.tsv"
GENOTYPE_QC_FILE="${RESULTS_DIR}/competitive_mapping_genotype_qc.tsv"
REFERENCE_QC_FILE="${RESULTS_DIR}/competitive_mapping_reference_qc.tsv"

# Counts header: one column per unique VP1 genotype
{
  printf "sample"
  for vp1_type in "${VP1_TYPES[@]}"
  do
    printf "\t%s" "${vp1_type}"
  done
  printf "\tTotal_unique_fragments\n"
} > "${COUNTS_FILE}"

# Raw proportions header: one column per unique VP1 genotype
{
  printf "sample"
  for vp1_type in "${VP1_TYPES[@]}"
  do
    printf "\t%s" "${vp1_type}"
  done
  printf "\n"
} > "${PROPORTIONS_FILE}"

# Filtered proportions header: one column per unique VP1 genotype
{
  printf "sample"
  for vp1_type in "${VP1_TYPES[@]}"
  do
    printf "\t%s" "${vp1_type}"
  done
  printf "\n"
} > "${FILTERED_PROPORTIONS_FILE}"

# Sample QC header
printf "sample\tinput_read_pairs\tunique_VP1_fragments\tnot_uniquely_assigned\tassigned_fraction\tmean_MAPQ\tmean_depth\tbreadth\n" \
  > "${QC_FILE}"

# Genotype QC header.
# Coverage metrics come from the best-covered representative for that genotype.
printf "sample\tVP1_type\treference_count\tbest_reference\tfragments\traw_proportion\tmean_depth\tbreadth_1x\tbreadth_5x\tbootstrap_support\tpasses_filter\n" \
  > "${GENOTYPE_QC_FILE}"

# Reference-level diagnostics are retained so we can see which representatives
# are attracting reads and which reference supplied the best coverage metrics.
printf "sample\tVP1_type\treference\tfragments\tmean_depth\tbreadth_1x\tbreadth_5x\n" \
  > "${REFERENCE_QC_FILE}"

echo "Processing samples..."

# Clinical samples
for number in 01 20 21 22 24
do
  directory="${SEQUENCE_DIR}"
  sample="sample_${number}"

  echo
  echo "Processing ${sample}..."

  # Find input files
  R1=""
  R2=""

  for extension in fastq.gz fq.gz fastq fq fasta.gz fa.gz fasta fa
  do
    if [[ -f "${directory}/${sample}_R1.${extension}" ]]
    then
      R1="${directory}/${sample}_R1.${extension}"
      break
    fi
  done

  for extension in fastq.gz fq.gz fastq fq fasta.gz fa.gz fasta fa
  do
    if [[ -f "${directory}/${sample}_R2.${extension}" ]]
    then
      R2="${directory}/${sample}_R2.${extension}"
      break
    fi
  done

  if [[ -z "${R1}" ]]
  then
    echo "ERROR: R1 file not found for ${sample}" >&2
    exit 1
  fi

  if [[ -z "${R2}" ]]
  then
    echo "ERROR: R2 file not found for ${sample}" >&2
    exit 1
  fi

  echo "R1: ${R1}"
  echo "R2: ${R2}"

  # Count input sequences
  input_reads_R1=$(
    seqkit stats -T "${R1}" \
      | awk 'NR == 2 {print $4}'
  )

  input_reads_R2=$(
    seqkit stats -T "${R2}" \
      | awk 'NR == 2 {print $4}'
  )

  if [[ "${input_reads_R1}" -ne "${input_reads_R2}" ]]
  then
    echo "ERROR: R1 and R2 contain different numbers of reads for ${sample}" >&2
    exit 1
  fi

  input_pairs="${input_reads_R1}"
  BAM="${BAM_DIR}/${sample}.multireference.sorted.bam"

  # Competitive alignment against every reference simultaneously
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

  # Genotype-level arrays.
  # Each index corresponds to one unique entry in VP1_TYPES.
  COUNTS=()
  GENOTYPE_MEAN_DEPTHS=()
  BREADTH_1X_VALUES=()
  BREADTH_5X_VALUES=()
  BEST_REFERENCES=()
  REFERENCE_COUNTS=()

  total_unique=0

  # Collapse reference-level assignments into one count per VP1 genotype.
  for vp1_type in "${VP1_TYPES[@]}"
  do
    safe_type="${vp1_type//./_}"
    TYPE_BAM="${GROUP_BAM_DIR}/${sample}.${safe_type}.unique.bam"

    GROUP_REFS=()

    for i in "${!REF_NAMES[@]}"
    do
      if [[ "${REF_TYPES[$i]}" == "${vp1_type}" ]]
      then
        GROUP_REFS+=("${REF_NAMES[$i]}")
      fi
    done

    reference_count="${#GROUP_REFS[@]}"
    REFERENCE_COUNTS+=("${reference_count}")

    # Extract qualifying primary alignments to ANY representative of this genotype.
    # A fragment can have only one retained primary R1 alignment, so summing across 
    # representatives does not create three genotype counts for one fragment.
    samtools view \
      -b \
      -f 2 \
      -q "${MIN_MAPQ}" \
      -F 2308 \
      "${BAM}" \
      "${GROUP_REFS[@]}" \
      > "${TYPE_BAM}"

    samtools index "${TYPE_BAM}"

    # One R1 record per properly paired qualifying fragment.
    genotype_count=$(
      samtools view \
        -c \
        -f 66 \
        "${TYPE_BAM}"
    )

    COUNTS+=("${genotype_count}")
    total_unique=$((total_unique + genotype_count))

    # Coverage is calculated separately for every representative.
    # Retain the representative with the highest breadth at MIN_DEPTH.
    # This avoids penalizing a genotype simply because several alternative references were included in the FASTA.

    best_reference="NA"
    best_mean_depth="0.00"
    best_breadth_1x="0.000000"
    best_breadth_5x="0.000000"
    best_reference_fragments=0

    for ref_name in "${GROUP_REFS[@]}"
    do
      reference_count_for_sample=$(
        samtools view \
          -c \
          -f 66 \
          "${TYPE_BAM}" \
          "${ref_name}"
      )

      reference_coverage=$(
        samtools depth \
          -aa \
          -q "${MIN_BASEQ}" \
          -Q "${MIN_MAPQ}" \
          -s \
          -r "${ref_name}" \
          "${TYPE_BAM}" \
        | awk \
            -v min_depth="${MIN_DEPTH}" '
            {
              positions++
              total_depth += $3

              if ($3 >= 1) {
                covered_1x++
              }

              if ($3 >= min_depth) {
                covered_min_depth++
              }
            }

            END {
              if (positions > 0) {
                mean_depth = total_depth / positions
                breadth_1x = covered_1x / positions
                breadth_min_depth = covered_min_depth / positions

                printf "%.2f\t%.6f\t%.6f",
                  mean_depth,
                  breadth_1x,
                  breadth_min_depth
              } else {
                printf "0.00\t0.000000\t0.000000"
              }
            }
          '
      )

      reference_mean_depth=$(
        printf "%s\n" "${reference_coverage}" \
          | cut -f1
      )

      reference_breadth_1x=$(
        printf "%s\n" "${reference_coverage}" \
          | cut -f2
      )

      reference_breadth_5x=$(
        printf "%s\n" "${reference_coverage}" \
          | cut -f3
      )

      printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "${sample}" \
        "${vp1_type}" \
        "${ref_name}" \
        "${reference_count_for_sample}" \
        "${reference_mean_depth}" \
        "${reference_breadth_1x}" \
        "${reference_breadth_5x}" \
        >> "${REFERENCE_QC_FILE}"

      is_better=$(
        awk \
          -v current_breadth_5x="${reference_breadth_5x}" \
          -v best_breadth_5x="${best_breadth_5x}" \
          -v current_breadth_1x="${reference_breadth_1x}" \
          -v best_breadth_1x="${best_breadth_1x}" \
          -v current_depth="${reference_mean_depth}" \
          -v best_depth="${best_mean_depth}" \
          -v current_fragments="${reference_count_for_sample}" \
          -v best_fragments="${best_reference_fragments}" \
          'BEGIN {
            if (current_breadth_5x > best_breadth_5x) {
              print 1
            } else if (current_breadth_5x == best_breadth_5x && current_breadth_1x > best_breadth_1x) {
              print 1
            } else if (current_breadth_5x == best_breadth_5x && current_breadth_1x == best_breadth_1x && current_depth > best_depth) {
              print 1
            } else if (current_breadth_5x == best_breadth_5x && current_breadth_1x == best_breadth_1x && current_depth == best_depth && current_fragments > best_fragments) {
              print 1
            } else {
              print 0
            }
          }'
      )

      if [[ "${is_better}" -eq 1 ]]
      then
        best_reference="${ref_name}"
        best_mean_depth="${reference_mean_depth}"
        best_breadth_1x="${reference_breadth_1x}"
        best_breadth_5x="${reference_breadth_5x}"
        best_reference_fragments="${reference_count_for_sample}"
      fi
    done

    GENOTYPE_MEAN_DEPTHS+=("${best_mean_depth}")
    BREADTH_1X_VALUES+=("${best_breadth_1x}")
    BREADTH_5X_VALUES+=("${best_breadth_5x}")
    BEST_REFERENCES+=("${best_reference}")
  done

  # Mean MAPQ across primary R1 proper-pair alignments before the MAPQ cutoff
  mean_mapq=$(
    samtools view \
      -f 66 \
      -F 2308 \
      "${BAM}" \
    | awk '
        {
          sum += $5
          n++
        }

        END {
          if (n > 0) {
            printf "%.2f", sum / n
          } else {
            printf "0.00"
          }
        }
      '
  )

  # Sample-level mean depth and breadth
  coverage_stats=$(
    samtools coverage \
      -q "${MIN_MAPQ}" \
      --rf 2 \
      --ff 2308 \
      "${BAM}" \
    | awk '
        $1 !~ /^#/ && $4 > 0 {
          reference_length = $3 - $2 + 1
          total_reference_length += reference_length
          total_covered_bases += $5
          total_depth += $7 * reference_length
        }

        END {
          if (total_reference_length > 0) {
            mean_depth = total_depth / total_reference_length
            breadth = total_covered_bases / total_reference_length
            printf "%.2f\t%.6f", mean_depth, breadth
          } else {
            printf "0.00\t0.000000"
          }
        }
      '
  )

  mean_depth=$(
    printf "%s\n" "${coverage_stats}" \
      | cut -f1
  )

  breadth=$(
    printf "%s\n" "${coverage_stats}" \
      | cut -f2
  )

  not_unique=$((input_pairs - total_unique))

  # Calculate raw genotype proportions from collapsed genotype counts
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

  # Write raw genotype counts
  {
    printf "%s" "${sample}"
    for count in "${COUNTS[@]}"
    do
      printf "\t%s" "${count}"
    done
    printf "\t%s\n" "${total_unique}"
  } >> "${COUNTS_FILE}"

  # Write raw genotype proportions
  {
    printf "%s" "${sample}"
    for proportion in "${PROPORTIONS[@]}"
    do
      printf "\t%s" "${proportion}"
    done
    printf "\n"
  } >> "${PROPORTIONS_FILE}"

  # Write sample QC
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "${sample}" \
    "${input_pairs}" \
    "${total_unique}" \
    "${not_unique}" \
    "${assigned_fraction}" \
    "${mean_mapq}" \
    "${mean_depth}" \
    "${breadth}" \
    >> "${QC_FILE}"

  # Bootstrap analysis
  BOOTSTRAP_PROPORTIONS="${BOOTSTRAP_TMP}/${sample}_bootstrap_proportions.tsv"

  {
    printf "replicate"
    for vp1_type in "${VP1_TYPES[@]}"
    do
      printf "\t%s" "${vp1_type}"
    done
    printf "\n"
  } > "${BOOTSTRAP_PROPORTIONS}"

  # Choose temporary sequence format
  if [[ "${R1}" == *.fastq ]] \
    || [[ "${R1}" == *.fastq.gz ]] \
    || [[ "${R1}" == *.fq ]] \
    || [[ "${R1}" == *.fq.gz ]]
  then
    BOOT_EXT="fastq"
  else
    BOOT_EXT="fasta"
  fi

  sample_number=$((10#${number}))

  for replicate in $(seq 1 "${BOOTSTRAPS}")
  do
    replicate_label=$(
      printf "%03d" "${replicate}"
    )

    seed=$((BOOTSTRAP_SEED_BASE + sample_number * 1000 + replicate))
    BOOT_R1="${BOOTSTRAP_TMP}/${sample}.bootstrap_${replicate_label}_R1.${BOOT_EXT}"
    BOOT_R2="${BOOTSTRAP_TMP}/${sample}.bootstrap_${replicate_label}_R2.${BOOT_EXT}"
    BOOT_BAM="${BOOTSTRAP_TMP}/${sample}.bootstrap_${replicate_label}.bam"

    # Resample paired reads
    python scripts/bootstrap_pairs.py \
      --r1 "${R1}" \
      --r2 "${R2}" \
      --num-pairs "${input_pairs}" \
      --seed "${seed}" \
      --out-r1 "${BOOT_R1}" \
      --out-r2 "${BOOT_R2}"

    # Remap bootstrap replicate
    bwa-mem2 mem \
      -t 8 \
      -R "@RG\tID:${sample}_bootstrap_${replicate_label}\tSM:${sample}\tPL:ILLUMINA" \
      "${REFERENCE}" \
      "${BOOT_R1}" \
      "${BOOT_R2}" \
      | samtools sort \
          -@ 4 \
          -o "${BOOT_BAM}"

    samtools index "${BOOT_BAM}"

    BOOT_COUNTS=()
    bootstrap_total=0

    # Collapse bootstrap counts by genotype exactly as in the original sample.
    for vp1_type in "${VP1_TYPES[@]}"
    do
      GROUP_REFS=()

      for i in "${!REF_NAMES[@]}"
      do
        if [[ "${REF_TYPES[$i]}" == "${vp1_type}" ]]
        then
          GROUP_REFS+=("${REF_NAMES[$i]}")
        fi
      done

      bootstrap_count=$(
        samtools view \
          -c \
          -q "${MIN_MAPQ}" \
          -f 66 \
          -F 2308 \
          "${BOOT_BAM}" \
          "${GROUP_REFS[@]}"
      )

      BOOT_COUNTS+=("${bootstrap_count}")
      bootstrap_total=$((bootstrap_total + bootstrap_count))
    done

    # Write one bootstrap proportion per genotype
    {
      printf "%s" "${replicate}"
      if [[ "${bootstrap_total}" -gt 0 ]]
      then
        for bootstrap_count in "${BOOT_COUNTS[@]}"
        do
          bootstrap_proportion=$(
            awk \
              -v count="${bootstrap_count}" \
              -v total="${bootstrap_total}" \
              'BEGIN {printf "%.6f", count / total}'
          )
          printf "\t%s" "${bootstrap_proportion}"
        done
      else
        for vp1_type in "${VP1_TYPES[@]}"
        do
          printf "\t0.000000"
        done
      fi
      printf "\n"
    } >> "${BOOTSTRAP_PROPORTIONS}"

    # Remove bootstrap BAM and sequence files
    rm -f "${BOOT_R1}"
    rm -f "${BOOT_R2}"
    rm -f "${BOOT_BAM}"
    rm -f "${BOOT_BAM}.bai"

    if [[ "${replicate}" -eq 1 ]] \
      || [[ $((replicate % 10)) -eq 0 ]]
    then
      echo "Completed bootstrap ${replicate}/${BOOTSTRAPS}"
    fi
  done

  # Summarize bootstrap confidence.
  # summarize_bootstrap_confidence.py does not need to change because the
  # bootstrap table now contains one unique column per genotype.
  CONFIDENCE_FILE="${CONFIDENCE_DIR}/${sample}_confidence.tsv"

  python scripts/summarize_bootstrap_confidence.py \
    --proportions "${PROPORTIONS_FILE}" \
    --bootstraps "${BOOTSTRAP_PROPORTIONS}" \
    --sample "${sample}" \
    --detection-threshold "${BOOTSTRAP_DETECTION_THRESHOLD}" \
    --output "${CONFIDENCE_FILE}"

  # Apply combined genotype filter
  FILTERED_COUNTS=()
  PASS_FLAGS=()
  BOOTSTRAP_SUPPORT_VALUES=()

  filtered_total=0

  for i in "${!VP1_TYPES[@]}"
  do
    vp1_type="${VP1_TYPES[$i]}"
    raw_proportion="${PROPORTIONS[$i]}"
    breadth_5x="${BREADTH_5X_VALUES[$i]}"

    bootstrap_support_percent=$(
      awk \
        -F '\t' \
        -v type="${vp1_type}" '
        NR > 1 && $1 == type {
          gsub("%", "", $4)
          print $4
          found = 1
        }

        END {
          if (!found) {
            print 0
          }
        }
        ' \
        "${CONFIDENCE_FILE}"
    )

    bootstrap_support=$(
      awk \
        -v support="${bootstrap_support_percent}" \
        'BEGIN {printf "%.6f", support / 100}'
    )

    BOOTSTRAP_SUPPORT_VALUES+=("${bootstrap_support}")

    passes_filter=$(
      awk \
        -v proportion="${raw_proportion}" \
        -v breadth="${breadth_5x}" \
        -v support="${bootstrap_support}" \
        -v min_proportion="${MIN_PROPORTION}" \
        -v min_breadth="${MIN_BREADTH_5X}" \
        -v min_support="${MIN_BOOTSTRAP_SUPPORT}" \
        'BEGIN {
          if (proportion >= min_proportion && breadth >= min_breadth && support >= min_support) {
            print 1
          } else {
            print 0
          }
        }'
    )

    PASS_FLAGS+=("${passes_filter}")

    if [[ "${passes_filter}" -eq 1 ]]
    then
      FILTERED_COUNTS+=("${COUNTS[$i]}")
      filtered_total=$((filtered_total + COUNTS[$i]))
    else
      FILTERED_COUNTS+=("0")
    fi
  done

  # Renormalize only passing genotypes
  FILTERED_PROPORTIONS=()

  if [[ "${filtered_total}" -gt 0 ]]
  then
    for filtered_count in "${FILTERED_COUNTS[@]}"
    do
      filtered_proportion=$(
        awk \
          -v count="${filtered_count}" \
          -v total="${filtered_total}" \
          'BEGIN {printf "%.6f", count / total}'
      )
      FILTERED_PROPORTIONS+=("${filtered_proportion}")
    done
  else
    for vp1_type in "${VP1_TYPES[@]}"
    do
      FILTERED_PROPORTIONS+=("0.000000")
    done
  fi

  # Write filtered genotype proportions
  {
    printf "%s" "${sample}"
    for filtered_proportion in "${FILTERED_PROPORTIONS[@]}"
    do
      printf "\t%s" "${filtered_proportion}"
    done
    printf "\n"
  } >> "${FILTERED_PROPORTIONS_FILE}"

  # Write genotype-level QC
  for i in "${!VP1_TYPES[@]}"
  do
    if [[ "${PASS_FLAGS[$i]}" -eq 1 ]]
    then
      filter_result="PASS"
    else
      filter_result="FAIL"
    fi

    bootstrap_support_percent=$(
      awk \
        -v support="${BOOTSTRAP_SUPPORT_VALUES[$i]}" \
        'BEGIN {printf "%.0f%%", support * 100}'
    )

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
      "${sample}" \
      "${VP1_TYPES[$i]}" \
      "${REFERENCE_COUNTS[$i]}" \
      "${BEST_REFERENCES[$i]}" \
      "${COUNTS[$i]}" \
      "${PROPORTIONS[$i]}" \
      "${GENOTYPE_MEAN_DEPTHS[$i]}" \
      "${BREADTH_1X_VALUES[$i]}" \
      "${BREADTH_5X_VALUES[$i]}" \
      "${bootstrap_support_percent}" \
      "${filter_result}" \
      >> "${GENOTYPE_QC_FILE}"
  done
  for i in "${!VP1_TYPES[@]}"
  do
    echo "${VP1_TYPES[$i]}: ${FILTERED_PROPORTIONS[$i]}"
  done
done

# Remove bootstrap temporary files
rm -rf "${BOOTSTRAP_TMP}"

echo "Competitive mapping completed."