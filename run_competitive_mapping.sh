#!/bin/bash

set -uo pipefail

# configuration
REFERENCE="data/references/reference_multi.fasta"
SEQUENCE_DIR="data/sequences/clinical/"
KRAKEN2_DB="${KRAKEN2_DB:-data/kraken2}"

THREADS=8
SORT_THREADS=4
FASTQC_THREADS=4

MIN_MAPQ=20
MIN_BASEQ=20
MIN_PROPORTION=0.001
MIN_DEPTH=5
MIN_BREADTH_5X=0.50

# Kraken2 used for metagenomic screen before VP1 genotyping.
# Confidence is a Kraken2 k-mer support threshold, not a probability.
KRAKEN2_CONFIDENCE="${KRAKEN2_CONFIDENCE:-0.10}"
NOROVIRUS_TAXID=142786

TRIMMOMATIC_ADAPTERS="${TRIMMOMATIC_ADAPTERS:-}"
TRIMMOMATIC_SLIDING_WINDOW="4:20"
TRIMMOMATIC_MINLEN=50

# output directories
RESULTS_DIR="data/results"
FASTQC_DIR="${RESULTS_DIR}/fastqc"
TRIMMOMATIC_DIR="${RESULTS_DIR}/trimmomatic"
TAXONOMY_DIR="${RESULTS_DIR}/taxonomy"
COMPETITIVE_DIR="${RESULTS_DIR}/competitive_mapping"
TMP_DIR="${RESULTS_DIR}/tmp"

DETECTION_FILE="${TAXONOMY_DIR}/norovirus_detection.tsv"
FILTERED_PROPORTIONS_FILE="${COMPETITIVE_DIR}/mapping_proportions_filtered.tsv"
QC_FILE="${COMPETITIVE_DIR}/mapping_qc.tsv"
REFERENCE_QC_FILE="${COMPETITIVE_DIR}/mapping_reference_qc.tsv"

# Start each run with one clean results directory.
rm -rf "${RESULTS_DIR}"
mkdir -p \
  "${FASTQC_DIR}" \
  "${TRIMMOMATIC_DIR}" \
  "${TAXONOMY_DIR}" \
  "${COMPETITIVE_DIR}" \
  "${TMP_DIR}"

# check required programs
for program in bwa-mem2 samtools seqkit fastqc trimmomatic kraken2
do
  if ! command -v "${program}" >/dev/null 2>&1
  then
    echo "ERROR: ${program} is not available in the current environment." >&2
    exit 1
  fi
done

if [[ ! -f "${REFERENCE}" ]]
then
  echo "ERROR: VP1 reference FASTA not found: ${REFERENCE}" >&2
  exit 1
fi

if [[ ! -d "${KRAKEN2_DB}" ]]
then
  echo "ERROR: Kraken2 database not found: ${KRAKEN2_DB}" >&2
  echo "Set KRAKEN2_DB to the database directory before running." >&2
  exit 1
fi

# index and parse VP1 reference FASTA
if ! bwa-mem2 index "${REFERENCE}"
then
  echo "ERROR: bwa-mem2 could not index ${REFERENCE}." >&2
  exit 1
fi

if ! samtools faidx "${REFERENCE}"
then
  echo "ERROR: samtools could not index ${REFERENCE}." >&2
  exit 1
fi

REF_NAMES=()
REF_TYPES=()
VP1_TYPES=()

while IFS= read -r ref_name
do
  vp1_type=$(printf "%s\n" "${ref_name}" | cut -d'|' -f2)

  if [[ -z "${vp1_type}" || "${vp1_type}" == "${ref_name}" ]]
  then
    echo "ERROR: Invalid reference header: ${ref_name}" >&2
    echo "Expected: GROUP_VP1_type|VP1.type|accession|rep#" >&2
    exit 1
  fi

  REF_NAMES+=("${ref_name}")
  REF_TYPES+=("${vp1_type}")
done < <(
  grep '^>' "${REFERENCE}" |
    sed 's/^>//' |
    awk '{print $1}'
)

while IFS= read -r vp1_type
do
  VP1_TYPES+=("${vp1_type}")
done < <(
  grep '^>' "${REFERENCE}" |
    sed 's/^>//' |
    awk '{print $1}' |
    cut -d'|' -f2 |
    awk '!seen[$0]++'
)

if [[ "${#REF_NAMES[@]}" -eq 0 ]]
then
  echo "ERROR: No references found in ${REFERENCE}." >&2
  exit 1
fi

# initialize outputs
{
  printf "sample"
  for vp1_type in "${VP1_TYPES[@]}"
  do
    printf "\t%s" "${vp1_type}"
  done
  printf "\n"
} > "${FILTERED_PROPORTIONS_FILE}"

printf "sample\traw_read_pairs\tmapping_read_pairs\tretained_pair_fraction\tunique_VP1_fragments\tassigned_fraction\tmean_MAPQ\tmean_depth\tbreadth\n" \
  > "${QC_FILE}"

printf "sample\tVP1_type\treference\tfragments\tmean_depth\tbreadth_1x\tbreadth_5x\n" \
  > "${REFERENCE_QC_FILE}"

printf "sample\tkraken2_norovirus_percent\tkraken2_clade_fragments\tkraken2_direct_fragments\tkraken2_total_minimizers\tkraken2_distinct_minimizers\ttaxonomic_screen\tfinal_detection\n" \
  > "${DETECTION_FILE}"

# count sequences in FASTQ or FASTA
count_sequences() {
  seqkit stats -T "$1" | awk 'NR == 2 {print $4}'
}

# read samples with one strict naming convention
read_samples() {
  local sample_file
  local file
  local filename
  local sample

  SAMPLES=()
  sample_file=$(mktemp "${TMP_DIR}/samples.XXXXXX")

  for file in "${SEQUENCE_DIR}"/*_R1.fastq.gz "${SEQUENCE_DIR}"/*_R1.fasta
  do
    [[ -e "${file}" ]] || continue
    filename=$(basename "${file}")

    case "${filename}" in
      *_R1.fastq.gz) sample="${filename%_R1.fastq.gz}" ;;
      *_R1.fasta) sample="${filename%_R1.fasta}" ;;
      *) continue ;;
    esac

    printf "%s\n" "${sample}" >> "${sample_file}"
  done

  if [[ ! -s "${sample_file}" ]]
  then
    rm -f "${sample_file}"
    echo "ERROR: No SAMPLE_R1.fastq.gz or SAMPLE_R1.fasta files found in ${SEQUENCE_DIR}." >&2
    exit 1
  fi

  while IFS= read -r sample
  do
    SAMPLES+=("${sample}")
  done < <(sort -u "${sample_file}")

  rm -f "${sample_file}"
}

# find the paired files for one sample
find_sample_files() {
  local sample="$1"
  local fastq_r1="${SEQUENCE_DIR}/${sample}_R1.fastq.gz"
  local fastq_r2="${SEQUENCE_DIR}/${sample}_R2.fastq.gz"
  local fasta_r1="${SEQUENCE_DIR}/${sample}_R1.fasta"
  local fasta_r2="${SEQUENCE_DIR}/${sample}_R2.fasta"

  R1=""
  R2=""
  INPUT_FORMAT=""

  if [[ -f "${fastq_r1}" || -f "${fastq_r2}" ]]
  then
    if [[ ! -f "${fastq_r1}" || ! -f "${fastq_r2}" ]]
    then
      echo "Skipping ${sample}: FASTQ R1/R2 pair is incomplete." >&2
      return 1
    fi

    if [[ -f "${fasta_r1}" || -f "${fasta_r2}" ]]
    then
      echo "Skipping ${sample}: both FASTQ and FASTA inputs exist for the same sample." >&2
      return 1
    fi

    R1="${fastq_r1}"
    R2="${fastq_r2}"
    INPUT_FORMAT="fastq"
    return 0
  fi

  if [[ -f "${fasta_r1}" || -f "${fasta_r2}" ]]
  then
    if [[ ! -f "${fasta_r1}" || ! -f "${fasta_r2}" ]]
    then
      echo "Skipping ${sample}: FASTA R1/R2 pair is incomplete." >&2
      return 1
    fi

    R1="${fasta_r1}"
    R2="${fasta_r2}"
    INPUT_FORMAT="fasta"
    return 0
  fi

  echo "Skipping ${sample}: matching R1/R2 files were not found." >&2
  return 1
}

# run FastQC on raw FASTQ reads
run_fastqc() {
  local sample="$1"

  echo "Running FastQC for ${sample}..."

  if ! fastqc \
    --quiet \
    --threads "${FASTQC_THREADS}" \
    --outdir "${FASTQC_DIR}" \
    "${R1}" \
    "${R2}"
  then
    echo "Skipping ${sample}: FastQC could not analyze the input files." >&2
    return 1
  fi

  return 0
}

# trim paired FASTQ reads; only paired survivors are used downstream
trim_sample() {
  local sample="$1"
  local unpaired_r1="${TMP_DIR}/${sample}_R1.trimmed.unpaired.fastq.gz"
  local unpaired_r2="${TMP_DIR}/${sample}_R2.trimmed.unpaired.fastq.gz"
  local log_file="${TRIMMOMATIC_DIR}/${sample}.trimmomatic.log"

  TRIMMED_R1="${TRIMMOMATIC_DIR}/${sample}_R1.trimmed.paired.fastq.gz"
  TRIMMED_R2="${TRIMMOMATIC_DIR}/${sample}_R2.trimmed.paired.fastq.gz"

  echo "Trimming ${sample}..."

  if [[ -n "${TRIMMOMATIC_ADAPTERS}" ]]
  then
    if [[ ! -f "${TRIMMOMATIC_ADAPTERS}" ]]
    then
      echo "Skipping ${sample}: adapter file not found: ${TRIMMOMATIC_ADAPTERS}" >&2
      return 1
    fi

    if ! trimmomatic PE \
      -threads "${THREADS}" \
      -phred33 \
      "${R1}" \
      "${R2}" \
      "${TRIMMED_R1}" \
      "${unpaired_r1}" \
      "${TRIMMED_R2}" \
      "${unpaired_r2}" \
      "ILLUMINACLIP:${TRIMMOMATIC_ADAPTERS}:2:30:10" \
      "SLIDINGWINDOW:${TRIMMOMATIC_SLIDING_WINDOW}" \
      "MINLEN:${TRIMMOMATIC_MINLEN}" \
      > "${log_file}" 2>&1
    then
      echo "Skipping ${sample}: Trimmomatic failed. See ${log_file}." >&2
      return 1
    fi
  else
    if ! trimmomatic PE \
      -threads "${THREADS}" \
      -phred33 \
      "${R1}" \
      "${R2}" \
      "${TRIMMED_R1}" \
      "${unpaired_r1}" \
      "${TRIMMED_R2}" \
      "${unpaired_r2}" \
      "SLIDINGWINDOW:${TRIMMOMATIC_SLIDING_WINDOW}" \
      "MINLEN:${TRIMMOMATIC_MINLEN}" \
      > "${log_file}" 2>&1
    then
      echo "Skipping ${sample}: Trimmomatic failed. See ${log_file}." >&2
      return 1
    fi
  fi

  rm -f "${unpaired_r1}" "${unpaired_r2}"
  return 0
}

# prepare the paired reads that enter taxonomy and mapping
prepare_reads() {
  local sample="$1"
  local raw_r1_count
  local raw_r2_count
  local mapping_r1_count
  local mapping_r2_count

  raw_r1_count=$(count_sequences "${R1}")
  raw_r2_count=$(count_sequences "${R2}")

  if [[ -z "${raw_r1_count}" || -z "${raw_r2_count}" || "${raw_r1_count}" -eq 0 || "${raw_r2_count}" -eq 0 ]]
  then
    echo "Skipping ${sample}: input files are empty or unreadable." >&2
    return 1
  fi

  if [[ "${raw_r1_count}" -ne "${raw_r2_count}" ]]
  then
    echo "Skipping ${sample}: R1 and R2 contain different numbers of reads." >&2
    return 1
  fi

  RAW_INPUT_PAIRS="${raw_r1_count}"

  if [[ "${INPUT_FORMAT}" == "fastq" ]]
  then
    if ! run_fastqc "${sample}"
    then
      return 1
    fi

    if ! trim_sample "${sample}"
    then
      return 1
    fi

    MAPPING_R1="${TRIMMED_R1}"
    MAPPING_R2="${TRIMMED_R2}"
  else
    echo "FASTA input detected for ${sample}; skipping FastQC and Trimmomatic."
    MAPPING_R1="${R1}"
    MAPPING_R2="${R2}"
  fi

  mapping_r1_count=$(count_sequences "${MAPPING_R1}")
  mapping_r2_count=$(count_sequences "${MAPPING_R2}")

  if [[ -z "${mapping_r1_count}" || -z "${mapping_r2_count}" || "${mapping_r1_count}" -eq 0 || "${mapping_r2_count}" -eq 0 ]]
  then
    echo "Skipping ${sample}: no paired reads remain for analysis." >&2
    return 1
  fi

  if [[ "${mapping_r1_count}" -ne "${mapping_r2_count}" ]]
  then
    echo "Skipping ${sample}: processed R1 and R2 contain different numbers of reads." >&2
    return 1
  fi

  MAPPING_INPUT_PAIRS="${mapping_r1_count}"
  RETAINED_PAIR_FRACTION=$(
    awk \
      -v retained="${MAPPING_INPUT_PAIRS}" \
      -v raw="${RAW_INPUT_PAIRS}" \
      'BEGIN {printf "%.6f", retained / raw}'
  )

  return 0
}

# screen the metagenome for Norovirus before genotype mapping
screen_norovirus() {
  local sample="$1"
  local report_file="${TAXONOMY_DIR}/${sample}.kraken2.report.tsv"
  local metrics

  echo "Running Kraken2 taxonomic screen for ${sample}..."

  if ! kraken2 \
    --db "${KRAKEN2_DB}" \
    --threads "${THREADS}" \
    --paired \
    --confidence "${KRAKEN2_CONFIDENCE}" \
    --report "${report_file}" \
    --report-minimizer-data \
    --memory-mapping \
    --output - \
    "${MAPPING_R1}" \
    "${MAPPING_R2}"
  then
    echo "Skipping ${sample}: Kraken2 taxonomic classification failed." >&2
    return 1
  fi

  metrics=$(
    awk \
      -F '\t' \
      -v taxid="${NOROVIRUS_TAXID}" '
        $7 == taxid {
          printf "%s\t%s\t%s\t%s\t%s", $1, $2, $3, $4, $5
          found = 1
        }
        END {
          if (!found) {
            printf "0.00\t0\t0\t0\t0"
          }
        }' \
      "${report_file}"
  )

  KRAKEN_NORO_PERCENT=$(printf "%s\n" "${metrics}" | cut -f1)
  KRAKEN_NORO_CLADE_FRAGMENTS=$(printf "%s\n" "${metrics}" | cut -f2)
  KRAKEN_NORO_DIRECT_FRAGMENTS=$(printf "%s\n" "${metrics}" | cut -f3)
  KRAKEN_NORO_TOTAL_MINIMIZERS=$(printf "%s\n" "${metrics}" | cut -f4)
  KRAKEN_NORO_DISTINCT_MINIMIZERS=$(printf "%s\n" "${metrics}" | cut -f5)

  if [[ "${KRAKEN_NORO_CLADE_FRAGMENTS}" -gt 0 && "${KRAKEN_NORO_DISTINCT_MINIMIZERS}" -gt 0 ]]
  then
    TAXONOMIC_SCREEN="CANDIDATE"
    echo "Norovirus taxonomic evidence found for ${sample}: ${KRAKEN_NORO_CLADE_FRAGMENTS} fragment(s), ${KRAKEN_NORO_DISTINCT_MINIMIZERS} distinct minimizer(s)."
    return 0
  fi

  TAXONOMIC_SCREEN="NO_EVIDENCE"
  echo "No Norovirus genus evidence found by Kraken2 for ${sample}."
  return 2
}

# write zero VP1 outputs when the taxonomic screen is negative
write_no_mapping_results() {
  local sample="$1"
  local vp1_type

  printf "%s\t%s\t%s\t%s\t0\t0.000000\t0.00\t0.00\t0.000000\n" \
    "${sample}" \
    "${RAW_INPUT_PAIRS}" \
    "${MAPPING_INPUT_PAIRS}" \
    "${RETAINED_PAIR_FRACTION}" \
    >> "${QC_FILE}"

  {
    printf "%s" "${sample}"
    for vp1_type in "${VP1_TYPES[@]}"
    do
      printf "\t0.000000"
    done
    printf "\n"
  } >> "${FILTERED_PROPORTIONS_FILE}"

  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "${sample}" \
    "${KRAKEN_NORO_PERCENT}" \
    "${KRAKEN_NORO_CLADE_FRAGMENTS}" \
    "${KRAKEN_NORO_DIRECT_FRAGMENTS}" \
    "${KRAKEN_NORO_TOTAL_MINIMIZERS}" \
    "${KRAKEN_NORO_DISTINCT_MINIMIZERS}" \
    "${TAXONOMIC_SCREEN}" \
    "NOT_SUPPORTED" \
    >> "${DETECTION_FILE}"
}

# map taxonomically screened samples against all VP1 references
map_sample() {
  local sample="$1"
  local bwa_log="${TMP_DIR}/${sample}.bwa.log"
  local sort_log="${TMP_DIR}/${sample}.samtools_sort.log"

  BAM="${TMP_DIR}/${sample}.multireference.sorted.bam"

  echo "Mapping ${sample} against the VP1 reference panel..."

  if ! bwa-mem2 mem \
    -t "${THREADS}" \
    -R "@RG\tID:${sample}\tSM:${sample}\tPL:ILLUMINA" \
    "${REFERENCE}" \
    "${MAPPING_R1}" \
    "${MAPPING_R2}" \
    2> "${bwa_log}" \
    | samtools sort \
        -@ "${SORT_THREADS}" \
        -o "${BAM}" \
        2> "${sort_log}"
  then
    echo "Skipping ${sample}: VP1 mapping failed." >&2
    cat "${bwa_log}" >&2 || true
    cat "${sort_log}" >&2 || true
    return 1
  fi

  if ! samtools index "${BAM}"
  then
    echo "Skipping ${sample}: BAM indexing failed." >&2
    return 1
  fi

  return 0
}

# calculate reference QC and final genotype proportions
calculate_mapping_metrics() {
  local sample="$1"
  local vp1_type
  local ref_name
  local reference_count_for_sample
  local reference_coverage
  local reference_mean_depth
  local reference_breadth_1x
  local reference_breadth_5x
  local best_breadth_5x
  local genotype_count
  local total_unique=0
  local count
  local proportion
  local mean_mapq
  local coverage_stats
  local mean_depth
  local breadth
  local assigned_fraction
  local i
  local passes_filter
  local filtered_total=0
  local filtered_count
  local filtered_proportion

  COUNTS=()
  PROPORTIONS=()
  GENOTYPE_BREADTH_5X=()
  FILTERED_COUNTS=()
  FILTERED_PROPORTIONS=()

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

    genotype_count=$(
      samtools view \
        -c \
        -q "${MIN_MAPQ}" \
        -f 66 \
        -F 2308 \
        "${BAM}" \
        "${GROUP_REFS[@]}"
    )

    COUNTS+=("${genotype_count}")
    total_unique=$((total_unique + genotype_count))
    best_breadth_5x="0.000000"

    for ref_name in "${GROUP_REFS[@]}"
    do
      reference_count_for_sample=$(
        samtools view \
          -c \
          -q "${MIN_MAPQ}" \
          -f 66 \
          -F 2308 \
          "${BAM}" \
          "${ref_name}"
      )

      reference_coverage=$(
        samtools depth \
          -aa \
          -q "${MIN_BASEQ}" \
          -Q "${MIN_MAPQ}" \
          -s \
          --require-flags 2 \
          --excl-flags 3844 \
          -r "${ref_name}" \
          "${BAM}" \
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
                printf "%.2f\t%.6f\t%.6f",
                  total_depth / positions,
                  covered_1x / positions,
                  covered_min_depth / positions
              } else {
                printf "0.00\t0.000000\t0.000000"
              }
            }'
      )

      reference_mean_depth=$(printf "%s\n" "${reference_coverage}" | cut -f1)
      reference_breadth_1x=$(printf "%s\n" "${reference_coverage}" | cut -f2)
      reference_breadth_5x=$(printf "%s\n" "${reference_coverage}" | cut -f3)

      printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "${sample}" \
        "${vp1_type}" \
        "${ref_name}" \
        "${reference_count_for_sample}" \
        "${reference_mean_depth}" \
        "${reference_breadth_1x}" \
        "${reference_breadth_5x}" \
        >> "${REFERENCE_QC_FILE}"

      if awk \
        -v current="${reference_breadth_5x}" \
        -v best="${best_breadth_5x}" \
        'BEGIN {exit !(current > best)}'
      then
        best_breadth_5x="${reference_breadth_5x}"
      fi
    done

    GENOTYPE_BREADTH_5X+=("${best_breadth_5x}")
  done

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
        -v total="${MAPPING_INPUT_PAIRS}" \
        'BEGIN {printf "%.6f", assigned / total}'
    )
  else
    for vp1_type in "${VP1_TYPES[@]}"
    do
      PROPORTIONS+=("0.000000")
    done
    assigned_fraction="0.000000"
  fi

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
        }'
  )

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
            printf "%.2f\t%.6f",
              total_depth / total_reference_length,
              total_covered_bases / total_reference_length
          } else {
            printf "0.00\t0.000000"
          }
        }'
  )

  mean_depth=$(printf "%s\n" "${coverage_stats}" | cut -f1)
  breadth=$(printf "%s\n" "${coverage_stats}" | cut -f2)

  for i in "${!VP1_TYPES[@]}"
  do
    passes_filter=$(
      awk \
        -v proportion="${PROPORTIONS[$i]}" \
        -v breadth="${GENOTYPE_BREADTH_5X[$i]}" \
        -v min_proportion="${MIN_PROPORTION}" \
        -v min_breadth="${MIN_BREADTH_5X}" \
        'BEGIN {
          if (proportion >= min_proportion && breadth >= min_breadth) {
            print 1
          } else {
            print 0
          }
        }'
    )

    if [[ "${passes_filter}" -eq 1 ]]
    then
      FILTERED_COUNTS+=("${COUNTS[$i]}")
      filtered_total=$((filtered_total + COUNTS[$i]))
    else
      FILTERED_COUNTS+=("0")
    fi
  done

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

  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "${sample}" \
    "${RAW_INPUT_PAIRS}" \
    "${MAPPING_INPUT_PAIRS}" \
    "${RETAINED_PAIR_FRACTION}" \
    "${total_unique}" \
    "${assigned_fraction}" \
    "${mean_mapq}" \
    "${mean_depth}" \
    "${breadth}" \
    >> "${QC_FILE}"

  {
    printf "%s" "${sample}"
    for filtered_proportion in "${FILTERED_PROPORTIONS[@]}"
    do
      printf "\t%s" "${filtered_proportion}"
    done
    printf "\n"
  } >> "${FILTERED_PROPORTIONS_FILE}"

  if [[ "${filtered_total}" -gt 0 ]]
  then
    FINAL_DETECTION="SUPPORTED_GII_VP1"
    echo "Norovirus GII VP1 signal supported for ${sample}."
  else
    FINAL_DETECTION="NOT_SUPPORTED"
    echo "Kraken2 found Norovirus evidence, but no GII VP1 genotype passed the mapping filters for ${sample}."
  fi

  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "${sample}" \
    "${KRAKEN_NORO_PERCENT}" \
    "${KRAKEN_NORO_CLADE_FRAGMENTS}" \
    "${KRAKEN_NORO_DIRECT_FRAGMENTS}" \
    "${KRAKEN_NORO_TOTAL_MINIMIZERS}" \
    "${KRAKEN_NORO_DISTINCT_MINIMIZERS}" \
    "${TAXONOMIC_SCREEN}" \
    "${FINAL_DETECTION}" \
    >> "${DETECTION_FILE}"
}

# read samples
read_samples

echo "Detected ${#SAMPLES[@]} samples in ${SEQUENCE_DIR}:"
for sample in "${SAMPLES[@]}"
do
  echo "  ${sample}"
done

# process samples
for sample in "${SAMPLES[@]}"
do
  echo
  echo "Processing ${sample}..."

  if ! find_sample_files "${sample}"
  then
    continue
  fi

  echo "R1: ${R1}"
  echo "R2: ${R2}"

  if ! prepare_reads "${sample}"
  then
    continue
  fi

  screen_norovirus "${sample}"
  screen_status=$?

  if [[ "${screen_status}" -eq 1 ]]
  then
    continue
  fi

  if [[ "${screen_status}" -eq 2 ]]
  then
    write_no_mapping_results "${sample}"
    continue
  fi

  if ! map_sample "${sample}"
  then
    rm -f "${BAM:-}" "${BAM:-}.bai" 2>/dev/null || true
    continue
  fi

  calculate_mapping_metrics "${sample}"

  rm -f "${BAM}" "${BAM}.bai"
done

# finish
rm -rf "${TMP_DIR}"
rm -f "${FASTQC_DIR}"/*_fastqc.zip 2>/dev/null || true

echo
echo "Analysis completed."
echo "Results: ${RESULTS_DIR}"
