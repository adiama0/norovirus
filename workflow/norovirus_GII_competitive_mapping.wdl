version 1.0

workflow NorovirusGIICompetitiveMapping {
  input {
    String samplename
    File read1
    File read2

    File vp1_reference

    Int min_mapq = 20
    Int min_baseq = 20
    Float min_proportion = 0.001
    Int min_depth = 5
    Float min_breadth_5x = 0.50

    Int mapping_cpu = 8
    Int mapping_sort_cpu = 4
    Int mapping_memory_gb = 32
    Int mapping_disk_gb = 200

    String mapping_docker = "hmartiniano/bwa-mem2-samtools@sha256:f549b194c53c85ba4cced3555132e2f3ae79d5e40c74e8b576c2bc09fcfe9648"
    String utility_docker = "ubuntu:24.04"
  }

  parameter_meta {
    samplename: "Sample identifier. In Terra, use a sample-table attribute such as this.sample_name."
    read1: "Paired-end R1 FASTQ.GZ file."
    read2: "Paired-end R2 FASTQ.GZ file."
    vp1_reference: "VP1 multi-reference FASTA. In Terra, a literal GCS path can be entered as \"gs://bucket/reference_multi.fasta\"."
    kraken2_db: "Kraken2 database packaged as one .tar.gz file. In Terra, enter a literal GCS path such as \"gs://bucket/kraken2_database.tar.gz\"."
  }

  call CompetitiveMapping {
    input:
      samplename = samplename,
      read1 = read1,
      read2 = read2,
      vp1_reference = vp1_reference,
      min_mapq = min_mapq,
      min_baseq = min_baseq,
      min_proportion = min_proportion,
      min_depth = min_depth,
      min_breadth_5x = min_breadth_5x,
      cpu = mapping_cpu,
      sort_cpu = mapping_sort_cpu,
      memory_gb = mapping_memory_gb,
      disk_gb = mapping_disk_gb,
      docker = mapping_docker
    }


  output {


    File mapping_reference_qc_tsv = CompetitiveMapping.mapping_reference_qc_tsv

    File mapping_proportions_filtered_tsv = CompetitiveMapping.mapping_proportions_filtered_tsv

    File? aligned_bam = CompetitiveMapping.aligned_bam
    File? aligned_bam_index = CompetitiveMapping.aligned_bam_index

    String genotype_proportions = CompetitiveMapping.genotype_proportions

    String final_detection = CompetitiveMapping.final_detection

    Int vp1_fragments = CompetitiveMapping.vp1_fragments

    Float assigned_fraction = CompetitiveMapping.assigned_fraction

    Float mean_mapq = CompetitiveMapping.mean_mapq

    Float mean_depth = CompetitiveMapping.mean_depth

    Float breadth = CompetitiveMapping.breadth
  }
    
meta {
    author: "Aron A. Diamond"
    email: "aron.diamond.a@gmail.com"
    description: "Sample-level workflow for detecting Norovirus in paired-end sequencing reads using Kraken2 and estimating Norovirus GII VP1 genotype proportions through competitive BWA-MEM2 mapping. The workflow performs FastQC, Trimmomatic read trimming, Kraken2 taxonomic screening, conditional competitive VP1 mapping, reference-level coverage assessment, genotype filtering, and reporting of genotype proportions and mapping quality metrics."
    }
}


task CompetitiveMapping {
  input {
    String samplename
    File read1
    File read2
    File vp1_reference

    Int min_mapq
    Int min_baseq
    Float min_proportion
    Int min_depth
    Float min_breadth_5x

    Int cpu
    Int sort_cpu
    Int memory_gb
    Int disk_gb
    String docker
  }

  command <<<
    set -euo pipefail

    cp "~{vp1_reference}" reference_multi.fasta

    bwa-mem2 index reference_multi.fasta
    samtools faidx reference_multi.fasta

    REF_NAMES=()
    REF_TYPES=()
    VP1_TYPES=()

    while IFS= read -r ref_name; do
      vp1_type=$(printf "%s\n" "${ref_name}" | cut -d'|' -f2)

      if [[ -z "${vp1_type}" || "${vp1_type}" == "${ref_name}" ]]; then
        echo "ERROR: invalid reference header: ${ref_name}" >&2
        exit 1
      fi

      REF_NAMES+=("${ref_name}")
      REF_TYPES+=("${vp1_type}")
    done < <(grep '^>' reference_multi.fasta | sed 's/^>//' | awk '{print $1}')

    while IFS= read -r vp1_type; do
      VP1_TYPES+=("${vp1_type}")
    done < <(
      grep '^>' reference_multi.fasta \
        | sed 's/^>//' \
        | awk '{print $1}' \
        | cut -d'|' -f2 \
        | awk '!seen[$0]++'
    )

    if [[ "${#REF_NAMES[@]}" -eq 0 ]]; then
      echo "ERROR: no VP1 references were found." >&2
      exit 1
    fi
    
    mkdir -p bams

    bwa-mem2 mem \
      -t ~{cpu} \
      -R "@RG\tID:~{samplename}\tSM:~{samplename}\tPL:ILLUMINA" \
      reference_multi.fasta \
      "~{read1}" \
      "~{read2}" \
      2> "~{samplename}.bwa.log" \
      | samtools sort \
          -@ ~{sort_cpu} \
          -o "bams/~{samplename}.multireference.sorted.bam" \
          2> "~{samplename}.samtools_sort.log"

    BAM="bams/~{samplename}.multireference.sorted.bam"
    samtools index "${BAM}"
      mapping_input_pairs=$(samtools view -c -f 64 -F 2304 "${BAM}")

    printf "sample\tVP1_type\treference\tfragments\tmean_depth\tbreadth_1x\tbreadth_5x\n" \
      > mapping_reference_qc.tsv

    COUNTS=()
    PROPORTIONS=()
    GENOTYPE_BREADTH_5X=()
    FILTERED_COUNTS=()
    FILTERED_PROPORTIONS=()

    total_unique=0

    for vp1_type in "${VP1_TYPES[@]}"; do
      GROUP_REFS=()

      for i in "${!REF_NAMES[@]}"; do
        if [[ "${REF_TYPES[$i]}" == "${vp1_type}" ]]; then
          GROUP_REFS+=("${REF_NAMES[$i]}")
        fi
      done

      genotype_count=$(
        samtools view \
          -c \
          -q ~{min_mapq} \
          -f 64 \
          -F 2308 \
          "${BAM}" \
          "${GROUP_REFS[@]}"
      )

      COUNTS+=("${genotype_count}")
      total_unique=$((total_unique + genotype_count))
      best_breadth_5x="0.000000"

      for ref_name in "${GROUP_REFS[@]}"; do
        reference_count_for_sample=$(
          samtools view \
            -c \
            -q ~{min_mapq} \
            -f 64 \
            -F 2308 \
            "${BAM}" \
            "${ref_name}"
        )

        reference_coverage=$(
          samtools depth \
            -aa \
            -q ~{min_baseq} \
            -Q ~{min_mapq} \
            -s \
            --require-flags 2 \
            --excl-flags 3844 \
            -r "${ref_name}" \
            "${BAM}" \
            | awk -v min_depth="~{min_depth}" '
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
          "~{samplename}" \
          "${vp1_type}" \
          "${ref_name}" \
          "${reference_count_for_sample}" \
          "${reference_mean_depth}" \
          "${reference_breadth_1x}" \
          "${reference_breadth_5x}" \
          >> mapping_reference_qc.tsv

        if awk -v current="${reference_breadth_5x}" -v best="${best_breadth_5x}" \
          'BEGIN {exit !(current > best)}'; then
          best_breadth_5x="${reference_breadth_5x}"
        fi
      done

      GENOTYPE_BREADTH_5X+=("${best_breadth_5x}")
    done

    if [[ "${total_unique}" -gt 0 ]]; then
      for count in "${COUNTS[@]}"; do
        proportion=$(awk -v count="${count}" -v total="${total_unique}" \
          'BEGIN {printf "%.6f", count / total}')
        PROPORTIONS+=("${proportion}")
      done

      assigned_fraction=$(awk -v assigned="${total_unique}" -v total="${mapping_input_pairs}" \
        'BEGIN {printf "%.6f", assigned / total}')
    else
      for vp1_type in "${VP1_TYPES[@]}"; do
        PROPORTIONS+=("0.000000")
      done
      assigned_fraction="0.000000"
    fi

    mean_mapq=$(
      samtools view -f 66 -F 2308 "${BAM}" \
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
        -q ~{min_mapq} \
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

    filtered_total=0

    for i in "${!VP1_TYPES[@]}"; do
      passes_filter=$(
        awk \
          -v proportion="${PROPORTIONS[$i]}" \
          -v breadth="${GENOTYPE_BREADTH_5X[$i]}" \
          -v min_proportion="~{min_proportion}" \
          -v min_breadth="~{min_breadth_5x}" \
          'BEGIN {
            if (proportion >= min_proportion && breadth >= min_breadth) {
              print 1
            } else {
              print 0
            }
          }'
      )

      if [[ "${passes_filter}" -eq 1 ]]; then
        FILTERED_COUNTS+=("${COUNTS[$i]}")
        filtered_total=$((filtered_total + COUNTS[$i]))
      else
        FILTERED_COUNTS+=("0")
      fi
    done

    if [[ "${filtered_total}" -gt 0 ]]; then
      for filtered_count in "${FILTERED_COUNTS[@]}"; do
        filtered_proportion=$(awk -v count="${filtered_count}" -v total="${filtered_total}" \
          'BEGIN {printf "%.6f", count / total}')
        FILTERED_PROPORTIONS+=("${filtered_proportion}")
      done
      final_detection="SUPPORTED_GII_VP1"
    else
      for vp1_type in "${VP1_TYPES[@]}"; do
        FILTERED_PROPORTIONS+=("0.000000")
      done
      final_detection="NOT_SUPPORTED"
    fi

    {
      printf "sample"
      for vp1_type in "${VP1_TYPES[@]}"; do
        printf "\t%s" "${vp1_type}"
      done
      printf "\n"

      printf "%s" "~{samplename}"
      for filtered_proportion in "${FILTERED_PROPORTIONS[@]}"; do
        printf "\t%s" "${filtered_proportion}"
      done
      printf "\n"
    } > mapping_proportions_filtered.tsv

    genotype_text=""

    for i in "${!VP1_TYPES[@]}"; do
      p="${FILTERED_PROPORTIONS[$i]}"

      if awk -v p="${p}" 'BEGIN {exit !(p > 0)}'; then
        if [[ -n "${genotype_text}" ]]; then
          genotype_text="${genotype_text};"
        fi
        genotype_text="${genotype_text}${VP1_TYPES[$i]}=${p}"
      fi
    done

    if [[ -z "${genotype_text}" ]]; then
      genotype_text="NONE"
    fi

    printf "%s\n" "${genotype_text}" > genotype_proportions.txt
    printf "%s\n" "${final_detection}" > final_detection.txt

    printf "%s\n" "${total_unique}" > vp1_fragments.txt
    printf "%s\n" "${assigned_fraction}" > assigned_fraction.txt
    printf "%s\n" "${mean_mapq}" > mean_mapq.txt
    printf "%s\n" "${mean_depth}" > mean_depth.txt
    printf "%s\n" "${breadth}" > breadth.txt

  >>>

  output {
    File aligned_bam = "bams/" + samplename + ".multireference.sorted.bam"
    File aligned_bam_index = "bams/" + samplename + ".multireference.sorted.bam.bai"

    File mapping_reference_qc_tsv = "mapping_reference_qc.tsv"
    File mapping_proportions_filtered_tsv = "mapping_proportions_filtered.tsv"

    String genotype_proportions = read_string("genotype_proportions.txt")
    String final_detection = read_string("final_detection.txt")

    Int vp1_fragments = read_int("vp1_fragments.txt")
    Float assigned_fraction = read_float("assigned_fraction.txt")
    Float mean_mapq = read_float("mean_mapq.txt")
    Float mean_depth = read_float("mean_depth.txt")
    Float breadth = read_float("breadth.txt")
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: memory_gb + " GB"
    disks: "local-disk " + disk_gb + " SSD"
  }
}
