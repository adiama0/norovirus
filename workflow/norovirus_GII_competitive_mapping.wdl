version 1.0

workflow NorovirusGIICompetitiveMapping {
  input {
    String samplename
    File read1
    File read2

    File vp1_reference
    File kraken2_db

    Float kraken2_confidence = 0.10
    Int norovirus_taxid = 142786

    Int min_mapq = 20
    Int min_baseq = 20
    Float min_proportion = 0.001
    Int min_depth = 5
    Float min_breadth_5x = 0.50

    Int trimmomatic_window = 4
    Int trimmomatic_window_quality = 20
    Int trimmomatic_minlen = 50

    Int fastqc_cpu = 4
    Int fastqc_memory_gb = 8
    Int fastqc_disk_gb = 30

    Int trim_cpu = 8
    Int trim_memory_gb = 16
    Int trim_disk_gb = 100

    Int kraken_cpu = 8
    Int kraken_memory_gb = 64
    Int kraken_disk_gb = 200

    Int mapping_cpu = 8
    Int mapping_sort_cpu = 4
    Int mapping_memory_gb = 32
    Int mapping_disk_gb = 200

    String fastqc_docker = "quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0"
    String trimmomatic_docker = "quay.io/biocontainers/trimmomatic:0.39--hdfd78af_2"
    String kraken2_docker = "quay.io/biocontainers/kraken2:2.17.1--pl5321h077b44d_0"
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

  call FastQC {
    input:
      samplename = samplename,
      read1 = read1,
      read2 = read2,
      cpu = fastqc_cpu,
      memory_gb = fastqc_memory_gb,
      disk_gb = fastqc_disk_gb,
      docker = fastqc_docker
  }

  call TrimReads {
    input:
      samplename = samplename,
      read1 = read1,
      read2 = read2,
      window = trimmomatic_window,
      window_quality = trimmomatic_window_quality,
      minlen = trimmomatic_minlen,
      cpu = trim_cpu,
      memory_gb = trim_memory_gb,
      disk_gb = trim_disk_gb,
      docker = trimmomatic_docker
  }

  call Kraken2Screen {
    input:
      samplename = samplename,
      read1 = TrimReads.trimmed_r1,
      read2 = TrimReads.trimmed_r2,
      kraken2_db = kraken2_db,
      confidence = kraken2_confidence,
      norovirus_taxid = norovirus_taxid,
      cpu = kraken_cpu,
      memory_gb = kraken_memory_gb,
      disk_gb = kraken_disk_gb,
      docker = kraken2_docker
  }

  if (Kraken2Screen.candidate) {
    call CompetitiveMapping {
      input:
        samplename = samplename,
        read1 = TrimReads.trimmed_r1,
        read2 = TrimReads.trimmed_r2,
        vp1_reference = vp1_reference,
        raw_read_pairs = TrimReads.raw_read_pairs,
        mapping_read_pairs = TrimReads.mapping_read_pairs,
        retained_pair_fraction = TrimReads.retained_pair_fraction,
        kraken_percent = Kraken2Screen.norovirus_percent,
        kraken_clade_fragments = Kraken2Screen.clade_fragments,
        kraken_direct_fragments = Kraken2Screen.direct_fragments,
        kraken_total_minimizers = Kraken2Screen.total_minimizers,
        kraken_distinct_minimizers = Kraken2Screen.distinct_minimizers,
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
  }

  if (!Kraken2Screen.candidate) {
    call NegativeMappingOutputs {
      input:
        samplename = samplename,
        vp1_reference = vp1_reference,
        raw_read_pairs = TrimReads.raw_read_pairs,
        mapping_read_pairs = TrimReads.mapping_read_pairs,
        retained_pair_fraction = TrimReads.retained_pair_fraction,
        kraken_percent = Kraken2Screen.norovirus_percent,
        kraken_clade_fragments = Kraken2Screen.clade_fragments,
        kraken_direct_fragments = Kraken2Screen.direct_fragments,
        kraken_total_minimizers = Kraken2Screen.total_minimizers,
        kraken_distinct_minimizers = Kraken2Screen.distinct_minimizers,
        docker = utility_docker
    }
  }

  output {
    File fastqc_r1_html = FastQC.r1_html
    File fastqc_r2_html = FastQC.r2_html

    File trimmed_r1 = TrimReads.trimmed_r1
    File trimmed_r2 = TrimReads.trimmed_r2
    File trimmomatic_log = TrimReads.log

    File kraken2_report = Kraken2Screen.report
    Float kraken2_norovirus_percent = Kraken2Screen.norovirus_percent
    Int kraken2_clade_fragments = Kraken2Screen.clade_fragments
    Int kraken2_distinct_minimizers = Kraken2Screen.distinct_minimizers
    String taxonomic_screen = Kraken2Screen.taxonomic_screen

    File norovirus_detection_tsv = select_first([
      CompetitiveMapping.norovirus_detection_tsv,
      NegativeMappingOutputs.norovirus_detection_tsv
    ])

    File mapping_qc_tsv = select_first([
      CompetitiveMapping.mapping_qc_tsv,
      NegativeMappingOutputs.mapping_qc_tsv
    ])

    File mapping_reference_qc_tsv = select_first([
      CompetitiveMapping.mapping_reference_qc_tsv,
      NegativeMappingOutputs.mapping_reference_qc_tsv
    ])

    File mapping_proportions_filtered_tsv = select_first([
      CompetitiveMapping.mapping_proportions_filtered_tsv,
      NegativeMappingOutputs.mapping_proportions_filtered_tsv
    ])

    String genotype_proportions = select_first([
      CompetitiveMapping.genotype_proportions,
      NegativeMappingOutputs.genotype_proportions
    ])

    String final_detection = select_first([
      CompetitiveMapping.final_detection,
      NegativeMappingOutputs.final_detection
    ])

    Int raw_read_pairs = TrimReads.raw_read_pairs
    Int mapping_read_pairs = TrimReads.mapping_read_pairs
    Float retained_pair_fraction = TrimReads.retained_pair_fraction

    Int vp1_fragments = select_first([
      CompetitiveMapping.vp1_fragments,
      NegativeMappingOutputs.vp1_fragments
    ])

    Float assigned_fraction = select_first([
      CompetitiveMapping.assigned_fraction,
      NegativeMappingOutputs.assigned_fraction
    ])

    Float mean_mapq = select_first([
      CompetitiveMapping.mean_mapq,
      NegativeMappingOutputs.mean_mapq
    ])

    Float mean_depth = select_first([
      CompetitiveMapping.mean_depth,
      NegativeMappingOutputs.mean_depth
    ])

    Float breadth = select_first([
      CompetitiveMapping.breadth,
      NegativeMappingOutputs.breadth
    ])

    meta {
    author: "Aron A. Diamond"
    email: "aron.diamond.a@gmail.com"
    description: "Sample-level workflow for detecting Norovirus in paired-end sequencing reads using Kraken2 and estimating Norovirus GII VP1 genotype proportions through competitive BWA-MEM2 mapping. The workflow performs FastQC, Trimmomatic read trimming, Kraken2 taxonomic screening, conditional competitive VP1 mapping, reference-level coverage assessment, genotype filtering, and reporting of genotype proportions and mapping quality metrics."
    }
  }
}


task FastQC {
  input {
    String samplename
    File read1
    File read2

    Int cpu
    Int memory_gb
    Int disk_gb
    String docker
  }

  command <<<
    set -euo pipefail

    ln -s "~{read1}" "~{samplename}_R1.fastq.gz"
    ln -s "~{read2}" "~{samplename}_R2.fastq.gz"

    fastqc \
      --quiet \
      --threads ~{cpu} \
      "~{samplename}_R1.fastq.gz" \
      "~{samplename}_R2.fastq.gz"
  >>>

  output {
    File r1_html = samplename + "_R1_fastqc.html"
    File r2_html = samplename + "_R2_fastqc.html"
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: memory_gb + " GB"
    disks: "local-disk " + disk_gb + " SSD"
  }
}


task TrimReads {
  input {
    String samplename
    File read1
    File read2

    Int window
    Int window_quality
    Int minlen

    Int cpu
    Int memory_gb
    Int disk_gb
    String docker
  }

  command <<<
    set -euo pipefail

    raw_r1=$(gzip -cd "~{read1}" | awk 'END {print NR / 4}')
    raw_r2=$(gzip -cd "~{read2}" | awk 'END {print NR / 4}')

    if [[ "${raw_r1}" -le 0 || "${raw_r2}" -le 0 ]]; then
      echo "ERROR: empty or unreadable FASTQ input for ~{samplename}" >&2
      exit 1
    fi

    if [[ "${raw_r1}" -ne "${raw_r2}" ]]; then
      echo "ERROR: R1 and R2 contain different read counts for ~{samplename}" >&2
      exit 1
    fi

    trimmomatic PE \
      -threads ~{cpu} \
      -phred33 \
      "~{read1}" \
      "~{read2}" \
      "~{samplename}_R1.trimmed.paired.fastq.gz" \
      "~{samplename}_R1.trimmed.unpaired.fastq.gz" \
      "~{samplename}_R2.trimmed.paired.fastq.gz" \
      "~{samplename}_R2.trimmed.unpaired.fastq.gz" \
      "SLIDINGWINDOW:~{window}:~{window_quality}" \
      "MINLEN:~{minlen}" \
      > "~{samplename}.trimmomatic.log" 2>&1

    mapping_r1=$(gzip -cd "~{samplename}_R1.trimmed.paired.fastq.gz" | awk 'END {print NR / 4}')
    mapping_r2=$(gzip -cd "~{samplename}_R2.trimmed.paired.fastq.gz" | awk 'END {print NR / 4}')

    if [[ "${mapping_r1}" -le 0 || "${mapping_r2}" -le 0 ]]; then
      echo "ERROR: no paired reads remain after trimming for ~{samplename}" >&2
      exit 1
    fi

    if [[ "${mapping_r1}" -ne "${mapping_r2}" ]]; then
      echo "ERROR: trimmed R1 and R2 contain different read counts for ~{samplename}" >&2
      exit 1
    fi

    awk -v raw="${raw_r1}" -v kept="${mapping_r1}" \
      'BEGIN {printf "%.6f\n", kept / raw}' \
      > retained_pair_fraction.txt

    printf "%s\n" "${raw_r1}" > raw_read_pairs.txt
    printf "%s\n" "${mapping_r1}" > mapping_read_pairs.txt

    rm -f \
      "~{samplename}_R1.trimmed.unpaired.fastq.gz" \
      "~{samplename}_R2.trimmed.unpaired.fastq.gz"
  >>>

  output {
    File trimmed_r1 = samplename + "_R1.trimmed.paired.fastq.gz"
    File trimmed_r2 = samplename + "_R2.trimmed.paired.fastq.gz"
    File log = samplename + ".trimmomatic.log"

    Int raw_read_pairs = read_int("raw_read_pairs.txt")
    Int mapping_read_pairs = read_int("mapping_read_pairs.txt")
    Float retained_pair_fraction = read_float("retained_pair_fraction.txt")
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: memory_gb + " GB"
    disks: "local-disk " + disk_gb + " SSD"
  }
}


task Kraken2Screen {
  input {
    String samplename
    File read1
    File read2
    File kraken2_db

    Float confidence
    Int norovirus_taxid

    Int cpu
    Int memory_gb
    Int disk_gb
    String docker
  }

  command <<<
    set -euo pipefail

    mkdir -p kraken_db
    tar -xzf "~{kraken2_db}" -C kraken_db

    hash_file=$(find kraken_db -type f -name hash.k2d -print -quit)

    if [[ -z "${hash_file}" ]]; then
      echo "ERROR: hash.k2d was not found after extracting the Kraken2 database." >&2
      exit 1
    fi

    db_dir=$(dirname "${hash_file}")

    kraken2 \
      --db "${db_dir}" \
      --threads ~{cpu} \
      --paired \
      --confidence ~{confidence} \
      --report "~{samplename}.kraken2.report.tsv" \
      --report-minimizer-data \
      --memory-mapping \
      --output /dev/null \
      "~{read1}" \
      "~{read2}"

    metrics=$(
      awk -F '\t' -v taxid="~{norovirus_taxid}" '
        $7 == taxid {
          printf "%s\t%s\t%s\t%s\t%s", $1, $2, $3, $4, $5
          found = 1
        }
        END {
          if (!found) {
            printf "0.00\t0\t0\t0\t0"
          }
        }
      ' "~{samplename}.kraken2.report.tsv"
    )

    percent=$(printf "%s\n" "${metrics}" | cut -f1)
    clade=$(printf "%s\n" "${metrics}" | cut -f2)
    direct=$(printf "%s\n" "${metrics}" | cut -f3)
    total_minimizers=$(printf "%s\n" "${metrics}" | cut -f4)
    distinct_minimizers=$(printf "%s\n" "${metrics}" | cut -f5)

    printf "%s\n" "${percent}" > norovirus_percent.txt
    printf "%s\n" "${clade}" > clade_fragments.txt
    printf "%s\n" "${direct}" > direct_fragments.txt
    printf "%s\n" "${total_minimizers}" > total_minimizers.txt
    printf "%s\n" "${distinct_minimizers}" > distinct_minimizers.txt

    if [[ "${clade}" -gt 0 && "${distinct_minimizers}" -gt 0 ]]; then
      printf "true\n" > candidate.txt
      printf "CANDIDATE\n" > taxonomic_screen.txt
    else
      printf "false\n" > candidate.txt
      printf "NO_EVIDENCE\n" > taxonomic_screen.txt
    fi
  >>>

  output {
    File report = samplename + ".kraken2.report.tsv"

    Float norovirus_percent = read_float("norovirus_percent.txt")
    Int clade_fragments = read_int("clade_fragments.txt")
    Int direct_fragments = read_int("direct_fragments.txt")
    Int total_minimizers = read_int("total_minimizers.txt")
    Int distinct_minimizers = read_int("distinct_minimizers.txt")
    Boolean candidate = read_boolean("candidate.txt")
    String taxonomic_screen = read_string("taxonomic_screen.txt")
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: memory_gb + " GB"
    disks: "local-disk " + disk_gb + " SSD"
  }
}


task CompetitiveMapping {
  input {
    String samplename
    File read1
    File read2
    File vp1_reference

    Int raw_read_pairs
    Int mapping_read_pairs
    Float retained_pair_fraction

    Float kraken_percent
    Int kraken_clade_fragments
    Int kraken_direct_fragments
    Int kraken_total_minimizers
    Int kraken_distinct_minimizers

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

    bwa-mem2 mem \
      -t ~{cpu} \
      -R "@RG\tID:~{samplename}\tSM:~{samplename}\tPL:ILLUMINA" \
      reference_multi.fasta \
      "~{read1}" \
      "~{read2}" \
      2> "~{samplename}.bwa.log" \
      | samtools sort \
          -@ ~{sort_cpu} \
          -o "~{samplename}.multireference.sorted.bam" \
          2> "~{samplename}.samtools_sort.log"

    BAM="~{samplename}.multireference.sorted.bam"
    samtools index "${BAM}"

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
          -f 66 \
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
            -f 66 \
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

      assigned_fraction=$(awk -v assigned="${total_unique}" -v total="~{mapping_read_pairs}" \
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

    printf "sample\traw_read_pairs\tmapping_read_pairs\tretained_pair_fraction\tunique_VP1_fragments\tassigned_fraction\tmean_MAPQ\tmean_depth\tbreadth\n" \
      > mapping_qc.tsv

    printf "%s\t%s\t%s\t%.6f\t%s\t%s\t%s\t%s\t%s\n" \
      "~{samplename}" \
      "~{raw_read_pairs}" \
      "~{mapping_read_pairs}" \
      "~{retained_pair_fraction}" \
      "${total_unique}" \
      "${assigned_fraction}" \
      "${mean_mapq}" \
      "${mean_depth}" \
      "${breadth}" \
      >> mapping_qc.tsv

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

    printf "sample\tkraken2_norovirus_percent\tkraken2_clade_fragments\tkraken2_direct_fragments\tkraken2_total_minimizers\tkraken2_distinct_minimizers\ttaxonomic_screen\tfinal_detection\n" \
      > norovirus_detection.tsv

    printf "%s\t%.2f\t%s\t%s\t%s\t%s\tCANDIDATE\t%s\n" \
      "~{samplename}" \
      "~{kraken_percent}" \
      "~{kraken_clade_fragments}" \
      "~{kraken_direct_fragments}" \
      "~{kraken_total_minimizers}" \
      "~{kraken_distinct_minimizers}" \
      "${final_detection}" \
      >> norovirus_detection.tsv

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

    rm -f "${BAM}" "${BAM}.bai"
  >>>

  output {
    File norovirus_detection_tsv = "norovirus_detection.tsv"
    File mapping_qc_tsv = "mapping_qc.tsv"
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


task NegativeMappingOutputs {
  input {
    String samplename
    File vp1_reference

    Int raw_read_pairs
    Int mapping_read_pairs
    Float retained_pair_fraction

    Float kraken_percent
    Int kraken_clade_fragments
    Int kraken_direct_fragments
    Int kraken_total_minimizers
    Int kraken_distinct_minimizers

    String docker
  }

  command <<<
    set -euo pipefail

    mapfile -t VP1_TYPES < <(
      grep '^>' "~{vp1_reference}" \
        | sed 's/^>//' \
        | awk '{print $1}' \
        | cut -d'|' -f2 \
        | awk '!seen[$0]++'
    )

    printf "sample\traw_read_pairs\tmapping_read_pairs\tretained_pair_fraction\tunique_VP1_fragments\tassigned_fraction\tmean_MAPQ\tmean_depth\tbreadth\n" \
      > mapping_qc.tsv

    printf "%s\t%s\t%s\t%.6f\t0\t0.000000\t0.00\t0.00\t0.000000\n" \
      "~{samplename}" \
      "~{raw_read_pairs}" \
      "~{mapping_read_pairs}" \
      "~{retained_pair_fraction}" \
      >> mapping_qc.tsv

    {
      printf "sample"
      for vp1_type in "${VP1_TYPES[@]}"; do
        printf "\t%s" "${vp1_type}"
      done
      printf "\n"

      printf "%s" "~{samplename}"
      for vp1_type in "${VP1_TYPES[@]}"; do
        printf "\t0.000000"
      done
      printf "\n"
    } > mapping_proportions_filtered.tsv

    printf "sample\tVP1_type\treference\tfragments\tmean_depth\tbreadth_1x\tbreadth_5x\n" \
      > mapping_reference_qc.tsv

    printf "sample\tkraken2_norovirus_percent\tkraken2_clade_fragments\tkraken2_direct_fragments\tkraken2_total_minimizers\tkraken2_distinct_minimizers\ttaxonomic_screen\tfinal_detection\n" \
      > norovirus_detection.tsv

    printf "%s\t%.2f\t%s\t%s\t%s\t%s\tNO_EVIDENCE\tNOT_SUPPORTED\n" \
      "~{samplename}" \
      "~{kraken_percent}" \
      "~{kraken_clade_fragments}" \
      "~{kraken_direct_fragments}" \
      "~{kraken_total_minimizers}" \
      "~{kraken_distinct_minimizers}" \
      >> norovirus_detection.tsv

    printf "NONE\n" > genotype_proportions.txt
    printf "NOT_SUPPORTED\n" > final_detection.txt

    printf "0\n" > vp1_fragments.txt
    printf "0.000000\n" > assigned_fraction.txt
    printf "0.00\n" > mean_mapq.txt
    printf "0.00\n" > mean_depth.txt
    printf "0.000000\n" > breadth.txt
  >>>

  output {
    File norovirus_detection_tsv = "norovirus_detection.tsv"
    File mapping_qc_tsv = "mapping_qc.tsv"
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
    cpu: 1
    memory: "2 GB"
    disks: "local-disk 10 SSD"
  }
}