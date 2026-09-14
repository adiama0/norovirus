# Norovirus GII VP1 Competitive Mapping Workflow

## Overview

This repository contains a workflow for detecting norovirus in mixed metagenomic sequencing samples and estimating the relative abundance of norovirus GII VP1 genotypes.

The workflow separates two questions:

1. **Is there taxonomic evidence for Norovirus in the sample?**
2. **If Norovirus is supported, which GII VP1 genotypes are present and in what relative proportions?**

The current pipeline uses Kraken2 as a broad metagenomic screen before competitive mapping against a curated multi-reference VP1 panel. The VP1 panel can contain more than one representative sequence for the same genotype when one medoid does not adequately represent the diversity observed during validation.

The current analysis is:

```text
paired input reads
        |
        v
validate R1/R2 and count reads
        |
        +---------------------------+
        |                           |
      FASTQ                       FASTA
        |                           |
        v                           |
      FastQC                        |
        |                           |
        v                           |
   Trimmomatic                      |
        |                           |
        +-------------+-------------+
                      |
                      v
             Kraken2 taxonomic screen
                      |
              Norovirus evidence?
                 /          \
               no            yes
               |              |
               v              v
      report no supported   BWA-MEM2 competitive
      GII VP1 signal        mapping to VP1 panel
                                  |
                                  v
                         reference-level QC
                                  |
                                  v
                         collapse references
                           by VP1 genotype
                                  |
                                  v
                       abundance + breadth
                              filtering
                                  |
                                  v
                        filtered genotype
                           proportions
```

---

# Goals of the pipeline

Wastewater and other metagenomic samples contain reads from many organisms. Norovirus may represent only a small fraction of the total library, and more than one norovirus genotype can be present at the same time.

The pipeline was designed to address several problems:

- distinguish taxonomic evidence for Norovirus from genotype assignment;
- reduce reference bias caused by representing an entire genotype with one sequence;
- allow multiple GII genotypes to be quantified in the same sample;
- avoid counting the same paired fragment twice;
- use depth and breadth to distinguish broad VP1 support from localized cross-mapping;
- report genotype proportions only after unsupported assignments are filtered.

The final proportions describe the relative composition of the **supported GII VP1 mixture**. They are not the fraction of the complete wastewater metagenome and are not an absolute viral-load measurement.

---

# Directory structure

```text
norovirus/
|
├── data/
│   ├── references/
│   │   └── reference_multi.fasta
│   |
│   ├── kraken2/
│   │   └── Kraken2 database files
│   |
│   ├── sequences/
│   │   └── clinical/
│   │       ├── SAMPLE_R1.fastq.gz
│   │       ├── SAMPLE_R2.fastq.gz
│   │       └── ...
│   |
│   └── results/
│       ├── fastqc/
│       │   └── *_fastqc.html
│       |
│       ├── trimmomatic/
│       │   ├── *_R1.trimmed.paired.fastq.gz
│       │   ├── *_R2.trimmed.paired.fastq.gz
│       │   └── *.trimmomatic.log
│       |
│       ├── taxonomy/
│       │   ├── *.kraken2.report.tsv
│       │   └── norovirus_detection.tsv
│       |
│       └── competitive_mapping/
│           ├── mapping_qc.tsv
│           ├── mapping_reference_qc.tsv
│           └── mapping_proportions_filtered.tsv
│
├── scripts/
│   └── 
│
├── env.yml
├── run_competitive_mapping.sh
└── README.md
```

`run_competitive_mapping.sh` removes and recreates `data/results/` at the beginning of a complete run. Copy previous results elsewhere before rerunning if they need to be retained.

---

# Software environment

Create the Conda environment with:

```bash
conda env create -f env.yml
conda activate NV
```

The current competitive-mapping workflow directly uses:

- BWA-MEM2
- samtools
- seqkit
- FastQC
- Trimmomatic
- Kraken2

Additional packages in the project environment are used for reference construction, mock generation, and downstream analysis.

Kraken2 requires both the Kraken2 software and a Kraken2 database. The default database path in the shell workflow is:

```text
data/kraken2
```

A different database can be supplied with:

```bash
KRAKEN2_DB=/path/to/kraken2_database bash run_competitive_mapping.sh
```

# Why use multi-reference competitive mapping?

A genotype is not represented by one invariant nucleotide sequence. Strains belonging to the same genotype can differ enough that a single reference recruits some strains much better than others.

Using only one representative can therefore cause reference bias:

```text
true strain close to reference
    -> efficient mapping

true strain distant from reference
    -> fewer reads map
    -> lower breadth
    -> increased cross-mapping risk
```

A multi-reference FASTA provides several plausible representatives for genotypes that require more sequence diversity.

Competitive mapping means that all VP1 references are indexed together and every read is evaluated against the entire panel in the same BWA-MEM2 alignment.

This is preferable to mapping independently to each genotype because independent mapping could allow the same fragment to be counted against more than one genotype.

After alignment, the workflow retains qualifying primary alignments and collapses all representatives belonging to the same genotype.

Example:

```text
GII.17 rep1 = 1000 fragments
GII.17 rep2 =  400 fragments
GII.17 rep3 =  100 fragments

GII.17 total = 1500 fragments
```

---

# Input data

The current workflow uses a strict paired-end naming convention.

FASTQ input:

```text
SAMPLE_R1.fastq.gz
SAMPLE_R2.fastq.gz
```

FASTA input:

```text
SAMPLE_R1.fasta
SAMPLE_R2.fasta
```

Inputs are read from:

```text
data/sequences/clinical/
```

R1 and R2 must contain the same number of sequences.

A sample cannot have both FASTQ and FASTA input files with the same sample name.

---

# Running the workflow

Run:

```bash
bash run_competitive_mapping.sh
```

The main reference is:

```text
data/references/reference_multi.fasta
```

The script indexes the reference with BWA-MEM2 and samtools before processing samples.

---

# Raw-read QC

## FASTQ input

Raw FASTQ pairs are analyzed with FastQC before trimming.

FastQC is used to inspect sequencing properties such as:

- per-base quality;
- sequence length;
- GC distribution;
- adapter content;
- duplicated or overrepresented sequences.

FastQC is a technical QC step. It does not determine whether Norovirus is present.

The current streamlined workflow runs FastQC on the raw reads only and does not automatically rerun FastQC after Trimmomatic.

If FastQC cannot process a sample, that sample is skipped while the workflow continues with the remaining samples.

## FASTA input

FASTA files do not contain Phred quality scores.

FASTA samples therefore skip FastQC and Trimmomatic and proceed directly to taxonomic screening.

---

# Trimmomatic preprocessing

FASTQ reads are trimmed before Kraken2 or BWA-MEM2 analysis.

Current quality-trimming settings are:

```text
SLIDINGWINDOW:4:20
MINLEN:50
```

`SLIDINGWINDOW:4:20` trims when the average quality in a four-base window drops below Q20.

`MINLEN:50` removes reads shorter than 50 nucleotides after trimming.

Only paired survivors are retained for downstream analysis. Unpaired Trimmomatic outputs are temporary and are deleted.

Adapter trimming is optional. If an adapter FASTA is known, it can be supplied with:

```bash
TRIMMOMATIC_ADAPTERS=/path/to/adapters.fa \
  bash run_competitive_mapping.sh
```

When supplied, the workflow adds:

```text
ILLUMINACLIP:<adapter file>:2:30:10
```

---

# Kraken2 Norovirus screen

After preprocessing, Kraken2 is used as a broad taxonomic screening step.

The purpose is to separate:

```text
Is there taxonomic evidence for Norovirus?
```

from:

```text
Which GII VP1 genotype best explains the reads?
```

The workflow currently uses:

```text
Norovirus taxid: 142786
Kraken2 confidence: 0.10
paired-end mode
report minimizer data
memory mapping
```

Kraken2 reports are written to:

```text
data/results/taxonomy/SAMPLE.kraken2.report.tsv
```

The workflow records:

- percentage assigned to the Norovirus clade;
- Norovirus clade fragments;
- fragments assigned directly to the Norovirus genus;
- total minimizers;
- distinct minimizers.

A sample is marked as a Kraken2 `CANDIDATE` when the Norovirus genus has both:

```text
clade fragments > 0
```

and:

```text
distinct minimizers > 0
```

These are current screening rules and are not a validated diagnostic limit of detection.

If no Norovirus evidence is found, VP1 competitive mapping is skipped and zero genotype proportions are written for that sample.

---

# Important current Kraken2-to-BWA behavior

Kraken2 currently acts as a **gate**, not a read extractor.

If Kraken2 detects Norovirus evidence, the current shell script maps the complete set of processed paired reads against the VP1 panel.

It does not currently extract only Kraken2-classified Norovirus or Caliciviridae reads before BWA-MEM2.

Current behavior:

```text
all processed reads
        |
        v
Kraken2 screen
        |
        +-- no Norovirus evidence --> no VP1 mapping
        |
        +-- Norovirus evidence --> all processed reads go to VP1 mapping
```

This implementation detail should be considered when interpreting very low-level competitive-mapping assignments.

---

# Competitive VP1 mapping

Samples that pass the Kraken2 screen are aligned against:

```text
data/references/reference_multi.fasta
```

with BWA-MEM2.

All reference sequences are present in the same BWA index, so they compete for each read during a single alignment.

The BAM file is sorted and indexed with samtools, used for QC and genotype calculations, and then deleted to reduce disk usage.

---

# Fragment counting and alignment filters

The workflow counts one R1 record from each qualifying paired fragment.

The main requirements are:

```text
properly paired
R1 / first read in pair
MAPQ >= 20
not unmapped
not secondary
not supplementary
```

The count command uses the equivalent of:

```text
-f 66
-q 20
-F 2308
```

`-f 66` requires:

```text
64 = first read in pair
2  = properly paired
```

Counting only R1 means one physical paired fragment contributes one count rather than two.

`-F 2308` excludes:

```text
4    unmapped
256  secondary alignment
2048 supplementary alignment
```

This prevents secondary or supplementary records from being counted as additional fragments.

---

# Reference-level coverage QC

Coverage is calculated independently for every reference in the multi-reference panel.

The workflow reports:

- fragment count;
- mean depth;
- breadth at >=1x;
- breadth at >=5x.

`breadth_1x` is the fraction of reference positions covered by at least one read.

`breadth_5x` is the fraction covered at depth >=5.

References belonging to the same genotype are **not** concatenated for breadth calculation. Concatenating several alternative representatives would artificially penalize a genotype simply because more representatives were added.

For genotype filtering, the highest `breadth_5x` observed among that genotype's representatives is used.

---

# Genotype proportions and filtering

After reference-level counts are obtained, counts from all representatives belonging to the same genotype are summed.

Raw within-GII genotype proportions are calculated as:

```text
genotype qualifying fragments
--------------------------------
total qualifying GII VP1 fragments
```

The current filter requires:

```text
raw genotype proportion >= 0.001
AND
best representative breadth_5x >= 0.50
```

Equivalent thresholds:

```text
minimum raw proportion = 0.1%
minimum breadth at >=5x = 50%
```

The abundance filter reduces very small competitive-mapping assignments.

The breadth filter helps reject cases where many reads pile up over only a small conserved region.

After filtering, passing genotype counts are renormalized so that the supported genotype proportions sum to approximately 1.0.

These thresholds are method-development parameters and should not be interpreted as a validated clinical detection threshold.

---

# Output files

## FastQC

```text
data/results/fastqc/
```

The workflow retains FastQC HTML reports and removes the FastQC ZIP files at the end of the run.

## Trimmomatic

```text
data/results/trimmomatic/
```

Contains:

```text
SAMPLE_R1.trimmed.paired.fastq.gz
SAMPLE_R2.trimmed.paired.fastq.gz
SAMPLE.trimmomatic.log
```

## `norovirus_detection.tsv`

```text
data/results/taxonomy/norovirus_detection.tsv
```

Columns:

```text
sample
kraken2_norovirus_percent
kraken2_clade_fragments
kraken2_direct_fragments
kraken2_total_minimizers
kraken2_distinct_minimizers
taxonomic_screen
final_detection
```

`taxonomic_screen` records whether Kraken2 found Norovirus evidence.

`final_detection` is:

```text
SUPPORTED_GII_VP1
```

when at least one genotype passes the mapping filters, otherwise:

```text
NOT_SUPPORTED
```

## `mapping_qc.tsv`

```text
data/results/competitive_mapping/mapping_qc.tsv
```

Columns:

```text
sample
raw_read_pairs
mapping_read_pairs
retained_pair_fraction
unique_VP1_fragments
assigned_fraction
mean_MAPQ
mean_depth
breadth
```

### `raw_read_pairs`

Number of paired fragments in the original input.

### `mapping_read_pairs`

Number of paired fragments entering Kraken2 and, when supported, BWA-MEM2 after preprocessing.

### `retained_pair_fraction`

```text
mapping_read_pairs / raw_read_pairs
```

### `unique_VP1_fragments`

Total number of qualifying primary paired-fragment assignments across the GII VP1 panel.

The term `unique` does not mean PCR duplicate removal.

### `assigned_fraction`

```text
unique_VP1_fragments / mapping_read_pairs
```

This measures how much of the processed sequencing library contributed to the VP1 competitive-mapping counts. It is not a probability that the sample is positive.

### `mean_MAPQ`

Mean BWA mapping quality among primary properly paired R1 alignments before the MAPQ abundance cutoff.

MAPQ measures alignment ambiguity within the supplied reference space. It is not a probability that a biological genotype call is correct.

### `mean_depth` and `breadth`

These are panel-level metrics across the multi-reference FASTA.

A pure genotype sample can therefore have excellent coverage over its true reference while still having relatively low panel-wide breadth. Genotype interpretation should rely primarily on the individual-reference metrics in `mapping_reference_qc.tsv`.

## `mapping_reference_qc.tsv`

```text
data/results/competitive_mapping/mapping_reference_qc.tsv
```

Columns:

```text
sample
VP1_type
reference
fragments
mean_depth
breadth_1x
breadth_5x
```

This file is used to inspect the behavior of every reference separately and is particularly important when multiple representatives are included for the same genotype.

## `mapping_proportions_filtered.tsv`

```text
data/results/competitive_mapping/mapping_proportions_filtered.tsv
```

Contains one column per VP1 genotype.

Only genotypes that pass the abundance and breadth criteria contribute to the final renormalized proportions.

The current streamlined workflow does not generate separate raw-proportion, mapping-count, genotype-QC, or bootstrap-confidence tables.

---

# Reference panel construction

## Nextstrain VP1 dataset

Reference selection began with the Nextstrain norovirus VP1 dataset.

The Nextstrain sequence collection can be downloaded with:

```bash
wget --show-progress \
  https://data.nextstrain.org/files/workflows/norovirus/sequences.fasta.zst \
  -O tmp/nextstrain_sequences.fasta.zst
```

VP1 metadata and the VP1 phylogenetic tree are obtained from the Nextstrain/Auspice data used by the reference-building workflow.

The sequences are separated by VP1 genotype/lineage so that each lineage is treated independently during reference selection.

Examples include:

```text
GII.1
GII.2
GII.3
GII.4
GII.6
GII.17
```

---

## Medoid selection

A medoid is selected for each VP1 lineage.

The medoid is the sequence with the smallest total patristic distance to the other sequences in that lineage. In other words, it is a phylogenetically central representative rather than an arbitrary accession.

Conceptually:

```text
all sequences in one VP1 lineage
            |
            v
calculate phylogenetic distances
            |
            v
identify the most central eligible genome
```

Formally:

```text
medoid = argmin_i sum_j d(i,j)
```

where `d(i,j)` is the patristic distance between two VP1 tree tips.

---

## Completeness filtering

Reference candidates are filtered for completeness before a final genome is retained as the lineage representative.

The reference-selection workflow has used a minimum genome-length requirement of:

```text
7000 nt
```

This reduces the chance of choosing a partial genome that does not contain a complete VP1 region.

Sequences within the lineage can still contribute to the phylogenetic distribution used to define the medoid, while sufficiently complete genomes are used as final reference candidates.

---

## VP1 extraction from annotated GenBank records

After a medoid accession is selected, the corresponding annotated GenBank record is downloaded.

The VP1 region is retained according to the annotation in the `.gb` file. The extracted VP1 coding region is written to an individual FASTA and later combined with the other lineage representatives.

Using VP1 specifically is important because this pipeline is intended to estimate **VP1 capsid genotypes**, not whole-genome strain proportions.

---

## Adding second and third representatives

The original reference strategy used one medoid per VP1 lineage.

During validation against samples with known genotype composition, some single-medoid references did not produce good coverage or did not recover the known lineage proportion accurately. In those cases, one medoid was not sufficient to represent the within-lineage sequence diversity relevant to the sample.

A second representative was therefore added for that genotype. If coverage or known-proportion recovery was still inadequate, a third representative was added.

For example:

```text
GII.17
├── rep1
├── rep2
└── rep3
```

The additional sequences remain separate references during BWA-MEM2 alignment, but their qualifying counts are collapsed back into one biological genotype.

This approach attempts to reduce reference bias without reporting each representative as a different genotype.

---

## Final multi-reference FASTA

All selected VP1 references are combined into:

```text
data/references/reference_multi.fasta
```

Headers use the format:

```text
>GROUP_VP1_type|VP1.type|accession|rep#
```

For example:

```text
>GROUP_GII_17|GII.17|PX470709|rep2
```

The fields represent:

```text
GROUP_GII_17   filesystem-safe group identifier
GII.17         VP1 genotype
PX470709       accession
rep2           representative number
```

`run_competitive_mapping.sh` uses the second pipe-delimited field to group multiple references back into one genotype.

---

# Validation with Bygul mock samples

The effectiveness of the workflow was tested using synthetic mock samples generated with Bygul.

Mock data are useful because the source genomes and expected proportions are known before the analysis is run. The observed competitive-mapping proportions can therefore be compared directly against a defined truth set.

These tests were used to evaluate:

- whether the correct genotype was recovered;
- whether observed proportions matched the known composition;
- whether one medoid adequately represented the lineage;
- whether adding a second or third representative improved coverage and proportion recovery;
- how reference choice influenced cross-mapping, depth, and breadth;
- how low-abundance genotypes behaved in mixtures.

An example Bygul command used for a two-genotype mixture is:

```bash
bygul simulate-proportions \
  source_fastas/GII2_fasta,source_fastas/GII17_fasta \
  --proportions 0.50,0.50 \
  --outdir sample_43 \
  --simulation_mode metagenomics \
  --readcnt 100000 \
  --error_rate 0.001
```

The source FASTAs are supplied in the same positional order as the proportions.

For example:

```text
source_fastas/GII2_fasta,source_fastas/GII17_fasta
0.50,0.50
```

represents an expected mixture of:

```text
GII.2  = 50%
GII.17 = 50%
```

`--simulation_mode metagenomics` creates a mixed metagenomic-style simulation from the supplied source genomes.

`--readcnt 100000` sets the requested simulation size.

`--error_rate 0.001` introduces a 0.1% simulated nucleotide error rate.

Other mock compositions can be generated by changing the source FASTAs and `--proportions` values while retaining the same general command structure.

The mock results were also used to identify cases where the original medoid alone did not sufficiently represent the lineage, motivating the addition of second or third reference representatives.

---

# Interpreting the workflow

The workflow combines several forms of evidence:

```text
technical read QC
        +
metagenomic taxonomic evidence
        +
competitive VP1 mapping
        +
mapping quality
        +
reference-specific depth
        +
reference-specific breadth
        =
supported GII genotype composition
```

No single metric should be interpreted as definitive by itself.

Important interpretation points include:

- a low VP1 `assigned_fraction` can still be compatible with real wastewater Norovirus because viral reads may represent a small fraction of the total metagenome;
- high fragment counts with narrow breadth can reflect localized cross-mapping;
- broad coverage with extremely low abundance can still represent low-level background or cross-mapping;
- MAPQ measures alignment confidence relative to the supplied references, not taxonomic certainty;
- Kraken2 results depend on the database and database version used;
- filtered GII proportions are relative mixture proportions, not absolute viral load.

---

# Current limitations

The current workflow is a method-development pipeline rather than a validated diagnostic assay.

Current limitations include:
- Kraken2 detection depends on database content and version.
- The abundance and breadth thresholds are empirical development settings.
- Relative genotype proportions are not equivalent to absolute virus concentration in the original sample.
- The workflow targets VP1 capsid genotypes and does not perform full polymerase/capsid dual typing.
- The current workflow does not calculate bootstrap confidence or formal confidence intervals.

Continued validation should use known-positive samples, environmental controls, dilution series, and synthetic mixtures with known truth.

---

## Developed by

Aron Asher Diamond, M.S.  
APHL Fellow, New Hampshire Public Health Laboratory