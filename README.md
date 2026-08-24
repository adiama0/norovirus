# Norovirus GII VP1 Competitive Mapping Workflow

## Overview

This repository contains a workflow for estimating the relative abundance of norovirus GII VP1 types from mixed sequencing samples.

The workflow:

1. Obtains norovirus sequence data and VP1 metadata.
2. Identifies one phylogenetic medoid reference for each individual GII VP1 type.
3. Restricts final medoid reference selection to sufficiently complete genomes.
4. Extracts VP1 from each selected medoid genome.
5. Combines all VP1 medoids into a single multifasta reference.
6. Competitively maps paired-end reads against all GII VP1 references simultaneously.
7. Calculates VP1-type-specific fragment counts and proportions.
8. Reports mapping QC metrics including assigned fraction and mean MAPQ.

The workflow was validated using mock metagenomic samples generated with Bygul with known VP1 genotype proportions.

---

# Workflow structure

The main workflow consists of two scripts:

```text
run_medoid_reference.sh
        |
        v
data/references/reference_multi.fasta
        |
        v
run_competitive_mapping.sh
        |
        v
data/competitive_mapping_results/
```

`run_medoid_reference.sh` creates the GII VP1 multifasta reference.

`run_competitive_mapping.sh` maps sequencing reads competitively against the multifasta and calculates VP1-type proportions.

---

# Directory structure

```text
norovirus/
|
├── data/
│   ├── norovirus_all_VP1.json
│   │
│   ├── references/
│   │   ├── medoids.tsv
│   │   ├── medoid_GII_*.gb
│   │   ├── medoid_GII_*_VP1.fasta
│   │   └── reference_multi.fasta
│   │
│   ├── sequences/
│   │   ├── sample_33/
│   │   ├── sample_34/
│   │   ├── ...
│   │   └── sample_50/
│   │
│   ├── sample_proportions.xlsx
│   │
│   └── competitive_mapping_results/
│       ├── bams/
│       ├── group_bams/
│       ├── competitive_mapping_counts.tsv
│       ├── competitive_mapping_proportions.tsv
│       └── competitive_mapping_qc.tsv
│
├── scripts/
│   ├── downsample_seqs.py
│   ├── find_tree_medoids.py
│   └── extract_medoid_vp1.py
│
├── run_medoid_reference.sh
└── run_competitive_mapping.sh
```

---

# Software requirements

The workflow currently uses:

- Bash
- Python 3
- Biopython
- pandas
- BarcodeForge
- BWA-MEM2
- samtools
- seqkit
- wget
- curl
- zstd
- Bygul

---

# Step 1: Build the VP1 medoid reference

Run:

```bash
bash run_medoid_reference.sh
```

The final reference is:

```text
data/references/reference_multi.fasta
```

---

## Download Nextstrain sequences

The workflow downloads the current norovirus sequence FASTA:

```bash
wget --show-progress \
  https://data.nextstrain.org/files/workflows/norovirus/sequences.fasta.zst \
  -O tmp/nextstrain_sequences.fasta.zst
```

The compressed file is decompressed to:

```text
tmp/nextstrain_sequences.fasta
```

---

## Extract VP1 metadata and tree

BarcodeForge extracts the VP1 metadata and phylogenetic tree from:

```text
data/norovirus_all_VP1.json
```

The metadata are written to:

```text
tmp/metadata.tsv
```

The VP1 tree is written to:

```text
tmp/VP1_tree.nwk
```

The tree is not reconstructed during this workflow. It is extracted from the existing Auspice JSON.

---

## Filter VP1 types

Only clean GII VP1 labels matching:

```text
GII.[number]
```

are retained.

For example:

```text
GII.1
GII.2
GII.3
GII.4
...
```

The filtered metadata are stored in:

```text
tmp/metadata_filtered.tsv
```

---

## Count samples per VP1 type

The number of sequences assigned to each VP1 type can be checked with:

```bash
awk -F'\t' '
NR > 1 {
  count[$2]++
}
END {
  for (type in count) {
    print type, count[type]
  }
}
' tmp/metadata_filtered.tsv \
| sort -t. -k2,2n
```

---

# Medoid selection

One medoid is selected independently for every VP1 type.

For example:

```text
GII.1  -> one GII.1 medoid
GII.2  -> one GII.2 medoid
GII.3  -> one GII.3 medoid
...
```

The medoid is the sequence minimizing the total patristic distance to all sequences of the same VP1 type in the VP1 phylogenetic tree.

Mathematically:

```text
medoid = sequence with the smallest total tree distance
         to every other sequence in that VP1 type
```

More formally:

```text
m = argmin sum d(i,j)
```

where `d(i,j)` is the patristic distance between two VP1 tree tips.

---

## Full-genome candidate requirement

All sequences assigned to a VP1 type contribute to determining the phylogenetic center.

However, only genomes meeting the minimum genome-length requirement are allowed to become the final medoid reference.

Current threshold:

```text
7000 nt
```

The relevant argument is:

```bash
--min-genome-length 7000
```

This avoids selecting partial sequences that may not contain a complete VP1 region.

Importantly:

```text
partial genomes
    -> contribute to determining phylogenetic centrality

full genomes >= 7000 nt
    -> are eligible to become the final reference
```

The medoid is therefore the most central sufficiently complete genome, rather than simply the most central sequence of any length.

---

# `find_tree_medoids.py` variables

## `selected_names`

All sequences belonging to a specific VP1 type.

These sequences contribute to the patristic-distance calculation.

Example:

```text
all GII.13 tree tips
```

---

## `candidate_names`

Sequences that are eligible to become the medoid.

Candidates must:

1. belong to the correct VP1 type;
2. exist in the VP1 tree;
3. meet the minimum genome-length requirement.

---

## `sequence_lengths`

Dictionary containing the length of each sequence in the downloaded Nextstrain FASTA.

Example:

```python
{
    "LC122751": 7547,
    "AB809989": 2184
}
```

---

## `number_selected`

Number of sequences belonging to the VP1 type that were found in the VP1 tree.

---

## `number_candidates`

Number of sufficiently complete genomes eligible to become the medoid.

---

## `mean_distance`

Mean patristic distance between the selected medoid and the other sequences belonging to that VP1 type.

A lower value indicates that the reference is closer, on average, to the other members of that VP1 type.

---

# `medoids.tsv`

Medoid selection results are written to:

```text
data/references/medoids.tsv
```

The table contains fields such as:

```text
group
medoid_accession
VP1_type
group_tree_tips
full_genome_candidates
medoid_genome_length
mean_patristic_distance
```

### `group`

Filesystem-safe label for the VP1 type.

Example:

```text
GII_13
```

corresponds to:

```text
GII.13
```

### `medoid_accession`

Accession selected as the phylogenetic medoid reference.

### `VP1_type`

VP1 genotype associated with the medoid.

### `group_tree_tips`

Number of sequences used to define that genotype's phylogenetic distribution.

### `full_genome_candidates`

Number of sequences meeting the full-genome length requirement.

### `medoid_genome_length`

Length of the selected medoid genome.

### `mean_patristic_distance`

Mean VP1-tree distance between the selected medoid and the other sequences assigned to the VP1 type.

---

# VP1 extraction

The GenBank record for each selected medoid accession is downloaded from NCBI.

The workflow first searches annotated CDS features for terms including:

```text
VP1
ORF2
major capsid
capsid protein
```

If an annotated VP1 CDS is found, that sequence is extracted.

If no usable VP1 annotation is present, the script searches for a complete forward-strand ORF within the expected VP1 size range:

```text
1500-1800 nt
```

VP1 sequences are checked for expected length and coding-frame consistency before being written to FASTA.

---

# Final reference multifasta

All individual VP1 medoids are concatenated into:

```text
data/references/reference_multi.fasta
```

Headers have the form:

```text
>GROUP_GII_13|GII.13|ACCESSION
```

For example:

```text
>GROUP_GII_13|GII.13|MW305678
```

The three fields represent:

```text
GROUP_GII_13
    filesystem-safe genotype identifier

GII.13
    VP1 type

MW305678
    medoid accession
```

---

# Step 2: Competitive mapping

Run:

```bash
bash run_competitive_mapping.sh
```

The multifasta reference is:

```bash
REFERENCE="data/references/reference_multi.fasta"
```

The reference is indexed with:

```bash
bwa-mem2 index "${REFERENCE}"
samtools faidx "${REFERENCE}"
```

---

# Reference variables

## `REF_NAMES`

Array containing the exact FASTA reference names.

Example:

```text
GROUP_GII_1|GII.1|ACCESSION
GROUP_GII_2|GII.2|ACCESSION
GROUP_GII_3|GII.3|ACCESSION
```

These exact strings are required by samtools.

---

## `VP1_TYPES`

Array containing the VP1 genotype associated with each reference.

Example:

```text
GII.1
GII.2
GII.3
```

`REF_NAMES` and `VP1_TYPES` use the same array position.

For example:

```text
REF_NAMES[0] = GROUP_GII_1|GII.1|ACCESSION
VP1_TYPES[0] = GII.1
```

---

# Sample input variables

Samples are currently processed from:

```text
sample_33
```

through:

```text
sample_50
```

For each sample:

```bash
directory="data/sequences/sample_${number}"
sample="sample_${number}"
```

Read files are expected to be:

```text
data/sequences/sample_33/sample_33_R1.fastq
data/sequences/sample_33/sample_33_R2.fastq
```

Variables:

```bash
R1="${directory}/${sample}_R1.fastq"
R2="${directory}/${sample}_R2.fastq"
```

---

# Competitive alignment

Reads are aligned simultaneously against every GII VP1 medoid:

```bash
bwa-mem2 mem \
  -t 8 \
  -R "@RG\tID:${sample}\tSM:${sample}\tPL:ILLUMINA" \
  "${REFERENCE}" \
  "${R1}" \
  "${R2}" \
| samtools sort \
    -@ 4 \
    -o "${BAM}"
```

Because every VP1 reference is in the same multifasta, BWA evaluates the references competitively.

The sorted BAM is stored at:

```text
data/competitive_mapping_results/bams/
```

---

# Counting fragments

For each VP1 reference:

```bash
samtools view \
  -c \
  -q 20 \
  -f 66 \
  -F 2308 \
  "${BAM}" \
  "${ref_name}"
```

is used to count qualifying fragments.

---

## `-q 20`

Require:

```text
MAPQ >= 20
```

Alignments below this mapping-confidence threshold are excluded from genotype abundance calculations.

---

## `-f 66`

Require SAM flags:

```text
64 = first read in pair
2  = properly paired
```

Therefore:

```text
66 = first read in pair + properly paired
```

Only R1 is counted so that a paired fragment is counted once rather than twice.

---

## `-F 2308`

Exclude:

```text
4    unmapped
256  secondary alignment
2048 supplementary alignment
```

This prevents secondary and supplementary records from being counted as additional fragments.

---

# Count variables

## `COUNTS`

Array containing the number of qualifying fragments assigned to each VP1 type.

---

## `count`

Number of qualifying fragments assigned to the current VP1 medoid.

---

## `total_unique`

Sum of qualifying fragments across all VP1 references.

Conceptually:

```text
total_unique =
GII.1 count +
GII.2 count +
GII.3 count +
...
```

The name is retained from earlier versions of the workflow.

Here, "unique" refers to fragments passing the competitive mapping and MAPQ filters. It does not refer to PCR duplicate removal.

Mock reads are simulated without PCR amplification, so duplicate marking is not performed.

Duplicate handling should be evaluated separately before applying the workflow to real sequencing libraries.

---

# Proportion calculation

For each genotype:

```text
proportion =
genotype fragment count /
total qualifying VP1 fragments
```

For example:

```text
GII.4 count = 9000
GII.2 count = 1000
total = 10000
```

gives:

```text
GII.4 = 0.90
GII.2 = 0.10
```

---

# QC variables

## `input_pairs`

Number of original paired fragments in the input FASTQ.

Calculated from R1:

```bash
awk 'END {print NR / 4}' "${R1}"
```

FASTQ contains four lines per read.

---

## `not_unique`

Calculated as:

```text
input_pairs - total_unique
```

This value includes all input fragments that did not contribute to the final genotype counts.

It can include:

- reads outside VP1;
- unmapped reads;
- improperly paired reads;
- secondary/supplementary records;
- reads below MAPQ 20.

Therefore, `Not_uniquely_assigned` should not be interpreted exclusively as reads that were ambiguous between two VP1 genotypes.

---

## `assigned_fraction`

Calculated as:

```text
total_unique / input_pairs
```

This measures the fraction of all simulated read pairs that contribute to the final VP1 abundance calculation.

Because the mock reads are generated from whole norovirus genomes while the reference contains only VP1, an assigned fraction near the fraction of the genome occupied by VP1 is expected for well-represented whole-genome mocks.

---

## `Mean_MAPQ`

Mean mapping quality among primary, properly paired R1 alignments before the MAPQ 20 abundance cutoff.

This provides a sample-level measure of alignment confidence.

MAPQ is not the same as the simulated sequencing base-error rate.

`--error_rate 0.001` controls sequencing errors introduced into the mock reads.

MAPQ describes BWA's confidence in the placement of a read relative to competing references.

A read generated with a base-error rate corresponding approximately to Q30 can still have:

```text
MAPQ 60
```

if its reference placement is unambiguous.

Conversely, an error-free read can have:

```text
MAPQ 0
```

if multiple references provide equally plausible alignments.

---

# Output files

## Competitive mapping counts

```text
data/competitive_mapping_results/competitive_mapping_counts.tsv
```

Example structure:

```text
sample  GII.1  GII.2  GII.3  ...  Total_unique_fragments
```

Values are raw qualifying fragment counts.

---

## Competitive mapping proportions

```text
data/competitive_mapping_results/competitive_mapping_proportions.tsv
```

Example:

```text
sample  GII.1  GII.2  GII.3 ...
```

Values are normalized proportions among all qualifying VP1 fragments.

Each sample should approximately sum to:

```text
1.0
```

when at least one VP1 fragment is assigned.

---

## Competitive mapping QC

```text
data/competitive_mapping_results/competitive_mapping_qc.tsv
```

Columns:

```text
sample
Input_read_pairs
Unique_VP1_fragments
Not_uniquely_assigned
Assigned_fraction
Mean_MAPQ
```

---

# Bygul mock samples

Mock samples were generated using Bygul.

General command:

```bash
bygul simulate-proportions \
  source_fastas/GII2_fasta,source_fastas/GII17_fasta \
  --proportions 0.50,0.50 \
  --outdir sample_43 \
  --simulation_mode metagenomics \
  --readcnt 100000 \
  --error_rate 0.001
```

---

## Bygul parameters

### Input FASTAs

```text
source_fastas/GII2_fasta,source_fastas/GII17_fasta
```

Comma-separated genomes used to construct the mock.

---

### `--proportions`

Example:

```bash
--proportions 0.50,0.50
```

Defines the expected abundance of each input genome.

The values correspond positionally to the input FASTAs.

For example:

```text
GII2_fasta,GII17_fasta
0.50,0.50
```

means:

```text
GII.2  = 50%
GII.17 = 50%
```

---

### `--outdir`

Output directory for the simulated sample.

Example:

```bash
--outdir sample_43
```

---

### `--simulation_mode metagenomics`

Simulates a mixed metagenomic sample from the supplied input genomes.

---

### `--readcnt 100000`

Generate:

```text
100,000 read pairs
```

for the sample.

---

### `--error_rate 0.001`

Simulated nucleotide sequencing error rate.

```text
0.001 = 0.1% error probability
```

This is conceptually different from BWA MAPQ.

---

# Expected mock compositions

Expected mock proportions are stored in:

```text
data/sample_proportions.xlsx
```

---

Developed by: Aron Asher Diamond, M.S, APHL Fellow for the New Hampshire Public Health Laboratory. 