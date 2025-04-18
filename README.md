---
title: "Documentation"
date: "2025-04-18"
author: "CORRE Guillaume"
output: 
    rmdformats::readthedown:
      keep_md: yes
      highlight: kate
      toc: 3
bibliography: references.bib
---



# Introduction

# Installation

In order to get a working environment, we recommend to clone the git repository :

``` bash
git clone https://github.com/gcorre/GNT_GuideSeq
```

Your folder architecture should look similar to :

``` bash
Path/to/guideseq/
├── 00-pipeline/
├── 01-envs/
├── 02-resources/
├── test/
```

## Conda environments

Install the miniconda environment manager :

``` bash
## https://www.anaconda.com/docs/getting-started/miniconda/install

wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

bash ~/Miniconda3-latest-Linux-x86_64.sh

source ~/.bashrc

conda update -n base -c conda-forge conda

conda install mamba -c conda-forge # faster packages manager

mamba install snakemake # worflow manager 
```

Create the environment in your favorite path (-p) from the `environment.yaml` file:

``` bash
mamba env create -p path/to/env/ -f 01-envs/environment.yaml
```

You should now have an environment named 'guideseq' containing all the required programs when running:

``` bash
#list environments
mamba env list 

# programs version in guideseq environment
mamba list -n guideseq
```

## Reference genomes

The pipeline uses the bowtie2 program to align reads on the reference genome (@langmead2018).

For efficient genome alignment, you can use a pre-built index for Bowtie2. These indices can be downloaded from the [Bowtie2 website](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml). Using a pre-built index saves computational resources and time, as it eliminates the need to build the index from scratch.

If you decide to build the index, we recommend to **not** use haplotypes, scaffolds and unplaced chromosomes to avoid unnecessary multihits alignments. Instead, use the "primary Assembly" version from ensembl or gencode for example. Manually remove unwanted chromosomes if necessary using the seqkit program:

``` bash
conda activate guideseq # load the 'guideseq' environment

seqkit grep  -r -p '^[chr]?[0-9XYM]+' Homo_sapiens.GRCh38.dna.primary_assembly.fa > your_clean_fasta

conda activate           # return to "base" environment
```

Then use the bowtie2 command line:

``` bash
bowtie2-build -@ {threads} {your_clean_fasta} {index_prefix}
```

The index prefix path will be used in the configuration file.

## Annotation

Off-Target site annotation is performed from a GTF file that can be downloaded from any source (ensembl, gencode ...). The GTF file will be processed to keep only gene and exon features for annotation.

Annotation and reference genome **must** use the same chromosome nomenclature (UCSC or ensembl style, ie: "chr1" or "1").

## Prediction tool

The pipeline uses the SWOffinder program (@yaish2024) to predict gRNA off-targets on the fly and annotate OT accordingly as predicted or not.

-   Install from the github repository : <https://github.com/OrensteinLab/SWOffinder>
-   Set the path to SWOffinder in the config file ([Prediction])

# Running the pipeline

In order to run the analysis, 4 elements are mandatory:

1.  The conda environment (see above for installation)
2.  The Sample Data Sheet
3.  The configuration file
4.  The input un-demultiplexed reads in fastq.gz format (R1,R2,I1,I2). *If working with libraries that were already demultiplexed, please read section .....*

## Prepare samples data sheet (SDS) {#prepare-samples-data-sheet-sds}

The sample data-sheet (SDS) is a simple delimited file ( ; ) that contains information about each sample to process in the run. An example is proposed in `./test/sampleInfo.csv`.

| sampleName | CellType | Genome | gRNA_name | gRNA_sequence | orientation | Cas | PAM_sequence | Cut_Offset | Protocole | index1 | index2 |
|------|------|------|------|------|------|------|------|------|------|------|------|
| VEGFA_s1_K562_pos | K562 | GRCh38 | VEGFAs1 | GGGTGGGGGGAGTTTGCTCC | positive | Cas9 | NGG | -4 | guideseq | AGGCAGAA | CTAAGCCT |
| VEGFA_s1_K562_neg | K562 | GRCh38 | VEGFAs1 | GGGTGGGGGGAGTTTGCTCC | negative | Cas9 | NGG | -4 | guideseq | TCCTGAGC | CTAAGCCT |

: Mandatory columns are:

-   **`sampleName`** : Sample name to use.
    -   Will be use to name output files and folders.
    -   Samples with the same `sampleName` will be merged before processing ([Merging samples])
    -   Samples name [**must not include**]{.underline} '-' in their name as this symbol is used in the pipeline for a special purpose. Instead, use "\_" or any other separator of your choice.
-   **`Genome`**: Uses to define which reference genome to use for each sample.
    -   This **must** be one of the values present in config file **`genome`** key (see [Reference genome]).
-   **`gRNA_name`** : name of the gRNA used.
-   **`gRNA_sequence`**: Sequence of the gRNA [**without**]{.underline} the PAM sequence.
-   **`orientation`** : which PCR orientation was chosen ["positive", "negative", "mix" if they share the same indexes]
    -   (this info is only used for metadata purposes, both PCR orientations are automatically processes by the pipeline).
-   **`Cas`**: name of the Cas used
-   **`PAM_sequence`**: Sequence of the PAM [default : "NGG"]
-   **`Cut_Offset`**: Distance from gRNA end where the cut occurs [default : -4]
-   **`type`**: Type of experiment ["guideseq" or "iguideseq"]. This value will define which sequence to trim ([Reads adapter & ODN trimming sequences](#reads-adapter-odn-trimming-sequences)).
-   **`index1`**: Sequence of index 1
-   **`index2`**: Sequence of index 2

Additional columns can be added for metadata annotation purpose.

SDS format is validated using the `snakemake.utils - validate` function.

### Merging samples

In certain scenarios, it may be beneficial to merge different samples from a library during the reads processing stage. This can be particularly useful when dealing with Multiple Replicates or +/-PCR orientations.

To achieve this, you can use the same sample name for multiple rows in the Sample Description Sheet (SDS). Samples that share the same name will be merged during the reads processing, provided they meet the following criteria:

1.  **`Reference Genome`**: The samples must use the same reference genome.

2.  **`gRNA`**: The samples must have the same gRNA, both in terms of sequence and name.

3.  **`Cas Protein`**: The samples must use the same Cas protein, with identical PAM (Protospacer Adjacent Motif) and Offset values.

4.  **`Protocole`**: The samples must use the same ODN (Oligodeoxynucleotide), defined as guideseq or iguideseq.

If any of the above conditions are not met, an error will be raised, and the pipeline will be stopped. This ensures that only compatible samples are merged, maintaining the integrity of the data processing workflow.

Using the test SDS above, if you want to merge both positive and negative libraries, give the same sample name to both rows. As they use the same genome, gRNA, Cas & method, they will be aggregated in a single library.

| sampleName | CellType | Genome | gRNA_name | gRNA_sequence | orientation | Cas | PAM_sequence | Cut_Offset | Protocole | index1 | index2 |
|------|------|------|------|------|------|------|------|------|------|------|------|
| VEGFA_s1_K562 | K562 | GRCh38 | VEGFAs1 | GGGTGGGGGGAGTTTGCTCC | positive | Cas9 | NGG | -4 | guideseq | AGGCAGAA | CTAAGCCT |
| VEGFA_s1_K562 | K562 | GRCh38 | VEGFAs1 | GGGTGGGGGGAGTTTGCTCC | negative | Cas9 | NGG | -4 | guideseq | TCCTGAGC | CTAAGCCT |

: UMI and reads counts will be breakdown to each PCR orientation in the final result table.

#### Working with already demultiplexed libraries

To be documented

## Prepare configuration file

The configuration file is a yaml formatted file with key-value dictionary used to fine-tune the pipeline behavior. Settings will apply to [**all samples of the run**]{.underline}.

An example of configuration file is proposed in `./test/guideSeq.yaml.`

[config file must be present in working directory when starting the pipeline.]{.underline}

### Metadata

Author and affiliation will be printed on the final report.

``` yaml
author: "Guillaume CORRE, PhD"
affiliation: "Therapeutic Gene Editing - GENETHON, INSERM U951, Paris-Saclay University, Evry, France"
```

### Path to Sample Data Sheet and read files

SDS path is relative to current folder (from where the pipeline is started). Absolute path can be used.

``` yaml
## Library informations
sampleInfo_path: "sampleInfo.csv"
read_path: "."  # path to reads if not current folder
R1: "Undetermined_S0_L001_R1_001.fastq.gz"
R2: "Undetermined_S0_L001_R2_001.fastq.gz"
I1: "Undetermined_S0_L001_I1_001.fastq.gz"
I2: "Undetermined_S0_L001_I2_001.fastq.gz"
```

### Reference genome

Genome name is user-defined but must be referenced using the exact same name in the Sample Data Sheet.

'Index' corresponds to the prefix used when bowtie2 index was built.

``` yaml
## path to references
genome:
  GRCh37:
    fasta: "/PATH_TO_REFERENCE/GRCh37/Sequences/GRCh37.primary_assembly.genome.fa"
    index: "/PATH_TO_REFERENCE/GRCh37/Indexes/Bowtie2/GRCh37.primary_assembly.genome" 
    annotation: "/PATH_TO_REFERENCE/GRCh37/Annotations/gencode.v19.annotation.gtf.gz"
  GRCh38:
    fasta: "/PATH_TO_REFERENCE/GRCh38/Sequences/GRCh38.primary_assembly.genome.fa"
    index: "/PATH_TO_REFERENCE/GRch38/Indexes/Bowtie2/GRCh38.primary_assembly.genome" 
    annotation: "/PATH_TO_REFERENCE/GRCh38/Annotations/gencode.v46.annotation.gtf.gz"
  GRCm39:
    fasta: "/PATH_TO_REFERENCE/GRCm39/Sequences/GRCm39.primary_assembly.genome.fa"
    index: "/PATH_TO_REFERENCE/GRCm39/Indexes/Bowtie2/GRCm39.primary_assembly.genome" 
    annotation: "/PATH_TO_REFERENCE/GRCm39/Annotations/gencode.vM36.annotation.gtf.gz"
```

### Reads filtering

After adaptor & ODN triming and before alignment to the reference genome, reads that are too short are discarded. This filter is applied if **any** of the mates size is below this threshold.

``` yaml
################################################
minLength: 25 ## Minimal read length after trimming, before alignment
################################################
```

### Alignment to reference genome

Define which aligner to use and the range of fragment size.

``` yaml
## Alignement 
################################################
aligner: "bowtie2"   ## Aligner to use (bowtie2 only for now)
minFragLength: 100         # Minimal fragment length after alignment
maxFragLength: 1500        # Maximal fragment length after alignment 
################################################
```

### Insertion sites calling

After alignment on the genome, reads are collapse to single insertion points and aggregated if they cluster in a distance smaller than `ISbinWindow` defined here. Filters can be applied to exclude UMI with few reads (`minReadsPerUMI`) or insertion sites with few UMIs (`minUMIPerIS`).

Here, you can also define if you want to tolerate bulges in the alignment between the gRNA and gDNA.

``` yaml
################################################
## Off targets calling
tolerate_bulges: "FALSE"           # whether to include gaps in the gRNA alignment (this will change the gap penalty during SW pairwise alignment)
max_edits_crRNA: 6              # filter clusters with less or equal than n edits in the crRNA sequence (edits = substitutions + INDELs)
ISbinWindow: 100                # insertion sites closer than 'ISbinWindow' will be clustered together
minReadsPerUMI: 0               # 0 to keep all UMIs, otherwise min number of reads per UMIs
minUMIPerIS: 0                  # 0 to keep all IS, otherwise min number of UMI per IS
slopSize: 50                    # window size (bp) around IS (both directions) to identify gRNA sequence (ie 50bp = -50bp to +50bp)
################################################
```

### UMI correction

Due to potential sequencing errors, additionnal UMIs may be detected and a correction step is required.

For each cluster of Off Targets, a similarity matrix between all UMIs detected is calculated and similar UMI collapsed together if the editing distance is smaller than `UMI_hamming_distance` . The Adjacency method described in UMI-tools (@Smith2017) is used by default (see <https://umi-tools.readthedocs.io/en/latest/the_methods.html> for details).

``` yaml
################################################
# post alignment
minMAPQ: 1                      # Min MAPQ score to keep alignments -> !!! multihits have a MAPQ close to 1. Value greater than 1 will discard Offtargets with exact same sequence.
UMI_hamming_distance: 1         # min distance to cluster UMI using network-based deduplication, use [0] to keep raw UMIs
UMI_deduplication: "Adjacency"  # method to correct UMI (cluster or Adjacency)
UMI_pattern: "NNWNNWNN"  
UMI_filter: "FALSE"               # If TRUE, remove UMIs that do no match the expected pattern [FALSE or TRUE]
################################################
```

### Reporting

``` yaml
################################################
# reporting
max_clusters: 100                 # max number of cluster alignments to report
minUMI_alignments_figure: 1       # filter clusters with more than n UMI in the report alignment figure (set to 0 to keep all clusters -> can be slow)
################################################
```

### Prediction

``` yaml
# Prediction
################################################
SWoffFinder:
  path: "/opt/SWOffinder" ## Path to SWoffinder on your server (downloaded from https://github.com/OrensteinLab/SWOffinder)
  # maxE: 6                 # Max edits allowed (integer). --> actually use the same value as max_edits_crRNA above for consistency between prediction and filtering of OT
  # maxM: 6                 # Max mismatches allowed without bulges (integer). --> actually use the same value as max_edits_crRNA above for consistency between prediction and filtering of OT
  maxMB: 4                # Max mismatches allowed with bulges (integer).
  maxB: 3                 # Max bulges allowed (integer).
  window_size: 100
################################################
```

### Reads adapter & ODN trimming sequences {#reads-adapter-odn-trimming-sequences}

Indicate which sequence will be trimmed from R1 & R2 reads ends depending on the PCR orientation and ODN used.

``` yaml
################################################
# Sequences for the trimming steps

guideseq:
  positive:
    R1_trailing: "GTTTAATTGAGTTGTCATATGT"
    R2_leading: "ACATATGACAACTCAATTAAAC"
    R2_trailing: "AGATCGGAAGAGCGTCGTGT"
  negative:
    R1_trailing: "ATACCGTTATTAACATATGACAACTCAA"
    R2_leading: "TTGAGTTGTCATATGTTAATAACGGTAT"
    R2_trailing: "AGATCGGAAGAGCGTCGTGT"



iguideseq:
  positive:
    R2_leading: "ACATATGACAACTCAATTAAACGCGAGC"
    R2_trailing: "AGATCGGAAGAGCGTCGTGT"
    R1_trailing: "GCTCGCGTTTAATTGAGTTGTCATATGT"
  negative:
    R1_trailing: "TCGCGTATACCGTTATTAACATATGACAACTCAA"
    R2_leading: "TTGAGTTGTCATATGTTAATAACGGTATACGCGA"
    R2_trailing: "AGATCGGAAGAGCGTCGTGT"
################################################
```

## Folder structure

In order to start a run:

-   create a new directory in the installation folder and `cd` into.

-   Then :

    -   move the sample data-sheet

    -   move the configuration file

    -   move the illumina sequencing `undeterminded_R1/R2/I1/I2.fastq.gz` files (undemultiplexed)

> Input fastq files should respect the following structure from original paper :
>
> -   R1 : contains fragment sequence starting in gDNA
>
> -   R2: Starts with ODN sequence followed by gDNA sequence and potential adaptor sequence
>
> -   i1 : Contains barcode 1 (usually 8 nucleotides)
>
> -   i2 : Contains barcode 2 and UMI sequence (usually 8 + 8 nucleotides)

Your folder architecture should look similar to :

``` bash
Path/to/guideseq/
├── 00-pipeline/
├── 01-envs/
├── 02-resources/
├── test/
├── My/folder/
    ├── guideSeq_GNT.yml
    ├── sampleInfo.csv
    ├── undertermined_R1.fastq.gz
    ├── undertermined_R2.fastq.gz
    ├── undertermined_I1.fastq.gz
    └── undertermined_I2.fastq.gz
```

From inside your analysis folder, run the command below after adjusting for number of CPU (-j) :

``` bash
snakemake -s ../00-pipeline/guideseq_gnt.smk \
  -j 24 \         ## number of threads used
  -k \            ## keep running on rule error
  --use-conda \   ## use conda environment
  -n  \            ## dry-run, will print rules without running them. Remove this argument if no error is returned
  --report-after-run --report run_report.html # produce report with run-time
  
  
# remove the -n argument to start the pipeline
```

Each rule takes a maximum of 6 threads (alignment, trimming) to speed up the data processing. Set `-j` as a multiple of 6 to process 2 (12threads) , 3 (18 threads) ... samples in parallel.

``` bash
## usefull snakemake arguments 
# see https://snakemake.readthedocs.io/en/stable/executing/cli.html#all-options

--quiet
--notemp # keep all intermediate files
--rerun-trigger {code,input,mtime,params,software-env} 
--rerun-incomplete

 --filegraph | sed -n '/digraph/,$p' | dot -Tpdf > dag_files.pdf #DAG representation of workflow
 --rulegraph | sed -n '/digraph/,$p' | dot -Tpdf > dag_rules.pdf 
```

# Output

If everything goes well, the pipeline should end successfully :

``` latex
Workflow finished, no error
 _____________
< You rock baby !! >
 -------------
        \   ^__^
         \  (♥♥)\_______
            (__)\       )\/\
                ||----w |
                ||     ||
```

Otherwise, an error will be raised and the origine of the problem reported by snakemake:

``` latex
< Houston, we have a problem >
 ----------------------------
       \   \_______
 v__v   \  \   O   )
 (xx)      ||----w |
 (__)      ||     ||  \/\
```

Upon pipeline completion, your folder should now look like (test dataset provided in ./test) :

``` latex
my_folder_name/
├── 00-demultiplexing
│   ├── demultiplexing_R1.log
├── 01-trimming
│   ├── VEGFAs1.odn.log
│   ├── VEGFAs1.trailing.log
├── 02-filtering
│   ├── VEGFAs1.filter.log
│   ├── VEGFAs1_R1.UMI.ODN.trimmed.filtered.fastq.gz
│   ├── VEGFAs1_R2.UMI.ODN.trimmed.filtered.fastq.gz
├── 03-align
│   ├── VEGFAs1_multi.txt
│   ├── VEGFAs1_R1.UMI.ODN.trimmed.unmapped.fastq.gz
│   ├── VEGFAs1_R2.UMI.ODN.trimmed.unmapped.fastq.gz
│   ├── VEGFAs1.UMI.ODN.trimmed.filtered.align.log
│   ├── VEGFAs1.UMI.ODN.trimmed.filtered.sorted.filtered.bam
│   ├── VEGFAs1.UMI.ODN.trimmed.filtered.sorted.filtered.bam.bai
├── 04-IScalling
│   ├── VEGFAs1.cluster_slop.bed
│   ├── VEGFAs1.cluster_slop.fa
│   ├── VEGFAs1.reads_per_UMI_per_IS.bed
│   ├── VEGFAs1.reads_per_UMI_per_IS_corrected.bed
│   ├── VEGFAs1.UMIs_per_IS_in_Cluster.bed
├── 05-Report
│   ├── VEGFAs1.rdata
│   ├── VEGFAs1.stat
│   ├── VEGFAs1_summary.tsv
│   ├── VEGFAs1_summary.xlsx
│   ├── report-files
│   │   ├── VEGFAs1_offtargets_dynamic_files/
│   │   ├── VEGFAs1_offtargets_dynamic.html
│   │   ├── VEGFAs1_offtargets.html
│   ├── report.html
│   └── report.rdata
├── 06-offPredict
│   ├── GRCh38_GGGTGGGGGGAGTTTGCTCCNGG.csv
│   └── GRCh38_GGGTGGGGGGAGTTTGCTCCNGG.txt
├── guideSeq_GNT.yml
├── sampleInfo.csv
├── my_folder_name_report.html
├── undertermined_R1.fastq.gz
├── undertermined_R2.fastq.gz
├── undertermined_I1.fastq.gz
└── undertermined_I2.fastq.gz
```

## Report

A general report is generated. It summarizes all main QC and key features obtained from the run using graphical representations and tables.

## Off-targets files

For each sample, an excel file with the complete OT information is generated. This file has many columns among which are of particular interest :

-   

-   

-   

-   

-   

-   

# Pipeline step-by-step explanations

## Demultiplexing

**Tool** : cutadapt

Undetermined fastq files are demultiplexed to `sampleName` fastq files according to barcodes present in the sample data sheet.

-   Barcode1 and barcode2 sequences are concatenated to build the demultiplexing index.

-   i1 and i2 fastq files are concatenated to a single fastq file (i3)

-   R1 and R2 fastq are demultiplexed according to i1+i2 sequence present in i3 fastq files.

-   UMI sequence is added to R1 and R2 read name for future UMI quantification.

    -   UMI is extracted from i3 first nucleotides according to the length of `UMI_pattern` variable in the configuration file.

## ODN trimming:

**Tool** : cutadapt

The leading ODN sequence is remove from R2 reads according to the `method` and `PCR orientation` defined in the sample data sheet for each sample.

Reads that **do not start** with the ODN sequence are discarded.

## Adaptor trimming:

**Tool** : cutadapt

The ODN and adaptor trailing sequences are trimmed from R1 and R2 reads respectively if present. If the DNA fragment is large compared to R1/R2 sequences length, those sequences may not be present.

## Read filtering:

**Tool** : cutadapt

After trimming, only read pairs with **both** mates longer than a `minLength` defined in the configuration file are selected for alignment.

## Genome alignment:

Reference sequences are specified in the configuration file. Genome index is build if it does not already exist.

``` yaml
## path to references
genome:
  human:
    fasta: "/media/References/Human/Genecode/GRch38/Sequences/GRCh38.primary_assembly.genome.fa"
    index: "/media/References/Human/Genecode/GRch38/Indexes/Bowtie2/GRCh38.primary_assembly.genome" # path to index created during the run if not existing yet
    annotation: "/media/References/Human/Genecode/GRch38/Annotations/gencode.v46.annotation.gtf.gz"
```

Trimmed reads pairs that passed the filtering steps are then aligned on the `reference genome` specified in the sample data sheet for each sample using the `aligner` specified in configuration file.

``` yaml
## Alignement 
################################################
aligner: "bowtie2"   ## Aligner to use (bowtie2 or bwa)
minFragLength: 100         # Minimal fragment length after alignment
maxFragLength: 1500        # Maximal fragment length after alignment 
################################################
```

Multihit alignments with low MAPQ score are discarded except if the Alignment score (AS tag) is equal to the score of the second best alignment score (XS tag) .

### Multi-hits management:

Multihits are reads with multiple possible equally good alignment positions in the genome. We choose to keep only a single random alignment for each read instead of reporting all possible positions.

On the long run, if multi read arise from the same cuting site, they will distribute randomly to all sites (explain more).

## Cutting site calling:

Following alignment on the reference genome, nuclease cutting sites are extracted from the start position of R2 read alignment.

Reads that align at the same cutting site with the same UMI are aggregated together keeping track of total reads count.

## UMI correction:

UMI sequences are corrected for potential sequencing error using the parameters defined in the configuration file

``` yaml
UMI_hamming_distance: 1         # min distance to cluster UMI using network-based deduplication, use [0] to keep raw UMIs
UMI_deduplication: "Adjacency"  # method to correct UMI (cluster or Adjacency)
UMI_pattern: "NNWNNWNN"  
UMI_filter: "FALSE"               # If TRUE, remove UMIs that do no match the expected pattern [FALSE or TRUE]
```

## Cut site clustering:

Cut sites than fall in the same window of `ISbinWindow` defined in the configuration file are clustered together.

## gRNA match:

For each cluster of cutting sites, the gRNA sequence defined in the sample data sheet for each sample is looked up in a window of +/- `slopSize` bp using the Swith-Watterman algorithm.

Gap tolerance can be accepted to detect bulges in gDNA or gRNA if the `tolerate_bulges` variable is set.

## Annotation of clusters:

Clusters of cutting sites are annotated using the gtf file specified in the configuration file for each organism.

A first step prepare the annotation file to extract only gene and exons features.

A second step annotate clusters to gene and position relative to those genes (exon/intron). Multiple annotations may be present for each cluster and are reported.

## Off target prediction:

for each gRNA sequence and each organism specified in the sample data sheet, a prediction of OTS is realized using the SWOffinder tool with parameters defined in the configuration file :

``` yaml
# Prediction
################################################
SWoffFinder:
  path: "/opt/SWOffinder" ## Path to SWoffinder on your server (downloaded from https://github.com/OrensteinLab/SWOffinder)
  maxE: 6                 # Max edits allowed (integer).
  maxM: 6                 # Max mismatches allowed without bulges (integer).
  maxMB: 6                # Max mismatches allowed with bulges (integer).
  maxB: 3                 # Max bulges allowed (integer).
  window_size: 100
################################################
```

## Reporting:

A short report is generated with main tables and graphical representations to better understand pipeline results.

# References
