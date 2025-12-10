# Cardamom De Novo Assembly for lncRNA Prediction Using Snakemake

This repository contains a fully automated and reproducible Snakemake workflow designed for the discovery of long non-coding RNAs (lncRNAs) in Elettaria cardamomum. The pipeline performs end-to-end processing starting from raw FASTQ reads to final high-confidence lncRNA predictions, integrating multiple tools and custom Python filtering logic.

# Key Features

Automated execution using Snakemake

Raw-read QC (FastQC + Fastp)

De novo transcriptome assembly using RNA-SPAdes

Assembly quality assessment with BUSCO

Redundancy removal using CD-HIT

Transcript quantification using Salmon

Coding potential prediction using CPC2 and PLEK

Final accurate filtering using a custom Python script

Reproducible, modular, and conda-managed environments

# Input Requirements
Input Reads

Paired-end FASTQ files retrieved from NCBI SRA.

The accession list must be provided in:

data/sra.txt


where each line contains an SRA accession ID (e.g., SRRXXXXXXX).

Automatic Downloading

The workflow can fetch FASTQ reads from SRA automatically (if enabled), or you can provide pre-downloaded FASTQ files.

# Workflow Overview

The pipeline contains the following major steps:

FastQC – Raw read quality assessment

Fastp – Adapter trimming and quality filtering

RNA-SPAdes – De novo transcriptome assembly

CD-HIT – Clustering transcripts to remove duplicates

BUSCO – Assembly completeness check

Salmon – Transcript abundance quantification

CPC2 – Coding potential calculation

PLEK – Coding potential prediction using k-mer features

Python Filtering Script – Final selection of confident lncRNAs

# Tools Used (with typical versions)
Snakemake	Workflow management	≥8.0

FastQC	Quality check	0.11.9

Fastp	Read trimming	0.23.4

RNA-SPAdes	De novo assembly	3.15.5

CD-HIT	Clustering transcripts	4.8.1

BUSCO	Assembly completeness	5

Salmon	Quantification	1.10

CPC2	Coding potential	0.1

PLEK	Coding potential	1.2

Python 3	Filtering script	≥3.8

# Reproducibility

All tools are installed via conda.

Pipeline is fully version-controlled.

Uses Snakemake to ensure deterministic execution.

# How to Run the Workflow
1. Install Snakemake with Conda

conda create -n cardamom-lncRNA snakemake -c bioconda -c conda-forge
conda activate cardamom-lncRNA

2. Clone this repository

git clone https://github.com/sajeev26/cardamom_denovo.git

cd cardamom_denovo

3. Add SRA accessions

Place SRA IDs in:

data/sra.txt

4. Run the workflow

snakemake --use-conda --cores 8


# Final Output

The main output file is:

results/lncRNA/confident_lncRNA.csv

This table contains final high-confidence lncRNA transcripts filtered using:

> CPC2 label

> PLEK label & score

> Coding probability

> ORF integrity

> Peptide length

> Fickett score

> Duplicate removal

The FASTA sequences of retained lncRNAs can also be extracted using an included Python script.

# For your queries and suggestions
Reach out to 'sajeevrajssr@gmail.com'
