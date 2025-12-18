# EOH IonTorrent WGS Pipeline  

**Nextflow DSL2 – Reference-based assembly pipeline**


## Introduction

This pipeline implements a **complete, reference-based WGS workflow for IonTorrent bacterial data**

Special attention has been given to:
- IonTorrent-specific error profiles (homopolymers, indels)
- generation of *ORF-preserving consensus genomes*
- maximization of loci successfully called by an allele caller
- full reproducibility via Docker containers

The pipeline is written in **Nextflow DSL2** and is designed to be:
- reproducible
- modular
- easy to run on HPC or standalone servers


## Global Workflow Overview

For each sample, the pipeline executes the following steps:

* **MASH** – selection of the closest reference genome  
* **Reference preparation** – copy and normalization of selected reference  
* **Read mapping** – Bowtie2  
* **Alignment processing** – SAM → BAM → sorted & indexed BAM  
* **Consensus reconstruction** – samtools + bcftools + vcfutils  
* **Consensus conversion** – FASTQ → FASTA  

Each step is executed in an isolated Docker container.


## Requirements

### Software
- **Nextflow ≥ 22.04**
- **Docker**

### Verify installation
```bash
nextflow -version
docker --version
docker info | head
```


## Directory Structure

### Project layout (recommended)

```
EOH_iontorrent/
├── EOH_iontorrent_pipeline_light.nf
├── nextflow.config
└── params.json

/path-to-reads/
/path-to-refs/
```


## Input Data

### Reads directory (`reads_dir`)

The pipeline expects **single-end IonTorrent WGS reads**.

Accepted formats:
- `.fastq`
- `.fq`
- `.fastq.gz`
- `.fq.gz`

All reads must be placed **directly in the directory**, not in subfolders.

Example:
```
reads/
├── sample1.fastq.gz
├── sample2.fastq.gz
```

**Important notes**
- File name (without extension) is used as the **sample ID**
- Avoid spaces and special characters in file names


### Reference directory (`refs_dir`)

This directory must contain **candidate reference genomes** in FASTA format.

Accepted extensions:
- `.fa`
- `.fasta`
- `.fna`
- optionally gzipped (`.gz`)

Example:
```
refs/
├── XXXXXXXXXXXX.1.fasta
├── YYYYYYYYYYYY.1.fna
```

The pipeline uses **MASH** to:
- sketch all reference genomes
- compute distances between reads and references
- select the **best reference per sample**


## Configuration (`params.json`)

### Example
```json
{
  "reads_dir": "path-to-reads",
  "refs_dir": "path-to-refs",
  "outdir": "results",
  "threads": 16,
  "memory_gb": 64,
  "genus_species": "listeria_monocytogenes"
}
```

### Parameter description

| Parameter | Description |
||
| `reads_dir` | Directory containing FASTQ reads |
| `refs_dir` | Directory with candidate reference genomes |
| `outdir` | Output directory |
| `threads` | Threads for multithreaded tools |
| `memory_gb` | Memory allocation (GB) |


## Detailed Process Description

### MASH – Reference selection

**Input**
- Trimmed reads in FASTQ format
- Folder with references as FASTA files

**Output**
- Distance table
- Selected reference genome code

**Purpose**
- Automatically select the most appropriate reference for each sample


### Read Mapping (Bowtie2)

**Input**
- Trimmed reads in FASTQ format
- Selected reference genome as fasta file

**Output**
- SAM alignment file

**Purpose**
- Align IonTorrent reads to the closest reference genome


### Alignment processing (SAMTOOLS)

**Input**
- SAM alignment file

**Output**
- Sorted and indexed BAM file

**Purpose**
- Prepare alignment for variant calling and consensus generation


### Consensus reconstruction (samtools + bcftools + vcfutils)

**Input**
- Sorted and indexed BAM file
- Reference genome as fasta file

**Output**
- Consensus as FASTQ format

**Purpose**
- Reconstruct a conservative consensus genome
- Minimize indel-induced frameshifts
- Preserve intact CDS for cgMLST analysis

This step is optimized for **IonTorrent WGS** and is intentionally conservative.


### FASTQ → FASTA conversion

**Input**
- Consensus as FASTQ format

**Output**
- Consensus as FASTA format

**Purpose**
- Produce final genome sequence 

## Output Structure

```
results/SAMPLE/
├── mash/
├── ref/
├── mapping/
└── consensus/
```


## Downloading the Project and Example Data from GitHub

### Clone the repository

```bash
git clone https://github.com/adipi71/EOH_iontorrent.git
cd EOH_iontorrent
```


### Example data included

```
example_refs_listeria.zip
```


### Unzip example data

```bash
unzip example_refs_listeria.zip
```

This creates:

```
dataset/
└── listeria/
```
Put in params.json 

```
  "refs_dir": "dataset/listeria",
  "reads_dir": "path-to-reads",
```

## Run the pipeline

**nextflow version**

```bash
nextflow run EOH_iontorrent_pipeline_light.nf -params-file params.json -profile standard -c no_container.config
```

**docker version**

```bash
docker run --rm -it -u 0:0 -v $(pwd):/work -v /MY/SAMPLE/DIR:/mnt/bioinfonas/ -w /work adipi71/eoh_mapping:latest
```

## Output Structure

For each sample:

```
outdir/SAMPLE/
├── mash/         # reference selection results
├── ref/          # chosen reference genome
├── mapping/      # SAM/BAM/VCF/FASTQ files
└── consensus/    # final consensus FASTA
```

## Dockerfile of adipi71/eoh_mapping

```Dockerfile
FROM continuumio/miniconda3:latest

RUN apt-get update && apt-get install -y --no-install-recommends \
    wget make gcc file \
    zlib1g-dev libbz2-dev libncurses5-dev \
    perl procps rsync ca-certificates \
  && rm -rf /var/lib/apt/lists/*

RUN set -eux; \
    wget -O samtools-0.1.19.tar.bz2 "https://sourceforge.net/projects/samtools/files/samtools/0.1.19/samtools-0.1.19.tar.bz2/download"; \
    file samtools-0.1.19.tar.bz2; \
    tar -xjf samtools-0.1.19.tar.bz2; \
    cd samtools-0.1.19; \
    make; \
    test -x samtools; \
    test -x bcftools/bcftools; \
    test -f bcftools/vcfutils.pl; \
    install -m 0755 samtools /usr/local/bin/samtools; \
    install -m 0755 bcftools/bcftools /usr/local/bin/bcftools; \
    install -m 0755 bcftools/vcfutils.pl /usr/local/bin/vcfutils.pl; \
    cd /; \
    rm -rf samtools-0.1.19 samtools-0.1.19.tar.bz2

RUN set -eux; \
    conda config --system --add channels conda-forge; \
    conda config --system --add channels bioconda; \
    conda config --system --set channel_priority strict; \
    conda install -y mash=2.3 bowtie2=2.5.4 nextflow; \
    conda clean -a -y; \
    which mash; which bowtie2; which nextflow; \
    nextflow -version

RUN pip install --no-cache-dir numpy biopython pandas xlrd

```

