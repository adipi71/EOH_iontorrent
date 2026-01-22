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
- **Docker version 24.0.4**

### Verify installation
```bash
docker --version
```


## Directory Structure

### Project layout (recommended)

```
EOH_iontorrent/
└── params.json

/path-to-reads/
/path-to-refs/
/outputdir/
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

### Optional: precomputed MASH sketch (refs_msh)

If you run the pipeline multiple times using the same reference database, you can skip the reference sketching step by providing a precomputed MASH sketch file (.msh).

Parameter

refs_msh – path to a precomputed refs.msh file

Important

If refs_msh is provided, refs_dir is still required.
The pipeline uses refs_dir to resolve and copy the selected reference FASTA after MASH selection.

If you provide refs_msh without refs_dir, the pipeline will fail immediately with a clear error message.

Typical location of the sketch produced by the pipeline:

## Configuration (`params.json`)

### Example
```json
{
  "reads_dir": "path/input/reads",
  "refs_dir": "path/input/refs",
  "outdir": "results",
  "threads": 16,
  "memory_gb": 64,
}
```
with msh 

```json
{
  "reads_dir": "path/input/reads",
  "refs_dir": "path/input/refs",
  "refs_msh": "path/input/refs.msh",
  "outdir": "results",
  "threads": 16,
  "memory_gb": 64
}
```

### Parameter description

| Parameter | Description |
||
| `reads_dir` | Directory containing FASTQ reads |
| `refs_dir` | Directory with candidate reference genomes |
| `refs_msh`  | *(Optional)* Precomputed MASH sketch (`.msh`) to skip reference sketching (**requires `refs_dir`**) |
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

```bash
docker run -it -u 0:0 --rm \
  -v /path/input/:/path/input/ \
  -v $(pwd):/work \ 
  adipi71/eoh_mapping:0.10 nextflow run /pipeline/EOH_iontorrent_pipeline_light.nf -params-file params.json
```
or if you don't want to use the params.json file

```bash
docker run -it -u 0:0 --rm \
  -v /path/input:/path/input/ \
  -v $(pwd):/work \ 
  adipi71/eoh_mapping:0.10  nextflow run  /pipeline/EOH_iontorrent_pipeline_light.nf \ 
  --reads_dir /path/input/reads \
  --refs_msh  /path/input/refs.msh \
  --refs_dir /path/input/listeria \
  --outdir outdir \
  --threds 16 --memory_gb 64 
```

### Explanation of the mounts

| Option | Meaning |
|------|--------|
| `-v /path/input:/path/input/` | Makes input/output data visible inside the container |
| `-v $(pwd):/work` | Stores Nextflow work directory and logs on the host |
| `-u 0:0` | Runs as root (avoids permission issues) |
| `--rm` | Removes the container when finished |


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
FROM continuumio/miniconda3:25.3.1-1
# Pin Nextflow here (stable modern version; Java 17 friendly)
ARG NXF_VER=25.10.2
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget make gcc file \
    zlib1g-dev libbz2-dev libncurses5-dev \
    perl procps rsync ca-certificates \
    openjdk-17-jre-headless \
  && rm -rf /var/lib/apt/lists/*
# ---- samtools 0.1.19 + bcftools + vcfutils.pl ----
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
# ---- conda tools (no nextflow here) ----
RUN set -eux; \
    conda config --system --add channels conda-forge; \
    conda config --system --add channels bioconda; \
    conda config --system --set channel_priority strict; \
    conda config --system --remove channels defaults || true; \
    conda create -y -n eoh \
        python=3.11 \
        mash=2.3 \
        bowtie2=2.5.4; \
    conda clean -a -y
ENV PATH=/opt/conda/envs/eoh/bin:$PATH
RUN pip install --no-cache-dir numpy biopython pandas xlrd
# ---- Nextflow standalone distribution (pinned) ----
# This is the "dist" nextflow executable. It should not trigger CAPSULE downloads like the conda launcher.
RUN set -eux; \
    wget -qO /usr/local/bin/nextflow \
      "https://github.com/nextflow-io/nextflow/releases/download/v${NXF_VER}/nextflow"; \
    chmod +x /usr/local/bin/nextflow; \
    nextflow -version
# ---- pipeline files ----
WORKDIR /pipeline
COPY EOH_iontorrent_pipeline_light.nf /pipeline/EOH_iontorrent_pipeline_light.nf
COPY nextflow.config                  /pipeline/nextflow.config
COPY entrypoint.sh                    /entrypoint.sh
RUN chmod +x /entrypoint.sh \
 && mkdir -p /work /tmp \
 && chmod 1777 /tmp
# Nextflow runtime dirs
ENV NXF_VER=${NXF_VER} \
    NXF_HOME=/work/.nextflow \
    NXF_WORK=/work/work \
    NXF_TEMP=/tmp \
    NXF_ANSI_LOG=false
ENTRYPOINT ["/entrypoint.sh"]
CMD ["bash"]

```

