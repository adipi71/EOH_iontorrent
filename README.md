# EOH IonTorrent WGS Pipeline  

**Nextflow DSL2 – Reference-based assembly pipeline**


## 1. Introduction

This pipeline implements a **complete, reference-based WGS workflow for IonTorrent bacterial data**

Special attention has been given to:
- IonTorrent-specific error profiles (homopolymers, indels)
- generation of *ORF-preserving consensus genomes*
- maximization of loci successfully called by chewBBACA
- full reproducibility via Docker containers

The pipeline is written in **Nextflow DSL2** and is designed to be:
- reproducible
- modular
- easy to run on HPC or standalone servers


## 2. Global Workflow Overview

For each sample, the pipeline executes the following steps:

1. **FASTP** – read trimming and quality control  
2. **MASH** – selection of the closest reference genome  
3. **Reference preparation** – copy and normalization of selected reference  
4. **Read mapping** – Bowtie2  
5. **Alignment processing** – SAM → BAM → sorted & indexed BAM  
6. **Consensus reconstruction** – samtools + bcftools + vcfutils  
7. **Consensus conversion** – FASTQ → FASTA  
8. **chewBBACA** – cgMLST allele calling  

Each step is executed in an isolated Docker container.


## 3. Requirements

### 3.1 Software
- **Nextflow ≥ 22.04**
- **Docker**

### 3.2 Verify installation
```bash
nextflow -version
docker --version
docker info | head
```


## 4. Directory Structure

### 4.1 Project layout (recommended)

```
EOH_iontorrent/
├── EOH_iontorrent_pipeline_final.nf
├── nextflow.config
└── params.json

/path-to-reads/
/path-to-refs/
```


## 5. Input Data

### 5.1 Reads directory (`reads_dir`)

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


### 5.2 Reference directory (`refs_dir`)

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


## 6. Configuration (`params.json`)

### 6.1 Example
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

### 6.2 Parameter description

| Parameter | Description |
||
| `reads_dir` | Directory containing FASTQ reads |
| `refs_dir` | Directory with candidate reference genomes |
| `outdir` | Output directory |
| `threads` | Threads for multithreaded tools |
| `memory_gb` | Memory allocation (GB) |
| `genus_species` | Species identifier for chewBBACA schema |


## 7. Detailed Process Description

### 7.1 FASTP – Read preprocessing

**Input**
- Raw reads in FASTQ format

**Output**
- Trimmed reads in FASTQ format
- HTML and JSON QC reports

**Purpose**
- Remove low-quality bases
- Generate quality metrics for traceability


### 7.2 MASH – Reference selection

**Input**
- Trimmed reads in FASTQ format
- Folder with references as FASTA files

**Output**
- Distance table
- Selected reference genome code

**Purpose**
- Automatically select the most appropriate reference for each sample


### 7.3 Read Mapping (Bowtie2)

**Input**
- Trimmed reads in FASTQ format
- Selected reference genome as fasta file

**Output**
- SAM alignment file

**Purpose**
- Align IonTorrent reads to the closest reference genome


### 7.4 Alignment processing (SAMTOOLS)

**Input**
- SAM alignment file

**Output**
- Sorted and indexed BAM file

**Purpose**
- Prepare alignment for variant calling and consensus generation


### 7.5 Consensus reconstruction (samtools + bcftools + vcfutils)

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


### 7.6 FASTQ → FASTA conversion

**Input**
- Consensus as FASTQ format

**Output**
- Consensus as FASTA format

**Purpose**
- Produce final genome sequence for chewBBACA


### 7.7 chewBBACA – Allele calling

**Input**
- Consensus as FASTA format
- Folder of Species-specific schema

**Output**
- Allele profiles as TSV file
- Statistics files 

**Purpose**
- cgMLST allele calling
- Generate profiles for downstream epidemiological analysis


### 7.8. chewBBACA Schema Selection

The schema is selected automatically based on `genus_species`.

| genus_species | Schema |
|-|-|
| listeria_monocytogenes | Pasteur cgMLST |
| escherichia_coli | INNUENDO wgMLST |
| salmonella_enterica | INNUENDO cgMLST |


## 8. Output Structure

```
results/SAMPLE/
├── fastp/
├── mash/
├── ref/
├── mapping/
├── consensus/
└── chewbbaca/
```


## 9. Downloading the Project and Example Data from GitHub

### 9.1 Clone the repository

```bash
git clone https://github.com/adipi71/EOH_iontorrent.git
cd EOH_iontorrent
```


### 9.2 Example data included

```
example_refs_listeria.zip
```


### 9.3 Unzip example data

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

### 9.5 Run the pipeline

```bash
nextflow run EOH_iontorrent_pipeline_final.nf   -params-file params.json
```

## 10. Output Structure

For each sample:

```
outdir/SAMPLE/
├── fastp/        # Trimmed reads in FASTQ format and QC reports
├── mash/         # reference selection results
├── ref/          # chosen reference genome
├── mapping/      # SAM/BAM/VCF/FASTQ files
├── consensus/    # final consensus FASTA
└── chewbbaca/    # allele calling results
```


## 11. Docker Images and Tool Versions

Exact tool versions are defined in `nextflow.config`.

fastp 1.0.1
mash 2.3
bowtie2-align version 2.1.0
samtools 0.1.19-44428cd
chewBBACA 2.8.5


