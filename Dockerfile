FROM continuumio/miniconda3:25.3.1-1
 
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
    conda install -y mash=2.3 bowtie2=2.5.4 nextflow=22.04.5; \
    conda clean -a -y; \
    which mash; which bowtie2; which nextflow; \
    nextflow -version

RUN pip install --no-cache-dir numpy biopython pandas xlrd