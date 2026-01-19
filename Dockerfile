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