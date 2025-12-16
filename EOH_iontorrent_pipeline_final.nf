nextflow.enable.dsl = 2

// ================== PARAMS ==================

params.threads   = (params.threads ?: 8)    as int
params.memory_gb = (params.memory_gb ?: 16) as int

process FASTP {
  tag { basename }
  publishDir "${params.outdir}", mode: 'copy', saveAs: { fn -> "${basename}/fastp/${fn}" }
  cpus { params.threads }
  memory { "${params.memory_gb} GB" }

  input:
    tuple val(basename), path(reads)

  output:
    tuple val(basename), path("*.trim.fq.gz"), emit: trimmed
    path "*.fastp.json"
    path "*.fastp.html"

  script:
  def bn = reads.getBaseName()
  def sample = bn.replaceAll(/(\.fastq|\.fq)(\.gz)?$/, '')
  """
  fastp -i ${reads} -o ${sample}.trim.fq.gz \
        --thread ${task.cpus} \
        --json ${sample}.fastp.json \
        --html ${sample}.fastp.html
  """
}

process MASH_SKETCH_REFS {
  tag "refs"
  publishDir "${params.outdir}/_global/mash", mode: 'copy'
  cpus 1

  input:
    path ref_files

  output:
    path "refs.msh", emit: refs_msh

  script:
  """
  printf "%s\n" ${ref_files} > refs.list
  [ -s refs.list ] || { echo "No FASTA files staged" >&2; exit 1; }
  mash sketch -o refs \$(cat refs.list)
  """
}

process MASH_SKETCH_READS {
  tag { basename }
  publishDir "${params.outdir}", mode: 'copy', saveAs: { fn -> "${basename}/mash/${fn}" }
  cpus 1

  input:
    tuple val(basename), path(trimmed)

  output:
    tuple val(basename), path("*.reads.msh"), emit: reads_msh

  script:
  def sample = trimmed.getBaseName().replaceAll(/\.trim$/, '')
  """
  mash sketch -o ${sample}.reads ${trimmed}
  """
}

process MASH_DIST {
  tag { basename }
  publishDir "${params.outdir}", mode: 'copy', saveAs: { fn -> "${basename}/mash/${fn}" }
  cpus 1

  input:
    tuple val(basename), path(reads_msh)
    path refs_msh

  output:
    tuple val(basename), path("mash_dist.tsv"), emit: dist

  script:
  """
  mash dist ${refs_msh} ${reads_msh} > mash_dist.tsv
  """
}

process PICK_REF {
  tag { basename }
  publishDir "${params.outdir}", mode: 'copy', saveAs: { fn -> "${basename}/mash/${fn}" }
  cpus 1

  input:
    tuple val(basename), path(dist)

  output:
    tuple val(basename), path("chosen_ref.txt"), emit: chosen_ref_txt
    tuple val(basename), path("best_dist.txt"),  emit: best_dist_txt

  script:
  """
  awk 'BEGIN{FS="\\t"; min=1e9}
       {
         if (NF>=3 && \$3+0<min) {min=\$3; ref=\$1}
       }
       END{
         if (ref=="") exit 1;
         print ref > "chosen_ref.txt";
         print min > "best_dist.txt"
       }' ${dist}
  """
}

process PREP_REF {
  tag { basename }
  publishDir "${params.outdir}", mode: 'copy', saveAs: { fn -> "${basename}/ref/${fn}" }
  cpus 1

  input:
    tuple val(basename),  path(chosen_ref_txt)
    tuple val(basename2), path(best_dist_txt)

  output:
    tuple val(basename), path("ref.fasta"),     emit: ref_fasta
    tuple val(basename), path("best_dist.txt")

  when:
    basename == basename2

  script:
  """
  set -euo pipefail

  REF=\$(cat ${chosen_ref_txt})
  D=\$(cat ${best_dist_txt})
  echo "\$D" > best_dist.txt

  BASE=\$(basename "\$REF")
  REAL=\$(find "${params.refs_dir}" -type f -name "\$BASE" | head -n 1 || true)
  [ -n "\$REAL" ] || { echo "Cannot resolve ref '\$REF' in ${params.refs_dir}" >&2; exit 1; }

  cp "\$REAL" ref.fasta
  """
}

process BOWTIE2_ALIGN_SE {
  tag { basename }
  publishDir "${params.outdir}", mode: 'copy', saveAs: { fn -> "${basename}/mapping/${fn}" }
  label 'mapping'
  cpus 8

  input:
    tuple val(basename), path(fq), path(ref_fa)

  output:
    tuple val(basename), path("*.sam"), path(ref_fa)

  script:
  def referencePath = ref_fa
  def reference     = ref_fa.getBaseName()
  def base          = fq.getBaseName().replaceAll(/(\.fastq|\.fq)(\.gz)?$/, '').replaceFirst(/\.trim$/, '')
  def base_ref      = "${base}_bowtie_${reference}"
  """
  set -euo pipefail
  bowtie2-build ${referencePath} ${reference}
  bowtie2 -p ${task.cpus} --very-fast -x ${reference} -U ${fq} -S ${base_ref}.sam 2>> ${base_ref}.log
  """
}

process SAMTOOLS_VCFUTILS_FASTQ {
  tag { basename }
  publishDir "${params.outdir}", mode: 'copy', saveAs: { fn -> "${basename}/mapping/${fn}" }
  label 'mapping'
  containerOptions '-u 0:0'
  cpus 8

  input:
  tuple val(basename), path(sam), path(referencePath)

  output:
  tuple val(basename), path("*.fq"), emit: fq

  script:
  def base_ref  = sam.getBaseName()

  // Iontorrent optimizations
  def maxDepth  = (params.max_depth  ?: 150) as int
  def minMapQ   = (params.min_mapq   ?: 10)  as int
  def minBaseQ  = (params.min_baseq  ?: 10)  as int
  def minDP     = (params.min_dp     ?: 8)   as int

  """
  set -euo pipefail

  samtools view -bS -o ${base_ref}.bam ${sam} 2>> ${base_ref}.log
  samtools sort ${base_ref}.bam ${base_ref}_sorted 2>> ${base_ref}.log
  samtools index ${base_ref}_sorted.bam 2>> ${base_ref}.log

  samtools mpileup \
    -d ${maxDepth} \
    -q ${minMapQ} -Q ${minBaseQ} \
    -uf ${referencePath} \
    ${base_ref}_sorted.bam \
    > ${base_ref}.bcf 2>> ${base_ref}.log

  bcftools view -cg ${base_ref}.bcf > ${base_ref}.var.flt.vcf 2>> ${base_ref}.log

  vcfutils.pl vcf2fq \
    -d ${minDP} \
    -D ${maxDepth} \
    ${base_ref}.var.flt.vcf \
    > ${base_ref}.fq 2>> ${base_ref}.log
  """
}

process FASTQ_TO_FASTA_SEQIO {
  tag { basename }
  publishDir "${params.outdir}", mode: 'copy', saveAs: { fn -> "${basename}/consensus/${fn}" }
  label 'mapping'
  cpus 1

  input:
  tuple val(basename), path(fq)

  output:
  tuple val(basename), path("*.consensus.fasta"), emit: consensus_fa

  script:
  def base = fq.getBaseName()
  """
  set -euo pipefail
  python3 - << 'PY'
  from Bio import SeqIO
  SeqIO.convert("${fq}", "fastq", "${base}.consensus.fasta", "fasta")
  PY
  """
}

process CHEWBBACA {
  tag { basename }
  publishDir "${params.outdir}", mode: 'copy', saveAs: { fn -> "${basename}/chewbbaca/${fn}" }
  cpus { Math.max(8, params.threads as int) }
  memory { "8 GB" }

  input:
    tuple val(basename), path(fa)
    val schema_path

  output:
    path "results/*/results_alleles.tsv",    emit: alleles_tsv
    path "results/*/results_statistics.tsv", emit: stats_tsv
    path "*_new_alleles.txt", optional: true, emit: new_alleles
    path "schema", emit: schema_dir

  script:
  def base = fa.getBaseName().replaceAll(/[^A-Za-z0-9._-]/, "_")
  """
  set -euo pipefail

  unzip ${schema_path} -d schema > /dev/null
  mkdir -p input
  cp ${fa} input/contigs.fa

  chewBBACA.py AlleleCall \
    -i input -g schema -o results \
    --cpu ${task.cpus} --force-continue --verbose

  d=\$(find results -maxdepth 1 -type d -name 'results_*' | head -n1)
  grep "\$(basename ${fa} | sed 's/_/-/g')" schema/*.fasta -A1 -h | grep -v "\\-\\-" > ${base}_new_alleles.txt || true
  """
}


// ================== WORKFLOW ==================

workflow {

  if( !params.reads_dir ) error "Missing param: reads_dir"
  if( !params.outdir )    error "Missing param: outdir"
  if( !params.refs_dir  ) error "Missing param: refs_dir"

  def rpat = [
    "${params.reads_dir}/*.fastq.gz",
    "${params.reads_dir}/*.fq.gz",
    "${params.reads_dir}/*.fastq",
    "${params.reads_dir}/*.fq"
  ]
  def rchs = rpat.collect { Channel.fromPath(it, followLinks:true, checkIfExists:false) }
  def ch_reads = rchs.drop(1).inject(rchs[0]) { a,b -> a.mix(b) }.unique().ifEmpty {
    error "No reads found under ${params.reads_dir}"
  }.map { f ->
    def basename = f.getBaseName().replaceAll(/(\.fastq|\.fq)(\.gz)?$/, '')
    [ basename, f ]
  }

  def pats = [
    "${params.refs_dir}/*.fa",
    "${params.refs_dir}/*.fasta",
    "${params.refs_dir}/*.fna",
    "${params.refs_dir}/*.fa.gz",
    "${params.refs_dir}/*.fasta.gz",
    "${params.refs_dir}/*.fna.gz"
  ]

  def chans = pats.collect { pat -> Channel.fromPath(pat, followLinks:true, checkIfExists:false) }
  def ch_ref_files = chans
    .inject(Channel.empty()) { acc, ch -> acc.mix(ch) }
    .unique()
    .ifEmpty {
      error "No FASTA found directly under ${params.refs_dir}\nExpected one of:\n  - " + pats.join("\n  - ")
    }

  def ch_ref_files_list = ch_ref_files.collect()
  def ch_trimmed      = FASTP(ch_reads).trimmed

  def ch_refs_msh_val = MASH_SKETCH_REFS(ch_ref_files_list).refs_msh.first()
  def ch_reads_msh    = MASH_SKETCH_READS(ch_trimmed).reads_msh
  def ch_dist         = MASH_DIST(ch_reads_msh, ch_refs_msh_val).dist
  def ch_picked       = PICK_REF(ch_dist)
  def ch_prep         = PREP_REF(ch_picked.chosen_ref_txt, ch_picked.best_dist_txt)

  def ch_sam          = BOWTIE2_ALIGN_SE( ch_trimmed.join(ch_prep.ref_fasta) )
  def ch_samtools = SAMTOOLS_VCFUTILS_FASTQ(ch_sam)
  def ch_consensus = FASTQ_TO_FASTA_SEQIO( ch_samtools ).consensus_fa
  
  if( !params.genus_species )
    error "chewBBACA: serve params.genus_species"

  def schema_key = params.schema ?: (
      params.genus_species?.toLowerCase() == 'listeria_monocytogenes' ? 'l_mono_chewie_1748_220623' :
      params.genus_species?.toLowerCase() == 'escherichia_coli'       ? 'e_coli_chewie_2360_210531' :
      params.genus_species?.toLowerCase() == 'salmonella_enterica'     ? 's_enterica_chewie_3255_210531' :
  null
  )
  if( !schema_key ) error "chewBBACA: specifica params.schema"

  def schema_path = (
      schema_key == 'l_mono_chewie_1748_220623'     ? '/schemas/Listeria_monocytogenes_Pasteur_cgMLST_2022-06-23T18_03_54.613576.zip' :
      schema_key == 'e_coli_chewie_2360_210531'     ? '/schemas/Escherichia_coli_INNUENDO_wgMLST_2021-05-31T14_24_05.304225.zip' :
      schema_key == 's_enterica_chewie_3255_210531' ? '/schemas/Salmonella_enterica_INNUENDO_cgMLST_2021-05-31T20_28_21.350919.zip' :
  null
  )
  if( !schema_path ) error "chewBBACA: schema non riconosciuto"


  def ch_chew = CHEWBBACA(ch_consensus, schema_path)
  def ch_chew_profiles = ch_chew.alleles_tsv.map { p -> tuple(p, file(schema_path)) }
}
