process BOWTIE2_ALIGN {
    container 'ghcr.io/bf528/bowtie2:latest'
    publishDir "${params.outdir}/alignments", mode: 'copy'
    label 'process_high'

    input:
    tuple val(sample_id), val(condition), val(replicate), path(reads), path(index_dir)

    output:
    tuple val(sample_id), val(condition), val(replicate), path("${sample_id}.bam"), emit: bam

    script:
    def genome_base = file(params.genome).baseName
    def bowtie_threads = task.cpus - 5  // Reserve threads for overhead
    """
    bowtie2 -p ${bowtie_threads} \
        -x ${index_dir}/${genome_base} \
        -U ${reads} \
        | samtools view -bS - > ${sample_id}.bam
    """

    stub:
    """
    touch ${sample_id}.bam
    """
}