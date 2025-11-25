#!/usr/bin/env nextflow

process SAMTOOLS_SORT {
    container 'ghcr.io/bf528/samtools:latest'
    publishDir "${params.outdir}/sorted_bams", mode: 'copy'
    label 'process_medium'

    input:
    tuple val(sample_id), val(condition), val(replicate), path(bam)

    output:
    tuple val(sample_id), val(condition), val(replicate), path("${sample_id}.sorted.bam"), emit: sorted_bam

    script:
    """
    samtools sort -@ ${task.cpus} -o ${sample_id}.sorted.bam ${bam}
    """

    stub:
    """
    touch ${sample_id}.sorted.bam
    """
}