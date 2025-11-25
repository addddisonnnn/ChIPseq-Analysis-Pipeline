#!/usr/bin/env nextflow

process SAMTOOLS_INDEX {
    container 'ghcr.io/bf528/samtools:latest'
    publishDir "${params.outdir}/sorted_bams", mode: 'copy'
    label 'process_low'

    input:
    tuple val(sample_id), val(condition), val(replicate), path(sorted_bam)

    output:
    tuple val(sample_id), val(condition), val(replicate), path("${sorted_bam}.bai"), emit: bam_index

    script:
    """
    samtools index -@ ${task.cpus} ${sorted_bam}
    """

    stub:
    """
    touch ${sorted_bam}.bai
    """
}