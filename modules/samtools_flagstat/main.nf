#!/usr/bin/env nextflow

process SAMTOOLS_FLAGSTAT {
    container 'ghcr.io/bf528/samtools:latest'
    publishDir "${params.outdir}/flagstat", mode: 'copy'
    label 'process_single'

    input:
    tuple val(sample_id), val(condition), val(replicate), path(sorted_bam)

    output:
    path("${sample_id}_flagstat.txt"), emit: flagstat

    script:
    """
    samtools flagstat ${sorted_bam} > ${sample_id}_flagstat.txt
    """

    stub:
    """
    touch ${sample_id}_flagstat.txt
    """
}