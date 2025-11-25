#!/usr/bin/env nextflow

process DEEPTOOLS_BAMCOVERAGE {
    container 'ghcr.io/bf528/deeptools:latest'
    publishDir "${params.outdir}/bigwig", mode: 'copy'
    label 'process_medium'

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    path("${sample_id}.bw"), emit: bigwig

    script:
    """
    bamCoverage -b ${bam} -o ${sample_id}.bw -p ${task.cpus}
    """

    stub:
    """
    touch ${sample_id}.bw
    """
}