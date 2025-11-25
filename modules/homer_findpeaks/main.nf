#!/usr/bin/env nextflow

process HOMER_FINDPEAKS {
    container 'ghcr.io/bf528/homer_samtools:latest'
    publishDir "${params.outdir}/homer/peaks", mode: 'copy'
    label 'process_medium'

    input:
    tuple val(ip_sample), val(replicate), path(ip_tagdir), path(input_tagdir)

    output:
    tuple val(replicate), path("${ip_sample}_peaks.txt"), emit: peaks

    script:
    """
    findPeaks ${ip_tagdir} -style factor -o ${ip_sample}_peaks.txt -i ${input_tagdir}
    """

    stub:
    """
    touch ${ip_sample}_peaks.txt
    """
}