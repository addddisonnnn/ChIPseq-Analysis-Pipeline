#!/usr/bin/env nextflow

process HOMER_MAKETAGDIR {
    container 'ghcr.io/bf528/homer_samtools:latest'
    publishDir "${params.outdir}/homer/tagdirs", mode: 'copy'
    label 'process_medium'

    input:
    tuple val(sample_id), val(condition), val(replicate), path(bam)

    output:
    tuple val(sample_id), val(condition), val(replicate), path("${sample_id}_tagdir"), emit: tagdir

    script:
    """
    makeTagDirectory ${sample_id}_tagdir ${bam}
    """

    stub:
    """
    mkdir -p ${sample_id}_tagdir
    touch ${sample_id}_tagdir/tagInfo.txt
    """
}