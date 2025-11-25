#!/usr/bin/env nextflow

process HOMER_POS2BED {
    container 'ghcr.io/bf528/homer_samtools:latest'
    publishDir "${params.outdir}/homer/peaks_bed", mode: 'copy'
    label 'process_single'

    input:
    tuple val(replicate), path(peaks_txt)

    output:
    tuple val(replicate), path("${peaks_txt.baseName}.bed"), emit: bed

    script:
    """
    pos2bed.pl ${peaks_txt} > ${peaks_txt.baseName}.bed
    """

    stub:
    """
    touch ${peaks_txt.baseName}.bed
    """
}