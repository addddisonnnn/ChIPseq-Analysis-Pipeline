#!/usr/bin/env nextflow

process BEDTOOLS_REMOVE {
    container 'ghcr.io/bf528/bedtools:latest'
    publishDir "${params.outdir}/peaks", mode: 'copy'
    label 'process_single'

    input:
    path(peaks)
    path(blacklist)

    output:
    path("filtered_peaks.bed"), emit: filtered_peaks

    script:
    """
    bedtools intersect -a ${peaks} -b ${blacklist} -v > filtered_peaks.bed
    """

    stub:
    """
    touch filtered_peaks.bed
    """
}