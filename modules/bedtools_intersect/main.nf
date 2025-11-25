#!/usr/bin/env nextflow

process BEDTOOLS_INTERSECT {
    container 'ghcr.io/bf528/bedtools:latest'
    publishDir "${params.outdir}/peaks", mode: 'copy'
    label 'process_single'

    input:
    path(bed1)
    path(bed2)

    output:
    path("reproducible_peaks.bed"), emit: reproducible_peaks

    script:
    """
    bedtools intersect -a ${bed1} -b ${bed2} -wa > reproducible_peaks.bed
    """

    stub:
    """
    touch reproducible_peaks.bed
    """
}