#!/usr/bin/env nextflow

process HOMER_ANNOTATEPEAKS {
    container 'ghcr.io/bf528/homer_samtools:latest'
    publishDir "${params.outdir}/homer/annotations", mode: 'copy'
    label 'process_medium'

    input:
    path(peaks)
    path(genome)
    path(gtf)

    output:
    path("annotated_peaks.txt"), emit: annotations

    script:
    """
    annotatePeaks.pl ${peaks} ${genome} -gtf ${gtf} > annotated_peaks.txt
    """

    stub:
    """
    touch annotated_peaks.txt
    """
}