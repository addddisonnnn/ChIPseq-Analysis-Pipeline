#!/usr/bin/env nextflow

process DEEPTOOLS_MULTIBWSUMMARY {
    container 'ghcr.io/bf528/deeptools:latest'
    publishDir "${params.outdir}/deeptools", mode: 'copy'
    label 'process_medium'

    input:
    path(bigwigs)

    output:
    path("multibw_summary.npz"), emit: matrix

    script:
    def bw_files = bigwigs.collect { it }.join(' ')
    """
    multiBigwigSummary bins -b ${bw_files} -o multibw_summary.npz -p ${task.cpus}
    """

    stub:
    """
    touch multibw_summary.npz
    """
}