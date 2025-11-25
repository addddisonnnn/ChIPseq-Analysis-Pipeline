#!/usr/bin/env nextflow

process DEEPTOOLS_COMPUTEMATRIX {
    container 'ghcr.io/bf528/deeptools:latest'
    publishDir "${params.outdir}/deeptools", mode: 'copy'
    label 'process_medium'

    input:
    path(bigwigs)
    path(bed)

    output:
    path("matrix.gz"), emit: matrix

    script:
    def bw_files = bigwigs.collect { it }.join(' ')
    """
    computeMatrix scale-regions \
        -S ${bw_files} \
        -R ${bed} \
        -b ${params.window} \
        -a ${params.window} \
        -o matrix.gz \
        -p ${task.cpus}
    """

    stub:
    """
    touch matrix.gz
    """
}