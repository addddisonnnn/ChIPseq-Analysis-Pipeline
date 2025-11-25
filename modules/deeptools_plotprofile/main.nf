#!/usr/bin/env nextflow

process DEEPTOOLS_PLOTPROFILE {
    container 'ghcr.io/bf528/deeptools:latest'
    publishDir "${params.outdir}/deeptools", mode: 'copy'
    label 'process_single'

    input:
    path(matrix)

    output:
    path("profile_plot.png"), emit: plot
    path("profile_data.tab"), emit: data

    script:
    """
    plotProfile -m ${matrix} \
        -o profile_plot.png \
        --outFileNameData profile_data.tab \
        --plotTitle "Signal Profile across Gene Body" \
        --perGroup
    """

    stub:
    """
    touch profile_plot.png
    touch profile_data.tab
    """
}