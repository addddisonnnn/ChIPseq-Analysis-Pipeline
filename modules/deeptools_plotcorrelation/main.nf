#!/usr/bin/env nextflow

process DEEPTOOLS_PLOTCORRELATION {
    container 'ghcr.io/bf528/deeptools:latest'
    publishDir "${params.outdir}/deeptools", mode: 'copy'
    label 'process_single'

    input:
    path(matrix)

    output:
    path("correlation_heatmap.png"), emit: plot
    path("correlation_matrix.tab"), emit: matrix_file

    script:
    """
    plotCorrelation -in ${matrix} \
        --corMethod spearman \
        --skipZeros \
        --plotTitle "Spearman Correlation of Read Counts" \
        --whatToPlot heatmap \
        --colorMap RdYlBu \
        --plotNumbers \
        -o correlation_heatmap.png \
        --outFileCorMatrix correlation_matrix.tab
    """

    stub:
    """
    touch correlation_heatmap.png
    touch correlation_matrix.tab
    """
}