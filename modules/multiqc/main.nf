#!/usr/bin/env nextflow

process MULTIQC {
    container 'ghcr.io/bf528/multiqc:latest'
    publishDir "${params.outdir}/multiqc", mode: 'copy'
    label 'process_single'

    input:
    path(qc_files)

    output:
    path("multiqc_report.html"), emit: report
    path("multiqc_data"), emit: data

    script:
    """
    multiqc .
    """

    stub:
    """
    mkdir -p multiqc_data
    touch multiqc_report.html
    touch multiqc_data/multiqc_data.json
    """
}