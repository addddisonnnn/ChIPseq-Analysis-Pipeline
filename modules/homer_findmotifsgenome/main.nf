#!/usr/bin/env nextflow

process HOMER_FINDMOTIFSGENOME {
    container 'ghcr.io/bf528/homer_samtools:latest'
    publishDir "${params.outdir}/homer/motifs", mode: 'copy'
    label 'process_high'

    input:
    path(peaks)
    path(genome)

    output:
    path("motif_output"), emit: motif_dir

    script:
    """
    findMotifsGenome.pl ${peaks} ${genome} motif_output -size 200 -mask
    """

    stub:
    """
    mkdir -p motif_output
    touch motif_output/homerResults.html
    """
}