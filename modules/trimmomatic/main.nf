#!/usr/bin/env nextflow

process TRIMMOMATIC {
    container 'ghcr.io/bf528/trimmomatic:latest'
    publishDir "${params.outdir}/trimmomatic", mode: 'copy'
    label 'process_medium'

    input:
    tuple val(sample_id), val(condition), val(replicate), path(reads)

    output:
    tuple val(sample_id), val(condition), val(replicate), path("${sample_id}_trimmed.fastq.gz"), emit: trimmed_reads
    path("${sample_id}_trimmomatic.log"), emit: log

    script:
    """
    trimmomatic SE -threads ${task.cpus} \
        ${reads} \
        ${sample_id}_trimmed.fastq.gz \
        ILLUMINACLIP:${params.adapter_fa}:2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 \
        2> ${sample_id}_trimmomatic.log
    """

    stub:
    """
    touch ${sample_id}_trimmed.fastq.gz
    touch ${sample_id}_trimmomatic.log
    """
}