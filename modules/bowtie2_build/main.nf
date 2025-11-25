#!/usr/bin/env nextflow

process BOWTIE2_BUILD {
    container 'ghcr.io/bf528/bowtie2:latest'
    publishDir "${params.refdir}", mode: 'copy'
    label 'process_high'

    input:
    path(genome_fasta)

    output:
    path("bowtie2_index"), emit: index

    script:
    def genome_base = genome_fasta.baseName
    """
    mkdir -p bowtie2_index
    bowtie2-build --threads ${task.cpus} ${genome_fasta} bowtie2_index/${genome_base}
    """

    stub:
    def genome_base = genome_fasta.baseName
    """
    mkdir -p bowtie2_index
    touch bowtie2_index/${genome_base}.1.bt2
    touch bowtie2_index/${genome_base}.2.bt2
    touch bowtie2_index/${genome_base}.3.bt2
    touch bowtie2_index/${genome_base}.4.bt2
    touch bowtie2_index/${genome_base}.rev.1.bt2
    touch bowtie2_index/${genome_base}.rev.2.bt2
    """
}