#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Import modules
include { FASTQC as FASTQC_RAW } from './modules/fastqc/main.nf'
include { TRIMMOMATIC } from './modules/trimmomatic/main.nf'
include { FASTQC as FASTQC_TRIMMED } from './modules/fastqc/main.nf'
include { BOWTIE2_BUILD } from './modules/bowtie2_build/main.nf'
include { BOWTIE2_ALIGN } from './modules/bowtie2_align/main.nf'
include { SAMTOOLS_SORT } from './modules/samtools_sort/main.nf'
include { SAMTOOLS_INDEX } from './modules/samtools_index/main.nf'
include { SAMTOOLS_FLAGSTAT } from './modules/samtools_flagstat/main.nf'
include { DEEPTOOLS_BAMCOVERAGE } from './modules/deeptools_bamcoverage/main.nf'
include { MULTIQC } from './modules/multiqc/main.nf'

include { DEEPTOOLS_MULTIBWSUMMARY } from './modules/deeptools_multibwsummary/main.nf'
include { DEEPTOOLS_PLOTCORRELATION } from './modules/deeptools_plotcorrelation/main.nf'
include { HOMER_MAKETAGDIR } from './modules/homer_maketagdir/main.nf'
include { HOMER_FINDPEAKS } from './modules/homer_findpeaks/main.nf'
include { HOMER_POS2BED } from './modules/homer_pos2bed/main.nf'
include { BEDTOOLS_INTERSECT } from './modules/bedtools_intersect/main.nf'
include { BEDTOOLS_REMOVE } from './modules/bedtools_remove/main.nf'
include { HOMER_ANNOTATEPEAKS } from './modules/homer_annotatepeaks/main.nf'

include { DEEPTOOLS_COMPUTEMATRIX } from './modules/deeptools_computematrix/main.nf'
include { DEEPTOOLS_PLOTPROFILE } from './modules/deeptools_plotprofile/main.nf'
include { HOMER_FINDMOTIFSGENOME } from './modules/homer_findmotifsgenome/main.nf'

workflow {
    // Create input channel from samplesheet
    // Both subsampled and full datasets are single-end
    // Parse sample name to extract condition and replicate info
    Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row -> 
            def sample_name = row.name
            // Remove _subset suffix if present
            def clean_name = sample_name.replaceAll('_subset$', '')
            def parts = clean_name.tokenize('_')
            def condition = parts[0]  // IP or INPUT
            def replicate = parts[1]  // rep1 or rep2
            
            tuple(
                clean_name,
                condition,
                replicate,
                file(row.path)
            )
        }
        .set { reads_ch }

    // ========== WEEK 1: QC, Alignment, and Coverage ==========

    // 1. FastQC on raw reads
    FASTQC_RAW(reads_ch)

    // 2. Trimmomatic - trim adapters and low quality reads
    TRIMMOMATIC(reads_ch)

    // 3. FastQC on trimmed reads
    FASTQC_TRIMMED(TRIMMOMATIC.out.trimmed_reads)

    // 4. Build Bowtie2 index (run once for reference genome)
    genome_fasta = Channel.fromPath(params.genome)
    BOWTIE2_BUILD(genome_fasta)

    // 5. Align trimmed reads to reference genome
    // Combine trimmed reads with bowtie2 index
    TRIMMOMATIC.out.trimmed_reads
        .combine(BOWTIE2_BUILD.out.index)
        .set { reads_with_index }
    
    BOWTIE2_ALIGN(reads_with_index)

    // 6. Sort BAM files
    SAMTOOLS_SORT(BOWTIE2_ALIGN.out.bam)

    // 7. Index sorted BAM files
    SAMTOOLS_INDEX(SAMTOOLS_SORT.out.sorted_bam)

    // 8. Calculate alignment statistics with flagstat
    SAMTOOLS_FLAGSTAT(SAMTOOLS_SORT.out.sorted_bam)

    // 9. Generate bigWig coverage tracks
    SAMTOOLS_SORT.out.sorted_bam
        .map { sample_id, condition, replicate, bam ->
            tuple(sample_id, bam)
        }
        .join(SAMTOOLS_INDEX.out.bam_index.map { sample_id, condition, replicate, bai ->
            tuple(sample_id, bai)
        })
        .map { sample_id, bam, bai ->
            tuple(sample_id, bam, bai)
        }
        .set { bam_with_index }

    DEEPTOOLS_BAMCOVERAGE(bam_with_index)
    
    // 10. Collect all QC outputs for MultiQC
    FASTQC_RAW.out.fastqc_zip
        .mix(FASTQC_TRIMMED.out.fastqc_zip)
        .mix(TRIMMOMATIC.out.log)
        .mix(SAMTOOLS_FLAGSTAT.out.flagstat)
        .collect()
        .set { multiqc_input }

    // 11. Run MultiQC to aggregate all QC metrics
    MULTIQC(multiqc_input)

    // ========== WEEK 2: Peak Calling and Annotation ==========
    
    // 12. Correlation analysis of bigWigs
    DEEPTOOLS_BAMCOVERAGE.out.bigwig
        .collect()
        .set { all_bigwigs }
    
    DEEPTOOLS_MULTIBWSUMMARY(all_bigwigs)
    
    // 13. Plot correlation between samples
    DEEPTOOLS_PLOTCORRELATION(DEEPTOOLS_MULTIBWSUMMARY.out.matrix)

    // 14. HOMER makeTagDirectory for all samples
    HOMER_MAKETAGDIR(SAMTOOLS_SORT.out.sorted_bam)

    // 15. Pair IP and INPUT samples for peak calling
    HOMER_MAKETAGDIR.out.tagdir
        .branch {
            ip: it[1] == 'IP'
            input: it[1] == 'INPUT'
        }
        .set { tagdirs_branched }

    tagdirs_branched.ip
        .map { sample_id, condition, replicate, tagdir -> 
            tuple(replicate, sample_id, tagdir)
        }
        .set { ip_tagdirs }

    tagdirs_branched.input
        .map { sample_id, condition, replicate, tagdir -> 
            tuple(replicate, tagdir)
        }
        .set { input_tagdirs }

    ip_tagdirs
        .join(input_tagdirs)
        .map { replicate, ip_sample, ip_tagdir, input_tagdir ->
            tuple(ip_sample, replicate, ip_tagdir, input_tagdir)
        }
        .set { paired_tagdirs }

    // 16. HOMER findPeaks - call peaks for each replicate pair
    HOMER_FINDPEAKS(paired_tagdirs)

    // 17. Convert peaks to BED format
    HOMER_POS2BED(HOMER_FINDPEAKS.out.peaks)

    // 18. Get reproducible peaks via bedtools intersect
    HOMER_POS2BED.out.bed
        .map { replicate, bed -> bed }
        .collect()
        .set { all_peak_beds }

    all_peak_beds
        .map { beds -> tuple(beds[0], beds[1]) }
        .set { peak_pair }

    BEDTOOLS_INTERSECT(peak_pair.map { it[0] }, peak_pair.map { it[1] })

    // 19. Remove blacklist regions from peaks
    blacklist = Channel.fromPath(params.blacklist)
    BEDTOOLS_REMOVE(BEDTOOLS_INTERSECT.out.reproducible_peaks, blacklist)

    // 20. Annotate filtered peaks to nearest genomic features
    genome_fasta_annot = Channel.fromPath(params.genome)
    gtf_annot = Channel.fromPath(params.gtf)
    HOMER_ANNOTATEPEAKS(BEDTOOLS_REMOVE.out.filtered_peaks, genome_fasta_annot, gtf_annot)

    // ========== WEEK 3: Signal Profiling and Motif Analysis ==========
    
    // 21. Filter bigwigs to only IP samples for gene body analysis
    DEEPTOOLS_BAMCOVERAGE.out.bigwig
        .filter { bw_file ->
            bw_file.name.contains('IP_') && !bw_file.name.contains('INPUT')
        }
        .collect()
        .set { ip_bigwigs }

    // 22. computeMatrix - calculate signal across gene regions
    ucsc_genes_bed = Channel.fromPath(params.ucsc_genes)
    DEEPTOOLS_COMPUTEMATRIX(ip_bigwigs, ucsc_genes_bed)

    // 23. plotProfile - visualize the signal profile across genes
    DEEPTOOLS_PLOTPROFILE(DEEPTOOLS_COMPUTEMATRIX.out.matrix)

    // 24. Motif enrichment analysis on filtered peaks
    genome_fasta_motif = Channel.fromPath(params.genome)
    HOMER_FINDMOTIFSGENOME(BEDTOOLS_REMOVE.out.filtered_peaks, genome_fasta_motif)
}

workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'SUCCESS' : 'FAILED' }"
    println "Execution duration: $workflow.duration"
}
