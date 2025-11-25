# ChIP-seq Analysis: RUNX1 Binding in MCF-7 Breast Cancer Cells
[![Nextflow](https://img.shields.io/badge/nextflow-%23EAD69C.svg?style=for-the-badge&logo=nextflow&logoColor=000000)](https://www.nextflow.io/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://img.shields.io/badge/DOI-10.1016%2Fj.bbagrm.2016.100-00A3E0)](https://doi.org/10.1016/j.bbagrm.2016.100)

A comprehensive re-analysis of ChIP-seq data from Barutcu et al. (2016) investigating RUNX1's role as a transcriptional regulator in breast cancer cells, implemented as a reproducible Nextflow pipeline.

## Project Overview
This project implements a reproducible Nextflow pipeline for ChIP-seq analysis, focusing on RUNX1 transcription factor binding in MCF-7 breast cancer cells. The analysis includes quality control, peak calling, motif discovery, and integration with gene expression data to identify direct RUNX1 targets.

**Key Findings:**
- Identified 6,462 high-confidence RUNX1 binding sites
- Discovered RUNX1's role as a master regulator coordinating multiple cancer hallmarks
- Validated direct transcriptional regulation of key breast cancer genes
- Revealed extensive involvement in cell cycle, DNA damage response, and metabolic pathways
## Prerequisites
- Nextflow >= 22.10.0
- Singularity/Apptainer or Docker
- Python 3.9 or higher

## Installation
```bash
git clone https://github.com/yourusername/chipseq-runx1-analysis.git
cd chipseq-runx1-analysis
```

## Running the Pipeline
```bash
nextflow run main.nf -profile singularity --samplesheet samplesheet.csv
```
## Project Strucutre

``` text
chipseq-runx1-analysis/
├── main.nf                          # Main Nextflow pipeline
├── modules/                         # Custom Nextflow modules
│   ├── fastqc/
│   ├── trimmomatic/
│   ├── bowtie2_align/
│   ├── homer_findpeaks/
│   └── ...
├── results/                         # Analysis outputs
│   ├── MultiQC_Report.pdf          # Comprehensive QC report
│   ├── profile_plot.png            # Signal coverage across genes
│   ├── known_homer_motif.png       # Motif enrichment results
│   ├── figure_2F_accurate.png      # Integration with expression data
│   └── ...
├── config/                         # Configuration files
└── docs/                          # Documentation
```
## Analysis Workflow
1. Quality Control & Preprocessing
  - FastQC: Quality assessment of raw and trimmed reads
  - Trimmomatic: Adapter trimming and quality filtering
  - MultiQC: Aggregated QC report generation
2. Genome Alignment
  - Bowtie2: Alignment to hg38 reference genome
  - SAMtools: BAM processing and statistics
3. Peak Calling & Analysis
  - HOMER: Peak calling with factor-style parameters
  - BEDTools: Peak intersection and blacklist filtering
  - deepTools: Correlation analysis and signal visualization
4. Motif & Functional Analysis
  - HOMER: De novo motif discovery and enrichment
  - GREAT/Enrichr: Gene ontology and pathway analysis
  - IGV: Genomic visualization of specific loci

## Key Results
### Quality Control Metrics
| Sample | Total Reads (M) | Mapping Rate | Duplication Rate | GC Content |
|--------|-----------------|--------------|------------------|------------|
| INPUT_rep1 | ~30 | 89.1% | 1-4% | 43-47% |
| INPUT_rep2 | ~10.7 | 74.3% | 1-4% | 43-47% |
| IP_rep1 | ~30 | 89.1% | 12-13% | 43-47% |
| IP_rep2 | ~29 | 74.3% | 12-13% | 43-47% |
### RUNX1 Binding Characteristics
- Promoter enrichment: Strong binding at transcription start sites
- Motif validation: Canonical RUNX motif highly enriched (p < 1e-50)
- Target genes: Direct regulation of 14.7-18.7% of differentially expressed genes
### Functional Enrichment
**Top Pathways Regulated by RUNX1:**
- E2F Targets & Cell Cycle Control
- DNA Damage Response & Repair
- Metabolic Reprogramming (Glycolysis, TCA Cycle)
- Estrogen Response Pathways
- mTORC1 Signaling
### Visualizations
The pipeline generates comprehensive visualizations including:
- Signal coverage profiles across gene bodies
- Correlation heatmaps between replicates
- Motif enrichment plots
- Genome browser tracks for specific loci
- Pathway enrichment bar graphs
### Configuration
**Profiles**
- singularity: Uses Singularity containers for reproducibility
- docker: Uses Docker containers
- local: Local execution with Conda/Mamba
**Parameters**
Key customizable parameters in nextflow.config:
- genome: Reference genome path
- blacklist: ENCODE blacklist regions
- gtf: Gene annotation file
- ucsc_genes: Gene coordinates for profiling

## Citation
Barutcu, A. R., Hong, D., Lajoie, B. R., McCord, R. P., van Wijnen, A. J., Lian, J. B., Stein, J. L., Dekker, J., Imbalzano, A. N., & Stein, G. S. (2016). RUNX1 contributes to higher-order chromatin organization and gene regulation in breast cancer cells. Biochimica et biophysica acta, 1859(11), 1389–1397. https://doi.org/10.1016/j.bbagrm.2016.08.003
```bibtex
@article{barutcu2016,
  title={RUNX1 contributes to higher-order chromatin organization and gene regulation in breast cancer cells},
  author={Barutcu, A Rasim and Hong, Deli and Lajoie, Bryan R and McCord, Rachel P and van Wijnen, Andre J and Lian, Jane B and Stein, Janet L and Dekker, Job and Imbalzano, Anthony N and Stein, Gary S},
  journal={Biochimica et Biophysica Acta (BBA)-Gene Regulatory Mechanisms},
  volume={1859},
  number={12},
  pages={1389--1397},
  year={2016},
  publisher={Elsevier}
}
```
